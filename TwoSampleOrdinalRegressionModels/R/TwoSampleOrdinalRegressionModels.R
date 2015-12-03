## This file contains routines in R that perform maximum likelihood 
## ordinal regression on two sample ROC data (a la Dorfman and everyone else)
## In addition to calculating the empirical AUC value and its uncertainty,
## this software can fit the following models to the data:

## PowerLaw: a power-law or bi-exponential model, e.g. Egan 1975
## BiGamma: a two-parameter bi-gamma model a la Dorfman et al. 1997
## BiNormal: the standard two-parameter bi-normal model a la Dorfman & Alf 1969
## BiNormalEqVar: an equal variance bi-normal model (1 parameter)
## CBM: the two-parameter "Contaminated Binormal Model", Dorfman et al. 2000
## PropROC or ChiSq: the two parameter bi-chi-squared or "PropROC" model, a la Metz & Pan 1999 and Hillis 2015
## BiLogistic:  A two-parameter logistic model
## BiLogisticEqVar:  A one-parameter logistic model
## DualBeta: The two parameter beta model
## ... More to come?

## Examples:
##   source("TwoSampleOrdinalRegressionModels.R")  # Read this file.
##   x1=round(log(rexp(100,1/5)))    # Make some 
##   x0=round(log(rexp(100)))        # fake data
##   fits=TwoSampleOrdinalRegression(x1,x0,models=c("PowerLaw","BiGamma","BiNormal","BiNormalEqVar","CBM", "PropROC", "BiLogistic","BiLogistic1"))
##   print(fits)   # print results of the fits
##   plot(fits)    # too much of a good thing...
## 
##   t1=1:8      # Make some fake ordinal categorical counts.
##   t0=8:1
##   fits=TwoSampleOrdinalRegression(counts1=t1,counts0=t0,
##                              models=c("BiNormalEqVar", "BiLogisticEqVar"))
##   print(fits)   # Compare the fits
##   plot(fits)

## Printing the fits yields a table that gives the name of the model,
## the the number of parameters in the model, the log-likelihood, the
## AUC and its estimated standard error, and the Akaike and Bayesian
## information criteria.
## *Note that the "number of parameters" does not include the categorical
## boundaries which are also fit in the model.  The number of these
## categories (and therefore parameters) may depend on the data.
## Likewise, the log-likelihood, the AIC, and the BIC depend on the number
## of parameters, and therefore may not be meaningful when compared with
## a different fitting program.  However, the number of parameters,
## log-likelihood, AIC, and BIC are comparable among the models
## fit here.
## The fits should almost always be plotted to make certain that they
## converged to something that is reasonable.

## Several things to note:
## 0. This code is really slow, particularly the bi-chi-squared/PropROC model.
## 1. Category boundaries are done in probability space.
## 2. Likelihood maximization is performed using the vanilla R
##    optimization routine "optim()", and some generic initial guess values.
##    This means that the maximization sometimes fails to converge.  
##    We work around this by performing the maximization many times with many
##    different initial parameter values.  Statistically at least
##    one of them will converge.  This is very inefficient, 
##    rather dumb, but robust.
## 3. This code was never intended for external consumption, and
##    so therefore it is not well organized and you should not
##    have high expectations.


logfactorial=function (n) {
  ## An excellent approximation to the log of factorial
  ## snagged from wikipedia.
  ret=n
  wh=n<20;   # Do not approximate small values
  ret[wh]=log(factorial(n[wh]))
  nn=n[!wh]  # Approximate large values.
  ret[!wh]=( nn*log(nn) -nn + log( nn*(1+4*nn*(1+2*nn)))/6+ log(pi)/2 )
  return(ret)
}

ROCxy<-function(neg,pos) {  
  ## Input neg= a vector of scores from signal absent images
  ##       pos = a vector from sig present images
  ## Output:  Empirical Points for an ROC curve.  Handles ties.
  lvls<-sort(unique(c(neg,pos))) #,decreasing=T)
  t0=table(factor(neg,levels=lvls))
  t1=table(factor(pos,levels=lvls))
  return(ROCxyCounts(list(t1,t0)))
}

ROCxyCounts=function(counts) {
  ## Input a list of two vectors.  One vector is the number of signal
  ## present observations in each category, and the other is the number
  ## of signal absent observations.
  ## Output:  Empirical Points for an ROC curve.  Handles ties.
  c0=cumsum(counts[[2]])
  c1=cumsum(counts[[1]])
  return(list(x=1-c(0,c0)/c0[length(c0)],y=1-c(0,c1)/c1[length(c1)]))
}


AUCcat=function(counts) {
  # Calculate empirical AUC estimate from count data.
  c1=counts[[1]];  c0=counts[[2]] 
  n=length(c0)
  cu0=cumsum(c0)
  a=(sum(.5*c1*c0) + sum(c1[2:n] * cu0[1:(n-1)]))/cu0[n]/sum(c1)
  return(unname(a))
}

AUC = function (p1,n1) { # Calculate empirical AUC from score data
  ## Inputs n1= a vector of scores from signal absent images
  ##        p1 = a vector from sig present images
  nop=length(p1)
  empa=(sum(rank(c(p1,n1))[0:nop])- nop*(nop+1)/2)/nop/length(n1)
  empa   # return empirical auc
}

MWWStatVar=function(y,x,y0=NA,x0=NA){
  # This function returns the Mann-Whitney statistic and its unbiased, 
  # non-negative variance estimate. Arguments y and x are the two samples.
  # If y0 and x0 are provided, the difference of two Mann-Whitney
  # statistics and the variance estimate of that difference is output.
  k <- function(x,y) (sign(y-x)+1)/2;  # The indicator kernel function
  m <- length(y); n <- length(x);  X <- outer(x,y,k);
  if (length(y)==length(y0) && length(x)==length(x0)) X<- X-outer(x0,y0,k)
  MSA <- m*var(rowMeans(X)) ;  MSB <- n*var(colMeans(X)) ;  
  MST <- var(as.vector(X))
  ev <- ((n*m-1)*MST-(n-1)*MSA-(m-1)*MSB)/((n-1)*(m-1))
  return(c(mean(X), (MSA+MSB-ev)/(m*n)))
}


tablit=function(ll) {
  ## Takes a list of two vectors (one vector of signal present scores 
  ## and one vector of signal absent scores) and 
  ## makes a table of frequencies of each decision 
  ## category suitable for sending to the likelihood function
  lvls= sort(unique(c(unlist(ll))))
  lapply(ll, function(z) table(factor(z,levels=lvls)))
}

##  Frank:  Combine adjacent bins for which there are 0 entries in the
##  the other bins, i.e. 1st concatenate runs, then concatentate other stuff.

tablit2=function(ll,binlim=12,minbinlim=2) {  
   ## Take tables of category frequencies, and rebin
   ## the data into no more than binlim bins
   ## This can make likelihood maximization faster and more reliable.
   ## minbin is the minumum number of counts each bin should have.
   ## binlim is the maximum number of bins allowed.
   ## In general for ML fitting, having more than 10 bins or only
   ## one count per bin is not helpful.

   nbins=length(ll[[1]])
   binelim=function(ll, elim) {
      ll=lapply(ll, function(il)
             c( il[0:(elim-1)], il[elim]+il[elim+1], 
           ifelse(nbins>elim+1, 
             list(il[(elim+2):nbins]),list(c()))[[1]]))       
   }
   i=1
   while(i< nbins) {
      if ((ll[[1]][i]==0 && ll[[1]][i+1]==0) ||
          (ll[[2]][i]==0 && ll[[2]][i+1]==0)) {
          ll=binelim(ll,i)
          nbins=length(ll[[1]])
          next
      } 
      i=i+1
   }
   minbin=min(ll[[1]]+ll[[2]])
   while(nbins > binlim  || minbin < minbinlim) {
      ## Reduce the number of bins by combining the bins with the 
      ## fewest elements.  very simplistic.
      lb=ll[[1]]+ll[[2]]   # Combine signal present & absent tables
      elim= which.min(lb) # which Bin to eliminate.
      if (elim==nbins) elim=elim-1
      if (elim>1) if (lb[elim-1]< lb[elim+1]) elim=elim-1
      ll=binelim(ll,elim)
      nbins=length(ll[[1]])  # Should be one less than it was before.
      minbin=min(ll[[1]]+ll[[2]])
   }
   return(ll)  # Return frequency data with fewer categories.
}

tblunroll=function(ll) {  
   ## Create score data given a table of frequencies.
   lapply(ll, function(il) {
      vs= as.numeric(names(il))
      if (!length(vs)) vs=1:length(il)
      inverse.rle( list(lengths=il,values=vs))
   })
}


mindistfit<-function(x,y) {
    ##  Least squares fit to a line.  An implementation of algorithm
    ## on Mathworld.
    n<- length(x)
      xm<- mean(x)
      ym<- mean(y);
      A<- ((sum(y^2)-ym*ym*n) - (sum(x^2)-n*xm*xm))/ (n*xm*ym-sum(x*y))/2;
      a<- -A+sqrt(A*A+1);
      b<- mean(y)-a*mean(x);
      c(a,b);  # return slope, intercept.
}

normabvals=function(counts,link=qnorm) {
    ## return approx estimates of binormal parameters
    c0=cumsum(counts[[2]])
    c1=cumsum(counts[[1]])
    q0=link(c0/c0[length(c0)])
    q1=link(c1/c1[length(c0)])
    wh <- is.finite(q0) & is.finite(q1)
    q0=q0[wh]
    q1=q1[wh]
    ab=mindistfit(q0,q1)
    return(c( -ab[2]/ab[1], 1/ab[1] ) )  #return mu, sigma estimates
}

binhint=function (counts) {  
  ## make reasonable initial guesses for the bin
  ## boundaries based on counts.
    nn=sum(counts)
    bins=cumsum(counts)/nn  # frequencies into probabilities
    names(bins)=NULL

    bins=bins[-length(bins)]   # Take the 1 off the end.

    ## Add small random wiggles to the bin boundaries
    ## that are near 0, or 1
    wh=c(F,(diff(bins)==0))
    bins[wh]=bins[wh]+runif(sum(wh),-1,1)/nn

    ## Take other 0s and 1s off.
    wh=bins>=1
    bins[wh]=1-3*runif(sum(wh))/nn
    wh=bins<=0
    bins[wh]=3*runif(sum(wh))/nn
    return(sort(bins))
}

## function that adds random variations to intial parameters
## so that multiple optimizations will have different initial
## starting points.
parredist=function(x) rgamma(1,scale=sqrt(x),shape=sqrt(x))

binredist=function(bins) {  
   ## Take existing bin boundaries, and perturb them a bit.
  # Don't let bins boundaries be <=0 or >=1
     bins[ bins<=0] = 1e-2
     bins[ bins>=1] = 1.-1e-2
     sort(pnorm(qnorm(bins)+rnorm(length(bins),sd=4/length(bins))))
#    sort(bins+2*runif(length(bins),min=-bins,max=1-bins)/length(bins))
#    bnew=rbinom(length(hint[-1]), size=300, prob=hint[-1])/300
#    sort((bnew-.5)*.99+.5)
}

TwoSampleOrdinalRegression=function(x1=NA,x0=NA,counts1=NA,counts0=NA,
                                    models=c("PowerLaw")) {
  ## Input: x1=signal present scores, x0= signal absent scores
  ## Alternatively, counts1 is the number of signal present observations
  ## in each ordinal category, and counts0 is the number of signal absent
  ## observations in those same categories.
  if (any(is.na(x1))+any(is.na(x0))+any(is.na(counts1))+any(is.na(counts0))>2)
    stop("Please provide more data or data without NAs")
  if (length(counts1)!=length(counts0) )
    stop("Number of observational categories must be equal")

  if (any(is.na(counts1))||any(is.na(counts0)) ) {
    countsA=tablit(list(x1,x0))  # make frequency table
  } else {
    countsA=list(counts1,counts0)
  }
  counts=tablit2(countsA)  # make table for likelihood maximization.
        
  fits=lapply(models, function(mod){
    modu=toupper(mod)
    ffunc=switch(modu,EMPIRICAL=empirical,
      POWERLAW=plawfit,BIGAMMA=bigamfit,BINORMAL=binormalfit,
      BINORMALEQVAR=eqVarbinormalfit,CBM=CBMfit,CHISQ=bichisqfit,
      PROPROC=bichisqfit,
      MAXSIGNAL=maxObsfit,BILOGISTIC=logisticfit2,
      BILOGISTICEQVAR=logisticfit1, DUALBETA=dualBetafit)
    f=try( ffunc(),silent=T)
    if (inherits(f,'try-error')) {
      warning(paste("Model",mod,"not recognized"))
      return(list())
    }
    # Find optimal paramters for this model
    opt=lmaxit(counts, f)

    ## get auc and variance thereof
    aucstd=aucvalcalc( opt, f)
    return(list(parameters=opt$par[1:f$nparams],opt=opt,aucstd=aucstd,
                f=f,counts=countsA,model=mod))
  })

  x1x0=tblunroll(countsA)
  aucvs=MWWStatVar(x1x0[[1]], x1x0[[2]])
  # Do empirical "fit"
  emp=list(list( opt=list(value=NA,par=NA), f=list(nparams=NA),
    aucstd=c(aucvs[1], sqrt(aucvs[2])),counts=countsA,  model="Empirical"))
  fits=c(fits,emp)
  

  class(fits) <- c("TSORfits",class(fits))
  return(fits)
}

print.TSORfits <-function(fits) {
     # Provide a little table
     tblchar="%15s %8s %8s %8s %8s %8s %8s\n"
     cat(sprintf(tblchar,"Model", "NParams*", "LogLik*",
                   "AUC", "Std.Err", "AIC*", "BIC*"))
     tmp <- sapply(fits, function(fit) {
        cat(sprintf("%15s %8i %8.3f %8.3f ",
                 fit$model,fit$f$nparams, fit$opt$value,fit$aucstd[1]))
        cat(sprintf("%8.3f %8.3f %8.3f\n",
                  fit$aucstd[2], 2*(length(fit$opt$par)-fit$opt$value),
               -2*fit$opt$value + length(fit$opt$par)*log(sum(unlist(fit$counts)))))

     })
}         

plot.TSORfits=function(fits,file=NA,new=T,pch=1) {
  ## Make a plot with the empirical data and the fits
  nfit=length(fits)
  if (is.character(file)) pdf(file, width=4,height=4)
  # plot empirical data.
  plotroc(fits[[nfit]]$counts,type='b',new=new,pch=pch) 
  tpfs=((1:299)/300)^1.5
  ccc= c( '#000066FF','#00006699','#00006666','#00006644',
    '#660000FF','#66000099','#66000066','#66000044',
    '#006600FF','#00660099','#00660066','#00660044'
    )
  ltype=c( 'solid','12','41','4212')
  lwd=c(1,3,5,7)

  llist=list(fits[[nfit]]$model)
  for (i in 1:(nfit-1)) {  #Loop over fits, plot each.
    fit=fits[[i]]
    # print(fit$model)
    cj=(i)%%(length(ccc))+1
    lj=(i)%%(length(ltype))+1
    probs=fit$f$probabilities( c(fit$opt$par[1:fit$f$nparams],tpfs))
    lines(1-probs$p0,1-probs$p1,col=ccc[cj],lty=ltype[lj],lwd=lwd[lj])
    llist=c(llist,paste(fit$model," Model", sep=""))
  }
  legend(.4,.4,llist,col=ccc,lty=c( ltype),lwd=c(lwd),
         pch=c(pch,rep(NA,length(ccc))))
  if (is.character(file)) dev.off()
  
}

testTSOR=function() {  # Test my code
  t1=1:8      # Make some fake ordinal categorical counts.
  t0=8:1
  
  fits=TwoSampleOrdinalRegression(counts1=t1,counts0=t0,
    models=c("PowerLaw", "BiGamma",  "BiNormal", "BiNormalEqVar", "CBM",
      "ChiSq", "PropROC", "MaxSignal", "BiLogistic",
      "BiLogisticEqVar", "DualBeta"))
#     models=c( "PowerLaw",  "MaxSignal", "BiLogistic"))
#      models=c( "PowerLaw",  "ChiSq", "BiLogistic")) 
  return(fits)
    x1=round(log(rexp(100,1/5)))
    x0=round(log(rexp(100)))
#    fits=TwoSampleOrdinalRegression(x1,x0,models=c("PowerLaw","BiGamma","BiNormal","BiNormalEqVar","CBM", "PropROC", "BiLogistic","BiLogistic1"))
    fits=TwoSampleOrdinalRegression(x1,x0,models=c("BiNormal", "MaxSignal"))
#    fits=TwoSampleOrdinalRegression(x1,x0,models=c("DualBeta"))
#    x1=c(rnorm(150),rnorm(150)*.2+1)
#    x0=rnorm(300)
#    fits=TwoSampleOrdinalRegression(x1,x0,models=c("CBM"))
#    return(fits[[1]][[2]])
     return(fits)
}

aucvalcalc=function( opt, f) {  # Calculate AUC and its standard error
  pars=opt$par
  auc=f$aucfrompars(pars)    # Calculate AUC


  sigma=try( solve(-opt$hessian))  # Estimate of covariance matrix.
  if (inherits(sigma,"try-error")) {   # It barfed.
    delsd=NA
  } else {
    ## Numerical Delta method variance estimate.
    del= min(abs(pars))/1.e6
    if (del<1.e-9) del=1.e-9
    delauc=(apply(diag(length(pars))*del,1, function(d) f$aucfrompars(pars+d))
            -auc)/del
    delsd=sqrt( delauc%*% sigma %*% delauc)
  }
  return(c(auc,delsd))
}

plawfit=function(){
  ## Return functions necessary for doing a power-law fit to the data.

  nparams=1;  # This is a one parameter model
  ## hint is an initial guess of the parameters over which we are maximizing.
  ## It consists of the model parameters and the FNF values of the
  ## category thresholds.
  hint <- function(counts) {
    ## Use empirical AUC to generate initial guess of power index in power law.
    auc=AUCcat(counts)   
    return(c(auc/(1-auc),binhint(counts[[1]])))
  }

  newhint <- function(ohint,counts) {  # Return a perturbed hint.
     ## newhint() returns a hint vector like the original, but randomly
     ## perturbed, including both the model parameters and the thresholds.
     ## Randomly perturb the power law value for different optimizer runs.
     newlambda=parredist(ohint[1]*2)+0.7
     newbins=binredist(ohint[-1])
     return(c(newlambda, newbins))
   }

  checkbadpars <- function(pars) { 
     ## This is called at each step of the optimizer to see if the
     ## model parameters are outside allowed limits.
     ## Return TRUE if parameters are bad.
     if (pars[1]<0) return(TRUE) # Check scale, shape
     if (bincheck(pars[-1])) return(TRUE)  # Bin thresholds
     return(FALSE)  # Not bad
  }

  FNF2TNF <- function(x, pars)  1-(1-x)^pars[1] # TNF as function of FNF.
  
  probabilities <- function(pars)  {
      ## This function takes the vector of parameters over which we are
      ## maximizing and returns a list of FNF values and TNF values
      ## of the category thresholds.
      ## Note that usually the vector is mostly the FNF values.
      sigps=c(0,pars[-1],1)
      list( p1=sigps, p0=FNF2TNF(sigps,pars))   # Power law model
    }

  aucfrompars1 <- function(pars) { # Return auc as a function of the parameters
    integrate(FNF2TNF,0,1,pars=pars)$value
  }
 
  aucfrompars <- function(pars) { # Return auc as a function of the parameters
    p=pars[1];
    return( p/(1+p))
  }
  

  return(list(nparams=nparams,hint=hint,newhint=newhint,
              checkbadpars=checkbadpars,FNF2TNF=FNF2TNF,
              probabilities=probabilities, aucfrompars=aucfrompars))
}

bigamfit=function(){
  ## Return functions necessary for doing a bi-normal fit to the data.
  ## See the Dorfman reference.

  nparams=2
  hint <- function(counts) {
    auc=AUCcat(counts) # Use the empirical AUC to generate initial power index
    # I need better estimate of the second parameter
    return(c(1/(auc/(1-auc)),1.0,binhint(counts[[1]])))
  }

  newhint <- function(ohint,counts) {  # Return a perturbed hint.
     ## randomly perturb the initial or current parameter guesses.
     newlambda=parredist(ohint[1])+0.5
     newshape=parredist(ohint[2])
     newbins=binredist(ohint[-(1:2)])
     return(c(newlambda,newshape, newbins))
  }

  # Return TRUE if parameters are bad.
  checkbadpars <- function(pars) { 
     ## This is called at each step of the optimizer to see if the
     ## model parameters are outside allowed limits.
     ## Return TRUE if parameters are bad.
     if (pars[1]<0 || pars[2]<0) return(TRUE) # Check scale, shape
     if (bincheck(pars[-(1:2)])) return(TRUE)  # Bin thresholds
     return(FALSE)  # Not bad
  }

  FNF2TNF<- function(x, pars) pgamma(qgamma(x, shape=pars[2],rate=pars[1]),
           shape=pars[2])  # TNF as function of FNF.
  
  probabilities <- function(pars)  {
    ## Return bin boundaries (probability values) 
    ## as function of the model (bigamma) and its parameters 
    ## for signal absent and signal present data
    sigps=c(0,pars[-(1:2)],1)
    list( p1=sigps, p0=FNF2TNF(sigps,pars))
  }
  
  aucfrompars <- function(pars)  # Return auc as a function of the parameters
    integrate(FNF2TNF,0,1,pars=pars)$value  # Inefficient
  

  return(list(nparams=nparams,hint=hint,newhint=newhint,
              checkbadpars=checkbadpars,FNF2TNF=FNF2TNF,
              probabilities=probabilities, aucfrompars=aucfrompars))

}

dualBetafit=function(){
  ## Return functions necessary for doing a dual beta fit.
  ## See Mossman & Peng 2015.

  nparams=2
  hint <- function(counts) {
    auc=AUCcat(counts) # Use the empirical AUC to generate initial power index
    # I need better estimate of the second parameter
    return(c(2,2,binhint(counts[[1]])))
  }

  newhint <- function(ohint,counts) {  # Return a perturbed hint.
     ## randomly perturb the initial or current parameter guesses.
     newlambda=parredist(ohint[1])+1.0
     newshape=parredist(ohint[2])+1.0
     newbins=binredist(ohint[-(1:2)])
     return(c(newlambda,newshape, newbins))
  }

  # Return TRUE if parameters are bad.
  checkbadpars <- function(pars) { 
     ## This is called at each step of the optimizer to see if the
     ## model parameters are outside allowed limits.
     ## Return TRUE if parameters are bad.
     if (pars[1]<1 || pars[2]<1) return(TRUE) # Check scale, shape
     if (bincheck(pars[-(1:2)])) return(TRUE)  # Bin thresholds
     return(FALSE)  # Not bad
  }

  FNF2TNF<- function(x, pars) 1-(1-x^(1/pars[1]))^pars[2]#TNF as function of FNF
  
  probabilities <- function(pars)  {
    ## Return bin boundaries (probability values) 
    ## as function of the model (bigamma) and its parameters 
    ## for signal absent and signal present data
    sigps=c(0,pars[-(1:2)],1)
    list( p1=sigps, p0=FNF2TNF(sigps,pars))
  }
  
  aucfrompars <- function(pars) { # Return auc as a function of the parameters
    ## There may be a more efficient expression for this.
    # integrate(FNF2TNF,0,1,pars=pars)$value # inefficient
    1-pars[2]* beta(pars[2],pars[1]+1) # Same results, faster.
    
  }

  return(list(nparams=nparams,hint=hint,newhint=newhint,
              checkbadpars=checkbadpars, FNF2TNF=FNF2TNF,
              probabilities=probabilities, aucfrompars=aucfrompars))

}

binormalfit=function(){
  ## Return functions necessary for doing a bi-normal fit to the data.
  ## See the Dorfman reference.

  nparams=2
  hint <- function(counts) {
      return( c(normabvals(counts) ,binhint(counts[[1]])))
  } 

  newhint <- function(hint,counts) {  # Return a perturbed hint.
     newlambda=parredist(hint[1])
     newshape=parredist(hint[2])
     newbins=binredist(hint[-(1:2)])
     return(c(newlambda,newshape, newbins))
  }

  # Return TRUE if parameters are bad.
  checkbadpars <- function(pars) { 
     if (pars[1]< -5. || pars[2]<0) return(TRUE) # Check scale, shape
     if (bincheck(pars[-(1:2)])) return(TRUE)  # Bin thresholds
     return(FALSE)  # Not bad
  }

  # Return background bin boundaries as function of parameters
  probabilities <- function(pars)  {
      sigps=c(0,pars[-(1:2)],1)
      list( p1=sigps, p0=pnorm(qnorm(sigps, mean=pars[1],sd=pars[2])))
  }

  aucfrompars <- function(pars) { # Return auc as a function of the parameters
    return( pnorm( pars[1]/sqrt(1+pars[2]^2)))  # Binormal AUC.
  }  

  return(list(nparams=nparams,hint=hint,newhint=newhint,
              checkbadpars=checkbadpars,
              probabilities=probabilities, aucfrompars=aucfrompars))
}

bichisqfit=function(){
  ## Return functions necessary for doing a bi-normal fit to the data.
  ## See the Dorfman reference.
  ## The model can be parameterized as the
  ## non-centrality parameters of two X^2 distn, each with 1 d.o.f.
  ## See Hillis, Statistics in Medicine, 2015.
  ## Actually the maximizer uses the parameters mu, sigma because
  ## they seem to be better behaved.
  ## These are transformed to theta,lambda parameters of Hillis later.
  nparams=2

  musig2thetalambda= function(musig) {
      theta= (musig[1]/(musig[2]^2-1))^2  # Using Hillis parameterization
      lambda= musig[2]^2
      return(c(theta,lambda))
  }

  hint <- function(counts) {
      musig=normabvals(counts)  # if the data were binormal,...
#      return( c( musig2thetalambda(musig),binhint(counts[[1]])))
      return( c( musig,binhint(counts[[1]])))
  } 

  newhint <- function(hint,counts) {  # Return a perturbed hint.
     newtheta=parredist(hint[1])
     newlambda=parredist(hint[2])
     newbins=binredist(hint[-(1:2)])
     return(c(newtheta,newlambda, newbins))
  }

  # Return TRUE if parameters are bad.
#  checkbadpars <- function(pars) {   # For theta, lambda model
#     if (pars[1]< 0. || pars[2]<0) return(TRUE) # Check scale, shape
#     if (pars[1] > 1000 || pars[2]>100) return(TRUE)
#     if (bincheck(pars[-(1:2)])) return(TRUE)  # Check bin thresholds
#     return(FALSE)  # Not bad
#  }

    # Return TRUE if parameters are bad.
  checkbadpars <- function(pars) {   # for mu, sigma paramters
     if (pars[1]< -5 || pars[2]<0) return(TRUE) # Check scale, shape
     if (bincheck(pars[-(1:2)])) return(TRUE)  # Bin thresholds
     ms=musig2thetalambda(pars)
     theta=ms[1]
     lambda=ms[2]
     if (theta< 0. || lambda<0) return(TRUE) # Check scale, shape
     ## if X-squared non-centrality is too large, bad stuff happens.
     # if (theta > 1000 || lambda>100) return(TRUE) 
     return(FALSE)  # Not bad
  }

  FNF2TNF <- function (p, pars) {  # TNF as function of FNF.
    #theta=pars[1]  #If we use Hillis parameters theta,lambda
    # lambda=pars[2]
    thetalambda=musig2thetalambda(pars)
    theta=thetalambda[1]
    lambda=thetalambda[2]
     #cat( "    ", theta, lambda, "\r ")
     if (abs(pars[2]-1)<0.05 && TRUE) {
        ## Inconveniently the PROPROC/Chi-squared model blows up at sigma=b=1.
        ## However, at sigma=b=1, it is the same as the binormal model.
        return(pnorm(qnorm(p, mean=pars[1],sd=pars[2])))
     }
     if (lambda > 1) return(pchisq( qchisq(p,1,theta*lambda)*lambda, 1, theta))
     return(1-pchisq(qchisq(1-p,1,theta*lambda)*lambda, 1, theta)) 
  }
  
  # Return background bin boundaries as function of parameters
  probabilities <- function(pars)  {
     sigps=c(0,pars[-(1:2)],1)
     list( p1=sigps,p0=FNF2TNF(sigps,pars))
  }

  aucfrompars <- function(pars) { # Return auc as a function of the parameters
    # Do this two different ways
    # Using Hillis equation 16.  I think this works, but requires mvtnorm.
    if (F) if (require(mvtnorm)) {
       ms=musig2thetalambda(pars)
       th=ms[1]
       lam=ms[2]
       u1=sqrt(th)*(lam-1)/sqrt(lam+1); u2=sqrt(th)*sqrt(lam+1);rho=(lam-1)/(lam+1)
       rhom=matrix(c(1,rho,rho,1),ncol=2)
       aucprop=pmvnorm(upper=c(u1,u2),corr=rhom)+pmvnorm(upper=-c(u1,u2),corr=rhom)

       if (lam<1)  aucprop=1-aucprop
       return(aucprop)

     }
     a2=(integrate(FNF2TNF,0,1,pars=pars)$value)   # The slow? way
     return(c(a2))
  }  

  return(list(nparams=nparams,hint=hint,newhint=newhint,
              checkbadpars=checkbadpars, FNF2TNF=FNF2TNF,
              probabilities=probabilities, aucfrompars=aucfrompars))
}

logisticfit2=function() {
  # The two parameter logistic model

  nparams=2
  hint <- function(counts) {
      return( c(normabvals(counts,link=qlogis) ,binhint(counts[[1]])))
  } 

  newhint <- function(hint,counts) {  # Return a perturbed hint.
     newmu=parredist(hint[1])
     newshape=parredist(hint[2])
     newbins=binredist(hint[-(1:2)])
     return(c(newmu,newshape, newbins))
  }

  # Return TRUE if parameters are bad.
  checkbadpars <- function(pars) { 
     if (pars[1]< -5. || pars[2]<0) return(TRUE) # Check scale, shape
     if (bincheck(pars[-(1:2)])) return(TRUE)  # Bin thresholds
     return(FALSE)  # Not bad
  }

  FNF2TNF <- function (p, pars) {  # TNF as function of FNF.
    plogis( qlogis(p,location=pars[1],scale=pars[2]))
  }
  

  # Return background bin boundaries as function of parameters
  probabilities <- function(pars)  {
      sigps=c(0,pars[-(1:2)],1)
      list( p1=sigps, p0=FNF2TNF(sigps,pars))
  }

  aucfrompars <- function(pars) { # Return auc as a function of the parameters
    (integrate(FNF2TNF,0,1,pars=pars)$value)   # The slow? way
  }  

  return(list(nparams=nparams,hint=hint,newhint=newhint,
              checkbadpars=checkbadpars,
              probabilities=probabilities, aucfrompars=aucfrompars))
}


logisticfit1=function() {
  # A two logistic distribution with equal scales.
  nparams=1
  hint <- function(counts) {
      return( c(normabvals(counts,link=qlogis)[1] ,binhint(counts[[1]])))
  } 
  newhint <- function(hint,counts) {  # Return a perturbed hint.
     newmu=parredist(hint[1])
     newbins=binredist(hint[-(1)])
     return(c(newmu, newbins))
  }
  checkbadpars <- function(pars) { 
     if (pars[1]< -5.) return(TRUE) # Check scale, shape
     if (bincheck(pars[-(1)])) return(TRUE)  # Bin thresholds
     return(FALSE)  # Not bad
  }

  FNF2TNF <- function (p, pars) {  # TNF as function of FNF.
    plogis( qlogis(p,location=pars[1]))
  }
  probabilities <- function(pars)  {
    sigps=c(0,pars[-1],1)
    list( p1=sigps, p0=FNF2TNF(sigps,pars))
  }
  aucfrompars <- function(pars) { # Return auc as a function of the parameters
    mu=pars[1]
    or=exp(-mu)
    (1-or*(1+mu))/(1-or)^2  # This was derived by David Brown.
  }  
  return(list(nparams=nparams,hint=hint,newhint=newhint,
              checkbadpars=checkbadpars,FNF2TNF=FNF2TNF,
              probabilities=probabilities, aucfrompars=aucfrompars))

}

maxObsfit=function(nsig=1) {
  ## This uses a maximum signal model.  It assumes that the score from
  ## a signal-absent observation is the maximum of N Gaussian noise
  ## realizations.  The score from a signal-present observation is the
  ## maximum of N Gaussian noise realizations or nsig actual signals
  ## (of magnitude sked) in the image.
  ## Parameters: N and mu, and the probabilities of the signal
  ## ABSENT probabilities.  This is different from the other functions.
  ## This function does not work well.  I think I need a better initial guess.

  nparams=2
  hint<- function(counts) {
    auc=AUCcat(counts)
    ## Initial guess: A very rough initial guess
    return( c( 20, qnorm(auc)*sqrt(2)*0.9 ,binhint(counts[[2]])))
  } 

  newhint <- function(hint,counts) {  # Return a perturbed hint.
       newbins=binredist(hint[-(1:2)])
       return(c(parredist(hint[1]),parredist(hint[2]), newbins))
  }

  # Return TRUE if parameters are bad.
  checkbadpars <- function(pars) {
    if (pars[1]< 1) return(TRUE) # Check N
    if (pars[2]<0) return(TRUE) # Check skeauc
    if (bincheck(pars[-(1:2)])) return(TRUE)  # Bin thresholds
    return(FALSE)  # Not bad
  }

  TNF2FNF <- function(abps,pars) {  # Take TNF, return FNF
    N=pars[1]
    sked=pars[2]
    xx=qnorm(abps^(1/N))
    pnorm(xx)^(N-nsig)*pnorm(xx-sked)^nsig
  }
  # Return signal present! bin boundaries as a function of parameters
  probabilities <- function(pars)  {
    abps=c(0,pars[-(1:2)],1)
    return(list( p1=TNF2FNF(abps,pars), p0=abps))
  }

  aucfrompars <- function(pars) { # Return auc as a function of the parameters
    ## There is probably a more efficient expression for this.
    1-integrate(TNF2FNF,0,1,pars=pars)$value
  }
  
  return(list(nparams=nparams,hint=hint,newhint=newhint,
              checkbadpars=checkbadpars,TNF2FNF=TNF2FNF,
              probabilities=probabilities, aucfrompars=aucfrompars))
}

eqVarbinormalfit=function() {
  ## Equal variance binormal model fit. (one parameter)
  nparams=1
  hint<- function(counts) {
    auc=AUCcat(counts)
    return( c( qnorm(auc)*sqrt(2),binhint(counts[[1]])))
  } 

  newhint <- function(hint,counts) {  # Return a perturbed hint.
     newlambda=parredist(hint[1])
     newbins=binredist(hint[-1])
     ret=c(newlambda, newbins)
     return(ret)
  }

  # Return TRUE if parameters are bad.
  checkbadpars <- function(pars) { 
     if (pars[1]< 0. ) return(TRUE) # Check mean
     if (bincheck(pars[-(1)])) return(TRUE)  # Bin thresholds
     return(FALSE)  # Not bad
  }

  # Return background bin boundaries as function of parameters
  probabilities <- function(pars)  {
      sigps=c(0,pars[-1],1)
      ret= list( p1=sigps, p0=pnorm(qnorm(sigps, mean=pars[1])))
#      cat(" returned Pvalues "); print(ret)
      return(ret)
  }
  ## Return auc as a function of the parameters
  aucfrompars <- function(pars) pnorm( pars[1]/sqrt(2))

  return(list(nparams=nparams,hint=hint,newhint=newhint,
              checkbadpars=checkbadpars,
              probabilities=probabilities, aucfrompars=aucfrompars))
}


CBMfit=function(){
  ## fit data to the contaminated binormal model (CBM)
  ## Note that this does not currently work, it is incomplete.  

  # parameter 1 is the separation of the two normals
  # parameter 2 is the fraction of the signal
  nparams=2
  
  hint<- function(counts) {
    auc=AUCcat(counts)  # Poor initial guess.
    hint= c( qnorm(auc)*sqrt(2) , 0.5 ,binhint(counts[[1]]))
  } 

  newhint <- function(hint,counts) {  # Return a perturbed hint.
    auc=AUCcat(counts)
    while(TRUE) {
      newalpha=runif(1,.1,1)  # Random guess.
      dd=(auc-.5*(1-newalpha))/newalpha
      if ( 0<= dd && dd <=1 ) break
    }
    newdp=sqrt(2)*qnorm((auc-.5*(1-newalpha))/newalpha)
    newbins=binredist(hint[-(1:2)])
    return(c(newdp,newalpha, newbins))
  }

  # Return TRUE if parameters are bad.
  checkbadpars <- function(pars) { 
     if (pars[1]< 0. || pars[1]>10. || pars[2]<0 || pars[2]>1) return(TRUE) 
     if (bincheck(pars[-(1:2)])) return(TRUE)  # Bin thresholds
     return(FALSE)  # Not bad
  }

  # Return  bin boundaries as function of parameters
  probabilities <- function(pars)  {
    ##  can't invert the signal present cumulative distn
    sigabs=c(0,pars[-(1:2)],1)
    dp=pars[1]
    alpha=pars[2]
    x=qnorm(sigabs)
    list( p1=pnorm(x)*(1-alpha)+alpha*(pnorm(x-dp)), p0=sigabs)
  }
  aucfrompars <- function(pars) {  # CBM auc.
    dp=pars[1]
    alpha=pars[2]
    (1-alpha)*.5 + alpha* pnorm( pars[1]/sqrt(2))
  }

  return(list(nparams=nparams,hint=hint,newhint=newhint,
              checkbadpars=checkbadpars,
              probabilities=probabilities, aucfrompars=aucfrompars))
}

bincheck=function(bins) {  
  ## Optimizer calls this function at every step.  
  ## Return TRUE if bins take nonsensical values
  if (any(!is.finite(bins))) return (TRUE) # nonfinite/NA bins not allowed
  if (any(diff(bins)<0)) return(TRUE)
  if (any(bins<0) || any(bins>1)) return(TRUE) # bins are probabilities.
  return(FALSE)
}

lmaxit=function(counts,mdl ) {# hint,newhint,checkbadpars,probabilities) {
  ## This function maximizes the multinomial likelihood for purposes
  ## of ordinal regression.
  ## This function should only be called from the appropriate fit
  ## function, which sets up the proper variables and functions.
  ## This function performs multiple optimization attempts with 
  ## different random initial parameter and bin values, and picks
  ## the maximum of all of those.

  sigcount=counts[[1]] ; backcount=counts[[2]]
  ## The constant term on the front of the log likelihood.
  fscale= logfactorial(sum(sigcount))-sum(logfactorial(sigcount)) +
        logfactorial(sum(backcount))-sum(logfactorial(backcount))

  ## The R optim() function can use different algorithms for 
  ## maximization.  We usually use BFGS, because it is fast, 
  ## but occasionally throw in another for variety.
#  meths= c("Nelder-Mead","BFGS", "CG", "L-BFGS-B","SANN")
  meths= c("Nelder-Mead","BFGS", "L-BFGS-B","SANN")
  nmeth=length(meths)
  lastoptvalue = -9e19;
  hint=mdl$hint(counts)

#  cat("\n-> ")
  opts= lapply(1:49, function(iiii) { # A loop: Repeat optimization many times
    first <<- TRUE    # global variable...
    if (iiii<=nmeth) {
      optstart=hint # Start from exact hint value the first times
    } else {  
      if (lastoptvalue> -1e4  && iiii%%(nmeth*2-1)==1) { 
        ## Occasionally start different optimizer from last optimal value.
        optstart=lastoptpar; #cat(" * ")
      } else {  ## start optimizer with somewhat randomized parameters.
        if (runif(1)>.5) {
          optstart=mdl$newhint(hint,counts)  # randomize from original guess
        } else {
          optstart=mdl$newhint(lastoptpar,counts) # from current best guess.
        }
      }
    }
#    print(c("o",iiii));
#    cat (" } ",optstart)
    if (FALSE)
      if (runif(1)<.09) {
        cat("\r                                        \r")
      } else { cat("*")}
    ## Call optimizer which will attempt to maximize the function llfunc.
    method=meths[(iiii-1)%%nmeth+1]
    oo=try(optim( optstart, llfunc, hessian=T, method=method,
      control=list(fnscale=-1,maxit=10000),
      sigcount=sigcount, backcount=backcount,fscale=fscale,
      checkbadpars=mdl$checkbadpars,probabilities=mdl$probabilities
      ))
    ## If it barfed,
    if (inherits(oo,"try-error")) oo=list(convergence=40,value=lastoptvalue)
 #   cat(oo$par[1]," ",oo$value,'\n')
    ## Check to see if optimizer converged.
    # cat (oo$value, method, oo$par[1:2],'\n')
    if (oo$convergence==0 && oo$value>-9e5) 
       if (oo$value>lastoptvalue)  {
         lastoptvalue <<- oo$value
         lastoptpar <<- oo$par
       }
#    if (oo$convergence==0 && oo$value < -9e7) {print(optstart)}
    return(oo)
  })
  #cat("\r                                        \r")

  ## Which of the optimizations produced the highest likelihood?
  maxes=sapply(opts, function(i) i$value) 
  goods=sapply(opts, function(i) return(i$convergence)) #ones that did converge
#  print(maxes-goods*9e99)
  mmm=which.max(maxes-goods*9e99)   # Pick highest likelihood that did converge.

#  print(mmm)
  return(opts[[mmm]])  #Return optimal optimzation 

}


llfunc=function(pars, sigcount, backcount,fscale,checkbadpars,probabilities) {
   ## This function calculates the multinomial log likelihood
   ## of a particular model.  
   ## optim() calls this function repeatedly, trying to maximize it.

   rbad=-9e9 # Return this value if it's really bad..

   if (checkbadpars(pars)) {  # Check for invalid parameter values
     if (first) {
       ## if the first function call gives bad parameters, debug...
       # print('bad parameters')
       # cat(" pars: ",pars,"=>",rbad,'\n')
       first<<- FALSE
     }
     return(rbad)
   }
   
   # The probabilities are dictated by the model and bin boundaries.
   probs=probabilities(pars)  
   for (p in probs) if ( bincheck(p)) { 
     # Ack! optim gave us bad bin values return a bad number.
     if (first && FALSE) {
       print('bad probabilities')
       print(pars)
       print(probs)
     first<<- FALSE
     }
     return(rbad)
   }

   ## Calculate the multinomial log likelihood
   wh1=(c(sigcount==0, backcount==0))
   tosum=c(sigcount*log(diff(sort(probs$p1))),
        backcount*log(diff(sort(probs$p0))))
   tosum[wh1]=0
   tosum[!is.finite(tosum)]=rbad   # Check for bad nonsensical values

   ret=(  sum(tosum,na.rm=T) + fscale)  # our return value, the log ML value
#      (logfactorial(sum(sigcount))-sum(logfactorial(sigcount)) +
#        logfactorial(sum(backcount))-sum(logfactorial(backcount)) )
#   cat("  ", ret, "     \r")
#   lines(1-probs$p0,1-probs$p1,lty=3,col=sample(1:10,1));#locator(1)
   if (first) {
#     cat(' ',ret,' ')
     if (ret< -9e6) {
       ## Again, print some debugging info if this is the optimizer's first call
       if (F) {
         print("badtosum")
         print(pars)
         print(probs)
         print(sigcount); print(backcount)
         print(tosum) ;
       }
       first <<- FALSE
     }
   }
   first <<- FALSE
   return(ret)   # return multinomial log likelihood
}

##############################################################
### Above here are the routines that do the ordinal regression
### and likelihood maximiztaion.
##############################################################

plotroc=function(counts,type='p',pch=1,cex=1.3,new=T) {
  toplot=ROCxyCounts(counts)

  par(mar=c(2.7, 2.7, .5, .5), mgp=c(1.5,.6,0),lab=c(10,10,7))

  if (new) {
    plot(toplot,  xlab="False Positive Fraction", xlim=c(0,1),ylim=c(0,1),
         ylab="True Positive Fraction",cex=cex,type=type,pch=pch)
  } else {
    points(toplot, cex=cex,type=type,pch=pch)
  }
}

