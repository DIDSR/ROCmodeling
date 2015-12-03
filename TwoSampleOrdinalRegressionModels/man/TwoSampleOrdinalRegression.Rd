\name{TwoSampleOrdinalRegression}
\alias{TwoSampleOrdinalRegression}
\title{Fitting Two Sample Data with Maximum Likelihood
  Ordinal Regression}

\description{
  This routine performs maximum likelihood
  ordinal regression on two samples of data for many
  published models.  This is useful for ROC analysis,
  calculating AUC and its error, and plotting ROC curves.
}

\usage{
TwoSampleOrdinalRegression(x1=NA, x0=NA, counts1=NA, counts0=NA,
                           models=c(model1,model2,...))
}

\arguments{
  \item{x1}{a vector of values from sample 1.  Generally this
	sample is assumed to have larger values than  sample 0.
	It can be specified in conjunction with \code{x0} instead of
	using the count arguments.}
  \item{x0}{a vector of values from sample 0.  Generally this
	sample is assumed to have smaller values than sample 1.
	It can be specified in conjunction with \code{x1} instead of
	using the count arguments.}
  \item{counts1}{a vector of values of length \eqn{N} indicating the
	number of observations in each of \eqn{N} ordered categories
	from  sample 1.
	It can be specified in conjunction with \code{counts0} instead of
	using the \code{x1}, \code{x0} arguments.}
  \item{counts0}{a vector of values of length \eqn{N} indicating the
	number of observations in each of \eqn{N} ordered categories
	from sample 0.  Proportionally it should have more
	observations in lower categories than \code{counts1}.
	It can be specified in conjunction with \code{counts1} instead of
	using the \code{x1}, \code{x0} arguments.}
  \item{models}{a vector of names of models that should be fit to the
	data.  Currently possible model names are
    "PowerLaw", "BiGamma", "BiNormal",
    "BiNormalEqVar", "CBM", "ChiSq", "PropROC", "MaxSignal",
	"BiLogistic", "BiLogisticEqVar", and "DualBeta".  These are explained
	below.
  }
  
}

\details{
  \code{TwoSampleOrdinalRegression} takes as input two
  vectors of observed values (\code{x1}, \code{x0}), or two
  vectors of frequencies or counts (\code{counts1}, \code{counts0})
  of the number of observations in each ordered category for
  each of two samples.  It also takes a vector of names of models
  that should be fit to the data via maximum likelihood ordinal
  regression.

  \code{TwoSampleOrdinalRegression}
  performs maximum likelihood
  ordinal regression on two samples of data for many
  published models. 
  It uses the methods
  described in Dorfman & Alf 1969 and
  Metz et al. 1998.
  The two samples are assumed to be selected and  monotonically
  transformed from two  distributions, \eqn{F} and \eqn{G}.
  Equivalently we are assuming a relationship
  between the two
  cumulative distribution functions, i.e. the ROC curve,
  \deqn{ TPF= ROC(FPF)= 1 - F(G^{-1}(1-FPF|k) | h). }
  The derivative of this equation is the likelihood ratio, \eqn{dF/dG}.
  
  
  Let \eqn{n}=\code{counts1} and \eqn{m}=\code{counts0}.
  Then \eqn{n_i}{n[i]} observations are assumed to
  be discretized values of \eqn{X} from distribution 1,
  \eqn{F(x|h)},
  when
  \eqn{t_{i-1} < X < t_i}{t[i-1] < X < t[i]}, where \eqn{t} are
  latent thresholds.
  \eqn{m_i}{m[i]} is the number of deviates in the \eqn{i}th
  ordered category assumed to be discretized values of \eqn{X}
  from distribution 0, \eqn{G(x|k)}, between the same thresholds.
  \eqn{h} and \eqn{k} collectively are the model parameters. 
  
  The probability of obtaining a deviate from \eqn{F} that is
  placed in the \eqn{i}th category is
  \eqn{p_i=F(t_{i+1}|h)-F(t_i|h)}{p[i]=F(t[i+1]|h)-F(t[i]|h)},
  and from \eqn{G} it is
  \eqn{q_i=G(t_{i+1}|k)-G(t_i|k)}{q[i]=G(t[i+1]|k)-G(t[i]|k)}.
  The multinomial likelihood of obtaining this set of data given
  the model parameters \eqn{h, k} and \eqn{t} thresholds is
  \deqn{ \left(\sum_i{n_i}\right)! \left(\sum_i{m_i}\right)! \prod_i \frac{q_i^{m_i}}{m_i!}\frac{p_i^{n_i}}{n_i!}}{
         sum(n)! sum(m)! prod_i ( q[i]^{m[i]} * p[i]^{n[i]}/m[i]!/n[i]!)}
  The values of \eqn{t}, \eqn{h}, and \eqn{k} that maximize
  this likelihood are the ML estimates.


  Possible models for \eqn{F} and \eqn{G} are:
  \enumerate{
    \item  "PowerLaw".  This model assumes that the underlying
    survival functions, or true and false positive fractions,
    are related by a power law.  This is the case if the
    distributions are Pareto or exponential.  The parameter
	of the fit is the power index.  See Egan 1975 or
	Samuelson & He 2014.
    \item "BiGamma".  This model assumes that the distributions
    are gamma distributions with equal shape parameters but
    different rate parameters.  See Dorfman et al. 1996
    for details.  The parameters in the model are (1) the
	rate parameter of \eqn{F}, and (2) the
	shape parameter of both distributions.
    \item "BiNormal".  The distributions are modeled as 
    two Gaussians with different means and variances.
	See Dorfman & Alf 1969.
	The parameters of the model are (1) the differences in means
	between the distributions, and (2) the ratios of the standard
	deviations between distributions.
    \item "BiNormalEqVar".  The distributions are modeled as
    two Gaussians with equal variance, but different means.
    This mean difference is the single model parameter.
    \item "CBM", the contaminated binormal model.  Under this model
    distribution G is a standard 
	Gaussian (mean=0, std=1),   and distribution \eqn{F}
    is a mixed model of two Gaussians, a standard Gaussian and
    a Gaussian with a non-zero mean.  See Dorfman & Berbaum 2000.
    The parameters are (1) the mean of the non-zero mean Gaussian
    and (2) the mixture fraction of distribution \eqn{F}.
    \item "ChiSq". This models the data as two chi-squared distributions,
    each with one degree of freedom and its own non-centrality
    parameter.  See Hillis 2015.  The computed parameters are
    not the non-centrality values, but the corresponding
	normal model parameters
    because they are more numerically stable.
    \item "PropROC". This model is the same as the "ChiSq" model.
    See Metz & Pan 1999.
    \item "BiLogistic".  This option assumes that both samples came from
    logistic distributions.  The parameters of this model are
    (1) the difference of location parameters and (2) the ratio
    of the scale parameters between the distibutions.
    \item "BiLogisticEqVar".  This is the same as the previous
    model, but the scale parameter of both logistic distributions is
    the same.
    The fit parameter is the difference in means of the
    logistic distributions.
    \item "DualBeta". See Mossman & Peng 2015.
	The parameters are the two power indices, \eqn{\alpha}
    and \eqn{\beta} for	the two beta distributions.
    \item "MaxSignal".  Under this model distribution \eqn{G}, is
    assumed to be the maximum of M standard Gaussian noise deviates.
    Distribution \eqn{F} is assumed to be the maximum
    of \eqn{M-1} standard Gaussian noise deviates or a Gaussian with unit
    variance and mean \eqn{\mu}.  \eqn{M} and \eqn{\mu} are
    the model parameters.
  }

  For each model the optimization is performed
  many times at different initial values.  This is done because
  local optima are sometimes encountered.   Therefore the
  code runs very slowly, but it is very robust.
  The code is relatively simple and those looking to add a
  new model for two sample ordinal regression can easily
  add their model to  this code.
  For additional questions, consult the code or email the
  maintainer.

}

\value{
  \code{TwoSampleOrdinalRegression}  returns an object of
  class \code{TSORfits}, which is really just a list of
  the fits to the models.  Each model fit is a list of parameters
  describing the fit and what functions were used,
  including the output of the optimizer \code{optim}
  on the multinomial likelihood function.

  This returned value may be passed to \code{print} for a summary
  of the fits, or passed to \code{plot} to make a plot the
  model ROC curves.

  When printed the software includes the name of the model,
  the number of parameters,  the maximum log-likelihood,
  the AUC value and its estimated standard error, the 
  Akaike information criterion (AIC), and Bayesian information
  criterion (BIC) for each model.
  For example, see Samuelson & He 2014.

  The number of parameters, log-likelihood, AIC and BIC are
  marked with asterisks (*).  These values should only be used
  to compare models from this same program,
  and not compared with the output
  of other ordinial regression programs,   unless you really know what
  you are doing.  This is because the number of parameters does not
  include the number of
  latent modeled thresholds, which are not full free parameters
  in the usual sense.  This number of thresholds
  may be different in this program  than it is in others.
  Currently the number of latent modeled thresholds
  is at most 12, though that can be easily changed in the code.

}

%\seealso{
%  \code{\link{save}}.
%}
\examples{
library(TwoSampleOrdinalRegressionModels)

## Make some fake data
x1=round(log(rexp(100,1/5)))
x0=round(log(rexp(100)))
## Fit a power law model and a two-parameter bi-normal model to data
fits=TwoSampleOrdinalRegression(x1,x0,models=c("PowerLaw", "BiNormal"))
print(fits)   # How good are the fits?  What is the AUC?
plot(fits)    # Look at the model ROC curves
print(fits[[2]]$aucstd)     # The binormal model AUC and error
print(fits[[1]]$parameters) # Print the power-law index.

t1=1:8      # Make some fake ordinal categorical counts.
t0=8:1
fits=TwoSampleOrdinalRegression(counts1=t1,counts0=t0,
                                models=c("BiNormalEqVar", "BiLogisticEqVar"))
}

\keyword{regression}

