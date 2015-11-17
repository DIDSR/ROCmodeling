# ROCmodeling

This is the code for an R package that contains a routine that performs maximum likelihood ordinal regression on two samples of data for many published models.  This is useful for ROC analysis, calculating AUC and its error, and plotting ROC curves.

BUILDING THE R PACKAGE:

To build the R source package from this source tree, execute the commands in makepackage, e.g.
sh makepackage


INSTALLING FROM THE PACKAGE:

To install the package in R, download the .tar.gz file, change to the folder or directory that contains that source package and execute the following function:
> install.packages("TwoSampleOrdinalRegressionModels_0.1.tar.gz", type="source")

For information about the package:
> library(help="TwoSampleOrdinalRegressionModels")

To use the package in R, load the library:
> library(TwoSampleOrdinalRegressionModels)

Read how to use the function:
> help(TwoSampleOrdinalRegression)


