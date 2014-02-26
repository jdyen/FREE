\name{FREE-package}
\alias{FREE-package}
\alias{FREE}
\docType{package}
\title{FUNCTION REGRESSION IN ECOLOGY AND EVOLUTION}
\description{
Package FREE is a set of helper functions for fitting models where response variables are functions rather than scalar values. These functions draw on several existing approaches in R, some of which were developed for functional data and some that were not. Extensions for models with functional predictors variables are in development. There is a focus on ecological and evolutionary models, with aids for model interpretation and model validation.
}
\details{
\tabular{ll}{
Package: \tab FREE\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-02-21\cr
License: \tab GPL-3\cr
}
Package FREE is centred around two main functions: FREEfit and FREEfitCV. FREEfit has a formula interface and takes functional response data (in matrix format, replicates in rows) and scalar predictor variables (in matrix format, replicates in rows). FREEfit fits a model using a method selected by the user (options: gamboost, fda, INLA, stan, BUGSjump, BUGSspline). These methods each have their own strengths and weaknesses, as discussed in Yen et al. (in preparation) Functional regression in ecology and evolution.
}
\author{
  Jian Yen <jdl.yen@gmail.com>
  Maintainer: Jian Yen <jdl.yen@gmail.com>
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\keyword{ package }
\keyword{ function valued data }
\keyword{ functional response }
\keyword{ longitudinal data analysis }
\keyword{ Bayesian }
\examples{
  \dontrun{
    data <- FREEdataSim()
    data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
    model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method="fda")
    summary(model)
    plot(model)
  }
}
