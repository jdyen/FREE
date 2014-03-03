\name{FREE-package}
\alias{FREE-package}
\alias{FREE}
\docType{package}
\title{Function Regression in Ecology and Evolution}
\description{
\pkg{FREE} is a set of helper functions for fitting models where response variables are functions rather than scalar values. These functions draw on several existing approaches in R, some of which were developed for functional data and some that were not. Extensions for models with functional predictors variables are in development. There is a focus on ecological and evolutionary models, with aids for model interpretation and model validation.
}
\details{
\tabular{ll}{
Package: \tab FREE\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-02-21\cr
License: \tab GPL-3\cr
}
\pkg{FREE} is centred around two main functions: \code{\link{FREEfit}} and \code{\link{FREEfitCV}}. \code{\link{FREEfit}} has a formula interface and takes functional response data (in matrix format, replicates in rows) and scalar predictor variables (in matrix format, replicates in rows). \code{\link{FREEfit}} fits a model using a method selected by the user (options: \code{gamboost}, \code{fda}, \code{INLA}, \code{stan}, \code{BUGSjump}, \code{BUGSspline}). These methods each have their own strengths and weaknesses, as discussed in \cite{Yen et al. (in preparation)}.

\code{\link{FREEfitCV}} performs k-fold cross validation on models, fitted using \code{\link{FREEfit}}. This is useful for model validation and for assessment and tuning of model parameters.
}
\author{
  Jian Yen \cr
  Maintainer: Jian Yen \email{jdl.yen@gmail.com}
}
\references{
  \cite{Yen et al. (in preparation) Functional regression in ecology and evolution}
}
\keyword{ package }
\keyword{ function valued data }
\keyword{ functional response }
\keyword{ longitudinal data analysis }
\keyword{ Bayesian }
\examples{
# Fit a functional regression model using the fda package and plot the results using default plot methods for objects of class FREEfit
data <- FREEdataSim()
data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method="fda")
summary(model)
plot(model)
}
