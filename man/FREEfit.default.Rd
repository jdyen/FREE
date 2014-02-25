\name{FREEfit.default}
\alias{FREEfit.default}
\title{FUNCTION REGRESSION IN ECOLOGY AND EVOLUTION}
\description{This is the default method for models of class \code{FREEfit}. Fits a regression model for functional responses and one or more scalar covariates. Several different methods are possible, drawing on several existing R packages. Three methods are available that are unique to this package.}
\usage{
FREEfit.default(y, x, bins = 1:ncol(y), method = c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"), ...)
}
\arguments{
  \item{y}{Matrix of response data, sites in rows.}
  \item{x}{Matrix of predictor variables, sites in rows}
  \item{bins}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{method}{String defining the method to be used. Seven options are available but only six of these are currently implemented. Defaults to gamboost. Methods BUGSspline and BUGSjump can be slow, as can gamboost, INLA and stan if inappropriate settings are used.}
  \item{\dots}{Additional arguments to be passed to FREEmethod. See ?FREEfda, ?FREEgamboost, ?FREEinla, ?FREEstan, ?FREEbugs, ?FREEbugsJump and ?FREEfree for details.}
}
\details{This is an internal method, typically called using FREEfit.formula(). If using directly, the data y and x must be in an appropriate format. y must have one row for each replicate, with a column for each argument value. Columns must correspond to the same argument value across all replicates. x must have one row per replicate, with columns corresponding to covariates or predictor variables. bins defaults to 1:ncol(y), but providing appropriate bin values (argument values, corresponding to the columns of y) is recommended.}
\value{
Object of class \code{FREEfit} with the following components:
  \item{fitted}{fitted y values}
  \item{observed}{observed y data}
  \item{coefs.mean}{mean values for all fitted coefficients (including the mean function)}
  \item{coefs.sd}{sd values for all fitted coefficients (including the mean function)}
  \item{r2}{r-squared value for the fitted model, based on Pearson R}
  \item{family}{family used for model fitting}
  \item{bins}{vector of argument values for the response function}
  \item{xIC}{DIC value from the fitted model}
  \item{method}{the method used to fit model}
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\seealso{
\code{\link{FREEfit.formula}}
\code{\link{FREEfit.default}}
\code{\link{FREEfda}}
\code{\link{FREEgamboost}}
\code{\link{FREEinla}}
\code{\link{FREEstan}}
\code{\link{FREEbugs}}
\code{\link{FREEbugsJump}}
\code{\link{FREEfree}}
\code{\link{FREEfitCV}}
}
\examples{
  data <- FREEdataSim()
  data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
  mf <- model.frame(formula=y ~ VAR1 + VAR2, data=data.list)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  x <- x[, -1]
  y <- model.response(mf)
  bins <- 1:ncol(y)
  model <- FREEfit.default(y, x, bins, method='fda')
  summary(model)
  plot(model)
  \dontrun{
  	par(mfrow=c(2, 2))
    plot(data$mean.real, type='l', xlab="Argument", ylab="Response", ylim=c(min(data$mean.real, model$coefs.mean[1, ]), max(data$mean.real, model$coefs.mean[1, ])))
    lines(model$coefs.mean[1, ], col='red')
    plot(data$beta.real[1, ], type='l', xlab="Argument", ylab="Beta", ylim=c(min(data$beta.real[1, ], model$coefs.mean[2, ]), max(data$beta.real[1, ], model$coefs.mean[2, ])))
    lines(model$coefs.mean[2, ], col='red')
    plot(data$beta.real[2, ], type='l', xlab="Argument", ylab="Beta", ylim=c(min(data$beta.real[2, ], model$coefs.mean[3, ]), max(data$beta.real[2, ], model$coefs.mean[3, ])))
    lines(model$coefs.mean[3, ], col='red')
  }
}
\keyword{function regression}
\keyword{Bayesian}
\keyword{ecology}
\keyword{evolution}