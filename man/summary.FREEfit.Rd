\name{summary.FREEfit}
\alias{summary.FREEfit}
\title{SUMMARISE FUNCTION REGRESSION MODEL}
\description{Summary method for objects of class \code{FREEfit}.}
\usage{
summary.FREEfit(object, ...)
}
\arguments{
  \item{object}{An object of class \code{FREEfit}}
  \item{\dots}{Additional arguments NOT CURRENTLY USED}
}
\details{
Provides a summary of a fitted \code{FREEfit} model.
}
\value{
A named list with the following elements:
\item{call}{the model call}
\item{method}{the method used in model fitting}
\item{coefficients}{the mean coefficient values}
\item{coefficient.sd}{the standard deviations for coefficient estimates}
\item{r2}{the r2 value for the fitted model}
\item{xIC}{The available information criterion value (AIC or DIC where possible)}
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  data <- FREEdataSim()
  data.list <- list(y=data$y.iid, VAR1=data$x$VAR1, VAR2=data$x$VAR2)
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method='fda')
  summary(model)
}