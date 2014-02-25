\name{residuals.FREEfit}
\alias{residuals.FREEfit}
\title{CALCULATE RESIDUALS FROM FUNCTION REGRESSION}
\description{Residuals method for objects of class \code{FREEfit}.}
\usage{
residuals.FREEfit(object, ...)
}
\arguments{
  \item{object}{An object of class \code{FREEfit}}
  \item{\dots}{Additional arguments NOT CURRENTLY USED}
}
\details{
Calculates fitted - observed for the fitted \code{FREEfit} model.
}
\value{
A matrix of residuals for the fitted model, with replicates in rows.
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  data <- FREEdataSim()
  data.list <- list(y=data$y.iid, VAR1=data$x$VAR1, VAR2=data$x$VAR2)
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method='fda')
  resid <- residuals(model)
  plot(resid ~ fitted(model))
}