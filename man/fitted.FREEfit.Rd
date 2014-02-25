\name{fitted.FREEfit}
\alias{fitted.FREEfit}
\title{FITTED VALUES FOR CLASS FREEfit}
\description{Return a matrix of fitted values from a fitted FREEfit model}
\usage{
fitted.FREEfit(object, ...)
}
\arguments{
  \item{object}{An object of class FREEfit, typically returned from a call to FREEfit}
  \item{\dots}{Additional arguments: not currently used...}
}
\value{Fitted values from an object of class FREEfit, in a matrix format with one row per site}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  data <- FREEdataSim()
  data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method="fda")
  fitted <- fitted(model)
  model.r2 <- cor(c(fitted(model)), c(data.list$y)) * cor(c(fitted(model)), c(data.list$y))
  plot(c(fitted) ~ c(data.list$y))
}