\name{coef.FREEfit}
\alias{coef.FREEfit}
\title{COEFFICIENTS FROM CLASS FREEfit}
\description{Extracts coefficients from a FREEfit object, returning mean and confidence/credible intervals. Coefficients are returned in a matrix, with rows corresponding to the mean function followed by the variables in the order they were input.}
\usage{
coef.FREEfit(object, ...)
}
\arguments{
  \item{object}{An object of class \code{FREEfit}}
  \item{\dots}{Additional arguments: not currently used}
}
\value{
  \item{mean}{Mean value for each coefficient, returned as a matrix with one row per variable}
  \item{upper}{Upper confidence/credible interval for each coefficient}
  \item{lower}{Lower confidence/credible interval for each coefficient}
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  data <- FREEdataSim()
  data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method="fda")
  fitted.betas <- coef(model)
  plot(fitted.betas$mean[1,], type='l')
  lines(fitted.betas$upper[1,], lty=2)
  lines(fitted.betas$lower[1,], lty=2)
}