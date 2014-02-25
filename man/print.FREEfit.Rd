\name{print.FREEfit}
\alias{print.FREEfit}
\title{PRINT FUNCTION REGRESSION MODEL}
\description{Print method for objects of class \code{FREEfit}.}
\usage{
print.FREEfit(x, ...)
}
\arguments{
  \item{x}{An object of class \code{FREEfit}}
  \item{\dots}{Additional arguments NOT CURRENTLY USED}
}
\details{Prints the main model details to the console.}
\value{
Printed output in console.
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  data <- FREEdataSim()
  data.list <- list(y=data$y.iid, VAR1=data$x$VAR1, VAR2=data$x$VAR2)
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method='fda')
  print(model)
}