\name{print.FREEfitCV}
\alias{print.FREEfitCV}
\title{PRINT CROSS VALIDATED FUNCTION REGRESSION MODEL}
\description{Print method for an object of class \code{FREEfitCV}.}
\usage{
print.FREEfitCV(x, ...)
}
\arguments{
  \item{x}{An object of class \code{FREEfitCV}}
  \item{\dots}{Additional arguments NOT CURRENTLY USED}
}
\details{
Prints a summary of a cross validation for a \code{FREEfit} model
}
\value{
Printed summary in the console
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  data <- FREEdataSim()
  data.list <- list(y=data$y.iid, VAR1=data$x$VAR1, VAR2=data$x$VAR2)
  modelCV <- FREEfitCV(y ~ VAR1 + VAR2, data=data.list, method='fda', n.cv=5)
  print(modelCV)
}