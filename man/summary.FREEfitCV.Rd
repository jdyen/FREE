\name{summary.FREEfitCV}
\alias{summary.FREEfitCV}
\title{SUMMARISE CROSS VALIDATED FUNCTION REGRESSION MODEL}
\description{Summary method for objects of class \code{FREEfitCV}.}
\usage{
summary.FREEfitCV(object, ...)
}
\arguments{
  \item{object}{An object of class \code{FREEfitCV}}
  \item{\dots}{Additional arguments NOT CURRENTLY USED}
}
\details{
Provides a summary of a model cross validation.
}
\value{
A named list with the following elements:
\item{call}{the model call}
\item{observed}{the observed response data}
\item{predicted}{the predicted values for the held out data}
\item{cv.r2}{the r2 value for the predicted versus observed data}
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  data <- FREEdataSim()
  data.list <- list(y=data$y.iid, VAR1=data$x$VAR1, VAR2=data$x$VAR2)
  modelCV <- FREEfitCV(y ~ VAR1 + VAR2, data=data.list, method='fda', n.cv=5)
  summary(modelCV)
}