\name{residuals.FREEfitCV}
\alias{residuals.FREEfitCV}
\title{CALCULATE RESIDUALS FOR CROSS VALIDATED FUNCTION REGRESSION MODEL}
\description{Residuals method for objects of class \code{FREEfitCV}.}
\usage{
residuals.FREEfitCV(object, ...)
}
\arguments{
  \item{object}{An object of class \code{FREEfitCV}}
  \item{\dots}{Additional arguments NOT CURRENTLY USED}
}
\details{Calculates residuals from a cross validated function regression model.}
\value{
Matrix of residuals with replicates in rows.
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  data <- FREEdataSim()
  data.list <- list(y=data$y.iid, VAR1=data$x$VAR1, VAR2=data$x$VAR2)
  modelCV <- FREEfitCV(y ~ VAR1 + VAR2, data=data.list, method='fda', n.cv=5)
  residuals(modelCV)
}