\name{predict.FREEfit}
\alias{predict.FREEfit}
\title{PREDICTIONS FOR FUNCTION REGRESSION MODELS}
\description{Predict method for objects of class \code{FREEfit}.}
\usage{
predict.FREEfit(object, newdata = NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{FREEfit}}
  \item{newdata}{A named list of newdata, with names corresponding the data used to fit the \code{FREEfit} model.}
  \item{\dots}{Additional arguments NOT CURRENTLY USED}
}
\value{
A matrix containing predicted values, with replicates in rows.
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  data <- FREEdataSim()
  data.list <- list(y=data$y.iid, VAR1=data$x$VAR1, VAR2=data$x$VAR2)
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method='fda')
  predict(model, newdata=list(VAR1=rnorm(10), VAR2=rnorm(10)))
}