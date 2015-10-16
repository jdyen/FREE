\name{methods}
\alias{print.FREEfit}
\alias{summary.FREEfit}
\alias{plot.FREEfit}
\alias{fitted.FREEfit}
\alias{coef.FREEfit}
\alias{residuals.FREEfit}
\alias{predict.FREEfit}
\alias{print.FREEfitCV}
\alias{summary.FREEfitCV}
\alias{plot.FREEfitCV}
\alias{residuals.FREEfitCV}
\title{Methods for \code{FREEfit} models}
\description{Methods for objects of class \code{FREEfit} and class \code{FREEfitCV}.}
\usage{
\method{print}{FREEfit}(x, \dots)
\method{summary}{FREEfit}(object, \dots)
\method{plot}{FREEfit}(x, \dots)
\method{fitted}{FREEfit}(object, \dots)
\method{coef}{FREEfit}(object, \dots)
\method{residuals}{FREEfit}(object, \dots)
\method{predict}{FREEfit}(object, newdata=NULL, \dots)
\method{print}{FREEfitCV}(x, \dots)
\method{summary}{FREEfitCV}(object, \dots)
\method{plot}{FREEfitCV}(x, \dots)
\method{residuals}{FREEfitCV}(object, \dots)
}
\arguments{
  \item{x}{an object of class \code{FREEfit} or class \code{FREEfitCV}}
  \item{object}{an object of class \code{FREEfit} or class \code{FREEfitCV}}
  \item{newdata}{a named list of newdata, with names corresponding the data used to fit the \code{FREEfit} model.}
  \item{\dots}{additional arguments to be passed to method}
}
\details{
Basic methods for \code{FREEfit} models. The plot method gives plots of rhat values and the log-likelihood plotted against iterations, which are particularly useful for checking model convergence and mixing.
}
\references{
  \cite{Yen et al. (2015) Function regression in ecology and evolution: FREE. Methods in Ecology and Evolution, 6: 17-26.}
}
\author{Jian Yen \email{jdl.yen@gmail.com>}}
\seealso{
\code{\link{FREEfit}}
}
\examples{
# Fit a function regression model using the fda package and demonstrate the default methods
data <- FREEdataSim()
data.list <- list(y=data$y.iid, VAR1=data$x$VAR1, VAR2=data$x$VAR2)
model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, n.iters=5)
print(model)
summary(model)
plot(model)
fitted(model)
coef(model)
residuals(model)
predict(model, newdata=list(VAR1=rnorm(10), VAR2=rnorm(10)))
}