\name{plotPretty}
\alias{plotPretty}
\title{Create plot of fitted model coefficients}
\description{Takes a fitted FREEfit model and creates a plot of all fitted function-valued coefficients}
\usage{
plotPretty(x, ...)
}
\arguments{
  \item{x}{fitted model of class \code{FREEfit}}
  \item{\dots}{additional arguments to be passed to \code{plot}}
}
\details{This function plots the fitted coefficients from a freefit model in a slightly nicer way than the standard plot function (and without all the extra information on convergence and mixing). This approach is not particularly flexible; looking at the source code for the function \code{FREE::plotPretty} will give some ideas of how to create custom plots of fitted FREEfit models}
\references{
  \cite{Yen et al. (2015) Function regression in ecology and evolution: FREE. Methods in Ecology and Evolution, 6: 17-26.}
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\seealso{
\code{\link{FREEfit}}
}
\examples{
data <- FREEdataSim(function.class="response", n=10)
data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, n.iters=15)
plotPretty(model)
}