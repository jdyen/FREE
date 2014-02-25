\name{plot.FREEfit}
\alias{plot.FREEfit}
\title{PLOT FUNCTIONAL REGRESSION MODEL}
\description{Plotting method for objects of class \code{FREEfit}.}
\usage{
plot.FREEfit(x, ...)
}
\arguments{
  \item{x}{An object of class \code{FREEfit}}
  \item{\dots}{Additional arguments to be passed to plot()}
}
\details{
Produces a plot of the fitted versus observed values, as well as plots of the intercept function and all covariate effects.
}
\value{
Plot of an object of class \code{FREEfit}.
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\seealso{
\code{\link{print.FREEfit}}
\code{\link{summary.FREEfit}}
\code{\link{predict.FREEfit}}
\code{\link{residuals.FREEfit}}
\code{\link{coef.FREEfit}}
\code{\link{fitted.FREEfit}}
}
\examples{
  data <- FREEdataSim()
  data.list <- list(y=data$y.iid, VAR1=data$x$VAR1, VAR2=data$x$VAR2)
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method='fda')
  plot(model)
}