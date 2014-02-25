\name{plot.FREEfitCV}
\alias{plot.FREEfitCV}
\title{PLOT CROSS VALIDATED FUNCTION REGRESSION RESULTS}
\description{Plotting method for objects of class \code{FREEfitCV}.}
\usage{
plot.FREEfitCV(x, ...)
}
\arguments{
  \item{x}{An object of class \code{FREEfitCV}}
  \item{\dots}{Additional arguments to be passed to plot()}
}
\details{
Produces a plot of predicted versus observed data returned from a cross validated function regression model.
}
\value{
Plot of an object of class \code{FREEfitCV}.
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  data <- FREEdataSim()
  data.list <- list(y=data$y.iid, VAR1=data$x$VAR1, VAR2=data$x$VAR2)
  modelCV <- FREEfitCV(y ~ VAR1 + VAR2, data=data.list, method='fda', n.cv=5)
  plot(modelCV)
}