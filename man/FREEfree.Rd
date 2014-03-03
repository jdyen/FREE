\name{FREEfree}
\alias{FREEfree}
\title{Bayesian function regression using purpose-built MCMC}
\description{NOT CURRENTLY IMPLEMENTED: \code{FREEfree} will implement a purpose-built MCMC approach for Bayesian functional regression, intended to replace the \code{\link{FREEfda}}, \code{\link{FREEgamboost}}, \code{\link{FREEinla}}, \code{\link{FREEstan}}, \code{\link{FREEbugs}}, \code{\link{FREEbugsJump}} methods currently implemented in \pkg{FREE}. This is a work in progress and is not expected to be completed until 2015.}
\usage{
FREEfree(y, x, bins, ...)
}
\arguments{
  \item{y}{Matrix of response data, sites in rows.}
  \item{x}{Matrix of predictor variables, sites in rows}
  \item{bins}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{\dots}{Additional arguments to be passed to \code{FREEfree}}
}
\details{NOT YET IMPLEMENTED}
\value{
An error: NOT YET IMPLEMENTED
}
\references{
  \cite{Yen et al. (in preparation) Functional regression in ecology and evolution}
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\seealso{
\code{\link{FREEfda}}, \code{\link{FREEgamboost}}, \code{\link{FREEinla}}, 
\code{\link{FREEstan}}, \code{\link{FREEbugs}}, \code{\link{FREEbugsJump}}, 
\code{\link{FREEfit}}
}
\examples{
\dontrun{
  # MODEL NOT YET IMPLEMENTED
  data <- FREEdataSim()
  data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method="FREE")
}
}
\keyword{Function regression}