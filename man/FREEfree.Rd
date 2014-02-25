\name{FREEfree}
\alias{FREEfree}
\title{BAYESIAN FUNCTION REGRESSION USING PURPOSE-BUILT MCMC}
\description{NOT CURRENTLY IMPLEMENTED: This will be a purpose-built MCMC approach for Bayesian functional regression, intended to replace the INLA, stan, fda, gamboost and BUGS methods inplemented in the FREE package. This is a work in progress and is not expected to be completed until 2015.}
\usage{
FREEfree(y, x, bins, ...)
}
\arguments{
  \item{y}{Matrix of response data, sites in rows.}
  \item{x}{Matrix of predictor variables, sites in rows}
  \item{bins}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{\dots}{Additional arguments to be passed to method \code{FREE}}
}
\details{NOT YET IMPLEMENTED}
\value{
An error: NOT YET IMPLEMENTED
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\seealso{
\code{\link{FREEfda}}
\code{\link{FREEgamboost}}
\code{\link{FREEinla}}
\code{\link{FREEstan}}
\code{\link{FREEbugs}}
\code{\link{FREEbugsJump}}
\code{\link{FREEfit}}
}
\examples{
  \dontrun{
    # MODEL NOT YET IMPLEMENTED
    data <- FREEdataSim()
    data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
    model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method="FREE", n.chains=2, n.iters=100)
  }
}
\keyword{Function regression}