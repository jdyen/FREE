\name{stanCodeDefault}
\alias{stanCodeDefault}
\title{MAKE CODE FILE FOR STAN MODEL}
\description{Internal use only: makes the code file for a stan MCMC model.}
\usage{
stanCodeDefault()
}
\details{Internal use only.}
\value{
A model file for a stan model fit
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  \dontrun{
    stan.file <- stanCodeDefault()
    stan.file
  }
}