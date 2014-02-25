\name{MakeCovMat}
\alias{MakeCovMat}
\title{CREATE A POSITIVE SEMI-DEFINITE COVARIANCE MATRIX}
\description{INTERNAL USE ONLY: creates a covariance matrix for use by \code{FREEdataSim}.}
\usage{
MakeCovMat(n = 10)
}
\arguments{
  \item{n}{Dimensions of the square covariance matrix}
}
\value{
A matrix of dimensions n by n
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  \dontrun{
    testMat <- MakeCovMat()
  }
}