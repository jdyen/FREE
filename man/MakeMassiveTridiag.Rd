\name{MakeMassiveTridiag}
\alias{MakeMassiveTridiag}
\title{MAKE LARGE PSEUDO-DIAGONAL MATRIX FROM ONE SMALLER MATRIX}
\description{Internal use only: repeats a smaller matrix several times on the main diagonal of a large matrix to make a massive pseudo-diagonal matrix. Typically called by \code{\link{FREEfit}} for method gamboost}
\usage{
MakeMassiveTridiag(tridiag, n.times)
}
\arguments{
  \item{tridiag}{Smaller matrix to be repeated}
  \item{n.times}{Number of times to repeat the smaller matrix}
}
\value{
A large pseudo-diagonal matrix.
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  small.mat <- MakeTridiag(c(1, 2, 2, 2, 1), upper=-1, lower=-1)
  big.mat <- MakeMassiveTridiag(small.mat, n.times=5)
  big.mat
}