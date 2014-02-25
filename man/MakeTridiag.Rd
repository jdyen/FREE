\name{MakeTridiag}
\alias{MakeTridiag}
\title{MAKE TRIDIAGONAL MATRIX}
\description{Internal use only: makes a tridiagonal matrix from vectors or scalars}
\usage{
MakeTridiag(diagon, upper, lower, nrow = NULL)
}
\arguments{
  \item{diagon}{Scalar or vector giving the main diagonal. If scalar then nrow must be !NULL.}
  \item{upper}{Scalar or vector giving the upper diagonal. If vector must have length \code{nrow - 1} or \code{length(diagon) - 1}.}
  \item{lower}{Scalar or vector giving the lower diagonal. If vector must have length \code{nrow - 1} or \code{length(diagon) - 1}.}
  \item{nrow}{Specifies the number of rows in the matrix. Only used if \code{length(diagon) == 1}.}
}
\value{
A tridiagonal matrix.
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  mat <- MakeTridiag(c(1, 2, 2, 2, 1), -1, 10)
  mat
}