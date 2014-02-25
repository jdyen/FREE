\name{MakeBUGSfile}
\alias{MakeBUGSfile}
\title{CREATE WINBUGS MODEL FILE}
\description{Internal function only: create a standard model file for WinBUGS 1.4 that can be used to fit a spline-based functional regression model. All settings here are defined by the function \code{\link{FREEfit}} with method BUGSjump}
\usage{
MakeBUGSfile(Q, filename = "FREEbugsJumpTemp.txt", order = 3, cont = order, ARmod = FALSE)
}
\arguments{
  \item{Q}{Number of variables}
  \item{filename}{String giving the filename}
  \item{order}{Order of the splines}
  \item{cont}{Degree of continuity of splines}
  \item{ARmod}{Logical: should errors be AR1 or IID? Defaults to TRUE, which gives AR1 errors.}
}
\details{Internal function only.}
\value{
A file saved to the working directory, that can be called by R2WinBUGS.
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  \dontrun{
    bugs.file <- "TestBUGSfile.txt"
    MakeBUGSfile(Q=2, filename=bugs.file, order=3, cont=3, ARmod=TRUE)
    # check the working directory for a text file called "TestBUGSfile.txt"
  }
}