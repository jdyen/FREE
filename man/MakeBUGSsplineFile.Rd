\name{MakeBUGSsplineFile}
\alias{MakeBUGSsplineFile}
\title{MAKE MODEL FILE FOR WINBUGS 1.4 FUNCTIONAL REGRESSION MODEL}
\description{Internal function only: this function makes a default model file for WinBUGS 1.4 that can be used to fit a functional regression model. Typically this function would be called by \code{\link{FREEfit}} with method BUGSspline.}
\usage{
MakeBUGSsplineFile(filename = "FREEbugsSplineTemp.txt", ARmod = FALSE)
}
\arguments{
  \item{filename}{String giving the name of the model file, which will be created in the working directory}
  \item{ARmod}{Logical, should the observed responses have AR1 or IID errors? Defaults to FALSE, which generates code for IID errors.}
}
\details{Internal function only.}
\value{
A .txt file in the working directory.}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  \dontrun{
    bugs.file <- "TestBUGSfile.txt"
    MakeBUGSsplineFile(filename=bugs.file, ARmod=FALSE)
    # check working directory for a file called TestBUGSfile.txt
  }
}