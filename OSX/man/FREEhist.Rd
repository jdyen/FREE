\name{FREEhist}
\alias{FREEhist}
\title{Convert observed values to a distribution using fixed number of bins}
\description{Takes a set of observed values for several cases/replicates and converts these values to a distribution using the \code{hist} function in R.}
\usage{
FREEhist(y, n.bins, case.var=NULL, data=NULL)
}
\arguments{
  \item{y}{Vector of response data, either a numeric vector or a list with each index representing one case}
  \item{n.bins}{Number of bins to be used when converting data to a distribution}
  \item{case.var}{If \code{y} is numeric; vector of case IDs identifying which observations belong to which case}
  \item{data}{Optional: list containing the variables \code{y} and \code{case.var}}
}
\details{This function is designed as a helper function when binned data are not available. Although binning data is not ideal, in many cases it is essential if we do not wish to assume a particular parametric form. Alternatives (e.g., kernel estimation) are possible, but the estimated distributions still need to be discretized prior to function regression analyses.}
\value{
List with the following components:
  \item{y.hist}{Binned version of \code{y}}
  \item{breaks.hist}{Midpoints of the bins used in data conversion}
}
\references{
  \cite{Yen et al. (in preparation) Functional regression in ecology and evolution}
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\seealso{
\code{\link{FREEfit}}
}
\examples{
\dontrun{
  data <- rnorm(10000, mean=0, sd=1)
  sites <- rep(1:10, each=1000)
  data.hist <- FREEhist(data, n.bins=10, case.var=sites)$y.hist
}
}