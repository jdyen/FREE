\name{FREEcvInternal}
\alias{FREEcvInternal}
\title{CROSS VALIDATION FUNCTION FOR FREEFIT MODELS}
\description{Internal use only: cross validation function for use within \code{FREEfitCV}.}
\usage{
FREEcvInternal(i, n.obs, n.cv, y, x, bins, method, verbose, ...)
}
\arguments{
  \item{i}{Index variable}
  \item{n.obs}{number of replicates}
  \item{n.cv}{number of folds}
  \item{y}{response data}
  \item{x}{predictor data}
  \item{bins}{argument values}
  \item{method}{method to be used}
  \item{verbose}{logical: print progress to console?}
  \item{\dots}{additional arguments to be passed to FREEmethod}
}
\value{
A list of predicted and observed values
}
\author{
Jian Yen <jdl.yen@gmail.com>
}
\examples{
  \dontrun{
    # FIX THIS UP
  }
}