\name{ConvertB2A}
\alias{ConvertB2A}
\title{CONVERT DATA FROM MATRIX TO VECTOR FORMAT}
\description{FREEgamboost and FREEinla require data in a vector format. These functions call ConvertB2A to perform this conversion.}
\usage{
ConvertB2A(y.matrix, X.matrix, bin.vector, coord.data = NULL)
}
\arguments{
  \item{y.matrix}{y.matrix is the matrix of response functions, with one row per site}
  \item{X.matrix}{X.matrix is the matrix of predictor variables, with one row per site}
  \item{bin.vector}{Single vector of the bin values}
  \item{coord.data}{Optional matrix of coordinates, with one row per site}
}
\details{This function performs a very specific data format conversion and should not be used beyond this.}
\value{
  \item{y.vector}{A vector of response values, with p rows per site, where each function has p observed values}
  \item{X.vector}{A matrix of predictor variables, where the values for each site are repeated p times}
  \item{sites.vector}{A vector of site identification codes}
  \item{bin.vector}{A vector of bin values, repeated once for each site}
  \item{coord.vector}{Optional matrix of coordinates, where the values for each site are repeated p times}
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y.matrix, X.matrix, bin.vector, coord.data = NULL) 
{
    if (nrow(y.matrix) != nrow(X.matrix)) {
        stop("y.matrix and X.matrix must have the same number of rows.....", 
            call. = FALSE)
    }
    if (ncol(y.matrix) != length(bin.vector)) {
        stop("y.matrix must have the same number of columns as there are bins.....", 
            call. = FALSE)
    }
    n.rows <- nrow(y.matrix)
    n.bins <- length(bin.vector)
    y.vector <- as.vector(t(y.matrix))
    sites.vector <- rep(1:n.rows, each = n.bins)
    X.vector <- X.matrix[sites.vector, ]
    X.vector <- as.matrix(X.vector)
    rownames(X.vector) <- 1:nrow(X.vector)
    colnames(X.vector) <- colnames(X.matrix)
    bin.vector <- rep(bin.vector, times = n.rows)
    if (!is.null(coord.data)) {
        coord.vector <- coord.data[sites.vector, ]
        rownames(coord.vector) <- 1:nrow(coord.vector)
        return(list(y.vector = y.vector, X.vector = X.vector, 
            sites.vector = sites.vector, bin.vector = bin.vector, 
            coord.vector = coord.vector))
    }
    else {
        return(list(y.vector = y.vector, X.vector = X.vector, 
            sites.vector = sites.vector, bin.vector = bin.vector))
    }
  }
}