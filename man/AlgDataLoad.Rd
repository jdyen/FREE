\name{AlgDataLoad}
\alias{AlgDataLoad}
\title{LOAD NAWQA DIATOM DATA}
\description{Load NAWQA diatom data for a given subset of sites. Log-scale or untransformed data are available.}
\usage{
AlgDataLoad(data.subset = NULL, log.scale = TRUE)
}
\arguments{
  \item{data.subset}{A vector identifying the sites to be returned}
  \item{log.scale}{Logical specifying whether returned response data should be log-transformed or not}
}
\details{Internal function for package development only - TO BE DELETED}
\value{
  A list with four items.
  \item{y.data}{The response data in matrix form, with one site per row}
  \item{X.data}{The predictor data for each site. Six variables returned: carbon, orthophosphate, nitrogen, light, elevation and drainage area}
  \item{bin.data}{The bin values for the data}
  \item{coord.data}{The coordinates of each site}
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data.subset = NULL, log.scale = TRUE) 
{
    work.dir.store <- getwd()
    if (Sys.info()[["sysname"]] == "Windows") {
        setwd("C:/Users/jdyen/Dropbox/JianFiles")
    }
    else {
        setwd("/Users/jian/Dropbox/JianFiles")
    }
    stdx <- function(x) {
        return((x - mean(x, na.rm = T))/sd(x, na.rm = T))
    }
    response <- read.csv("y_alg_log_test.csv", row.names = 1)
    response <- as.matrix(response)
    preds.fish <- read.csv("x_alg_test.csv", row.names = 1)
    breakpoints <- read.csv("breakpoints_test.csv", row.names = 1)
    if (is.null(data.subset)) {
        data.subset <- 1:nrow(response)
    }
    else {
        data.subset <- data.subset
    }
    preds.fish[data.subset, c(1:7, 10:11)] <- apply(preds.fish[data.subset, 
        c(1:7, 10:11)], 2, stdx)
    predictors <- as.matrix(preds.fish[data.subset, c(2, 3, 6, 
        7, 10, 11)])
    coords <- as.matrix(preds.fish[data.subset, c(8, 9)])
    if (log.scale) {
        size.class <- c(as.matrix(breakpoints[2, ]))
    }
    else {
        size.class <- c(as.matrix(breakpoints[1, ]))
    }
    size.class <- size.class[1:ncol(response)] + diff(size.class[1:2])/2
    y.data <- response[data.subset, ]
    X.data <- predictors
    bin.data <- size.class
    coord.data <- coords
    colnames(X.data) <- c("CARBON", "PHOS", "NITRO", "LIGHT", 
        "ELEV", "DRAIN")
    colnames(coord.data) <- c("LATITUDE", "LONGITUDE")
    X.data <- sweep(X.data, 2, apply(X.data, 2, max), "/")
    if (log.scale) {
        y.data <- log(y.data + 1)
    }
    setwd(work.dir.store)
    return(list(y.data = y.data, X.data = X.data, bin.data = bin.data, 
        coord.data = coord.data))
  }
}
\keyword{dataset}
\keyword{diatoms}