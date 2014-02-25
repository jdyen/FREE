FishDataLoad <-
function(data.subset=NULL, log.scale=TRUE){
  # Internal function for NAWQA data only - not for general use
  work.dir.store <- getwd()
  if (Sys.info()[["sysname"]] == "Windows") {
    setwd("C:/Users/jdyen/Dropbox/JianFiles")
  } else {
    setwd("/Users/jian/Dropbox/JianFiles")
  }
  stdx <- function(x){
    return((x - mean(x, na.rm=T)) / sd(x, na.rm=T))
  }
  response <- read.csv("y_fish_log_test.csv", row.names=1)
  response <- as.matrix(response)
  preds.fish <- read.csv("x_fish_test.csv", row.names=1)
  breakpoints <- read.csv("breakpoints_test.csv", row.names=1)
  if (is.null(data.subset)) {
    data.subset <- 1:nrow(response)
  } else {
    data.subset <- data.subset
  }
  preds.fish[data.subset, c(1:6, 9:10)] <- apply(preds.fish[data.subset, c(1:6, 9:10)],
                                                 2, stdx)
  predictors <- as.matrix(preds.fish[data.subset, c(2, 3, 6, 9, 10)])
  coords <- as.matrix(preds.fish[data.subset, c(7, 8)])
  if (log.scale) {
    size.class <- c(as.matrix(breakpoints[4, ]))
  } else {
    size.class <- c(as.matrix(breakpoints[3, ]))
  }
  size.class <- size.class[1:ncol(response)] + diff(size.class[1:2]) / 2
  y.data <- response[data.subset, ]
  X.data <- predictors
  bin.data <- size.class
  coord.data <- coords
  colnames(X.data) <- c("CARBON", "PHOS", "LIGHT", "ELEV", "DRAIN")
  colnames(coord.data) <- c("LATITUDE", "LONGITUDE")
  X.data <- sweep(X.data, 2, apply(X.data, 2, max), "/")
  if (log.scale) {
    y.data <- log(y.data + 1)
  }
  setwd(work.dir.store)
  return(list(y.data=y.data, X.data=X.data, bin.data=bin.data, coord.data=coord.data))
}
