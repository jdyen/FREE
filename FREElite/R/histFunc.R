hist_func <- function(y, breaks) {
  return(hist(y, breaks=breaks, plot=FALSE)$counts)
}