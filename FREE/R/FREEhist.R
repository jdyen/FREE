FREEhist <- function(y, n.bins=10, case.var=NULL, data=NULL) {
  if (!is.null(data$y)) {
    y <- data$y
  } else {
    y <- y
  }
  if (is.null(y)) {
    stop("Variable y must be provided...", call.=FALSE)
  }
  if (!is.null(case.var)) {
    if (!is.null(data$case.var)) {
      case.var <- data$case.var
    } else {
      case.var <- case.var
    }
  }
  if (is.null(case.var) & !is.list(y)) {
    stop("Sorting variable case.var must be provided if y is not a list...", call.=FALSE)
  }
  if (is.list(y)) {
  	min.bin <- min(sapply(y, min, na.rm=TRUE), na.rm=TRUE)
  	max.bin <- max(sapply(y, max, na.rm=TRUE), na.rm=TRUE)
    breaks <- seq(min.bin, max.bin, length={n.bins+1})
    n.obs <- length(y)
    out <- matrix(sapply(y, hist_func, breaks=breaks), nrow=n.obs, byrow=TRUE)
    mids <- hist(y[[1]], breaks=breaks, plot=FALSE)$mids
  } else {
    if (is.vector(y)) {
      min.bin <- min(y, na.rm=TRUE)
      max.bin <- max(y, na.rm=TRUE)
      breaks <- seq(min.bin, max.bin, length={n.bins+1})
      n.obs <- length(unique(case.var))
      if (is.null(case.var)) {
        stop("Sorting variable case.var must be provided for vector form of y...", call.=FALSE)
      }
      out <- tapply(y, case.var, hist_func, breaks=breaks)
      out <- matrix(unlist(out), nrow=n.obs, byrow=TRUE)
      mids <- hist(y[which(case.var==unique(case.var)[1])], breaks=breaks, plot=FALSE)$mids
    } else {
      stop("y must be a vector or list...", call.=FALSE)
    }
  }
  return(list(y.hist=out, breaks.hist=mids))
}