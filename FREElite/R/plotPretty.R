plotPretty <- function(x, ...) {
  if (class(x) != "FREEfit") {
    stop("plotPretty works with FREEfit models only", call.=FALSE)
  }
  old.par <- par()
  old.par$cin <- NULL
  old.par$cra <- NULL
  old.par$csi <- NULL
  old.par$cxy <- NULL
  old.par$din <- NULL
  old.par$page <- NULL  
  vals <- coef(x)
  if (is.null(nrow(vals$mean))) {
    vals$mean <- matrix(vals$mean, nrow=1)
    vals$upper <- matrix(vals$upper, nrow=1)
    vals$lower <- matrix(vals$lower, nrow=1)
  }
  n.plots <- nrow(vals$mean)
  dim.plots <- ceiling(n.plots / 2)
  par(mfrow=c(dim.plots, 2))
  if (x$method == "scalar") {
    par(mfrow=c(1, 1))
    plot(vals$mean ~ x$bins, type='l', bty='l', xlab="Argument",
         ylab="Effect on response",
         ylim=c(vals$lower, vals$upper), las=1, cex.lab=1.25, ...)
    polygon(c(x$bins, rev(x$bins)), c(vals$lower, rev(vals$upper)), col="gray75",
            border=NA)
    lines(vals$mean ~ x$bins)
    lines(c(min(x$bins), max(x$bins)), c(0, 0), lty=2)
  } else {
    y.lab <- c("Mean", rep("Effect on response", {nrow(vals$mean) - 1}))
    for (i in 1:nrow(vals$mean)) {
      if (!is.null(vals$upper) & !is.null(vals$lower)) {
        y.min <- min(vals$lower[i, ], na.rm=TRUE)
        y.max <- max(vals$upper[i, ], na.rm=TRUE)
      } else {
        y.min <- min(vals$mean[i, ], na.rm=TRUE)
        y.max <- max(vals$mean[i, ], na.rm=TRUE)
      }
      plot(vals$mean[i, ] ~ x$bins, type='l', bty='l', xlab="Argument",
           ylab=y.lab[i], ylim=c(y.min, y.max), las=1, cex.lab=1.25, ...)
      if (!is.null(vals$upper) & !is.null(vals$lower)) {
        polygon(c(x$bins, rev(x$bins)), c(vals$lower[i, ], rev(vals$upper[i, ])),
                col="gray75", border=NA)
        lines(vals$mean[i, ] ~ x$bins)
      }
      lines(c(min(x$bins), max(x$bins)), c(0, 0), lty=2)
    }
  }
  par(old.par)
}