plot.FREEfit <-
function(x, ...){
  vals <- coef(x)
  n.plots <- nrow(vals$mean) + 1
  dim.plots <- ceiling(n.plots / 2)
  par(mfrow=c(dim.plots, 2))
  plot(c(x$fitted), c(x$observed), bty='l', xlab="Fitted", ylab="Observed", las=1, ...)
  y.lab <- c("Mean", rep("Beta", {nrow(vals$mean) - 1}))
  for (i in 1:{n.plots - 1}) {
    if (!is.null(vals$upper) & !is.null(vals$lower)) {
      y.min <- min(vals$lower[i, ], na.rm=TRUE)
      y.max <- max(vals$upper[i, ], na.rm=TRUE)
    } else {
      y.min <- min(vals$mean[i, ], na.rm=TRUE)
      y.max <- max(vals$mean[i, ], na.rm=TRUE)
    }
    plot(vals$mean[i, ] ~ x$bins, type='l', bty='l', xlab="Argument", ylab=y.lab[i],
         ylim=c(y.min, y.max), las=1, ...)
    if (!is.null(vals$upper) & !is.null(vals$lower)) {
      lines(vals$upper[i, ] ~ x$bins, lty=2)
      lines(vals$lower[i, ] ~ x$bins, lty=2)
    }
    if (i > 1) {
      lines(c(min(x$bins), max(x$bins)), c(0, 0), lty=3)
    }
  }
  if (is.null(vals$upper) | is.null(vals$lower)) {
    cat("Estimates of parameter uncertainty not available...\n")
  }
}
