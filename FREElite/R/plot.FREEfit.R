plot.FREEfit <-
function(x, ...){
  vals <- coef(x)
  if (is.null(nrow(vals$mean))) {
    vals$mean <- matrix(vals$mean, nrow=1)
    vals$upper <- matrix(vals$upper, nrow=1)
    vals$lower <- matrix(vals$lower, nrow=1)
  }
  n.plots <- nrow(vals$mean) + 4
  dim.plots <- ceiling(n.plots / 2)
  par(mfrow=c(dim.plots, 2))
  plot(c(x$fitted), c(x$observed), bty='l', xlab="Fitted", ylab="Observed", las=1, ...)
  if (x$method == "scalar") {
    y.lab <- rep("Beta", {nrow(vals$mean)})
  } else {
    y.lab <- c("Mean", rep("Beta", {nrow(vals$mean) - 1}))
  }
  plot(x$llik_all[[1]], type='l', ylim=c(min(unlist(x$llik_all)), max(unlist(x$llik_all))),
       bty='l', xlab="Iteration", ylab="Log likelihood", las=1)
  for (i in 2:length(x$llik_all)) {
    lines(x$llik_all[[i]], col=i)
  }
  hist(unlist(x$rhats), main="", xlab="Rhat values")
  for (i in 1:nrow(vals$mean)) {
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
    if ((i > 1) | (x$method == "scalar")) {
      lines(c(min(x$bins), max(x$bins)), c(0, 0), lty=3)
    }
  }
  if (is.null(vals$upper) | is.null(vals$lower)) {
    cat("Estimates of parameter uncertainty not available...\n")
  }
  plot(x$sigma2.mean, 1, pch=16, cex=2, xlim=c(min(x$sigma2.mean - 1.96 * x$sigma2.sd,
                                                   x$sigma2_gamma.mean - 1.96 * x$sigma2_gamma.sd),
                                               max(x$sigma2.mean + 1.96 * x$sigma2.sd,
                                                   x$sigma2_gamma.mean + 1.96 * x$sigma2_gamma.sd)),
       ylim=c(0, length(x$sigma2_gamma.mean) + 2),
       las=1, bty='l', xlab="Standard deviations", yaxt="n")
  lines(c(x$sigma2.mean - 1.96 * x$sigma2.sd, x$sigma2.mean + 1.96 * x$sigma2.sd), c(1, 1))
  for (i in seq(along=x$sigma2_gamma.mean)) {
    points(x$sigma2_gamma.mean[i], (i + 1), pch=16, cex=2)
    lines(c(x$sigma2_gamma.mean[i] - 1.96 * x$sigma2_gamma.sd[i],
            x$sigma2_gamma.mean[i] + 1.96 * x$sigma2_gamma.sd[i]),
          c(i + 1, i + 1))
  }
  axis(side=2, at=c(1:(length(x$sigma2_gamma.mean) + 1),
       labels=c("Residual", paste("Group", 1:length(x$sigma2_gamma.mean), sep=""))))
  par(mfrow=c(1, 1))
}
