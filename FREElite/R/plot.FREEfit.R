##' @method plot FREEfit
##' @export
plot.FREEfit <-
function(x, ...){
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
  n.plots <- nrow(vals$mean) + 4
  if (is.null(x$fp.sd.mean)) {
    n.plots <- n.plots - 1
  }
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
  if (!is.null(x$fp.sd.mean)) {
    sd.mean <- c(sqrt(x$sigma2.mean), x$fp.sd.mean)
    sd.mean <- ifelse(sd.mean < 0, 0, sd.mean)
    sd.upper <- c(sqrt(max(x$sigma2.mean + 1.96 * x$sigma2.sd, 0)), x$fp.sd.mean + 1.96 * x$fp.sd.sd)
    sd.upper <- ifelse(sd.upper < 0, 0, sd.upper)
    sd.lower <- c(sqrt(max(x$sigma2.mean - 1.96 * x$sigma2.sd, 0)), x$fp.sd.mean - 1.96 * x$fp.sd.sd)
    sd.lower <- ifelse(sd.lower < 0, 0, sd.lower)
    par(mar=c(5.1, 7.1, 4.1, 2.1))
    plot(sd.mean[1], 1, pch=16, cex=2, xlim=c(0, max(sd.upper)),
         ylim=c(0, length(sd.mean) + 1),
         las=1, bty='l', xlab="Standard deviations", yaxt="n", ylab="")
    lines(c(sd.lower[1], sd.upper[1]), c(1, 1))
    lines(c(0, 0), c(0, length(sd.mean) + 1), lty=2)
    for (i in 2:length(sd.mean)) {
      points(sd.mean[i], i, pch=16, cex=2)
      lines(c(sd.lower[i], sd.upper[i]),
            c(i, i))
    }
    axis(side=2, at=c(1:length(sd.mean)),
         labels=c("Residual", paste("Group", 1:length(x$fp.sd.mean), sep="")),
         las=1)
    mtext("Variance component", side=2, adj=0.5, line=5, cex=0.75)
  }
  par(old.par)
}
