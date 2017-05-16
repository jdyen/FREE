##' @method plot FREEfitCV
##' @export
plot.FREEfitCV <-
function(x, ...){
  par(mfrow=c(1,1))
  plot(x$observed, x$predicted,
       bty='l', las = 1,
       xlab = "Observed", ylab = "Predicted",
       pch = 16, col = grey(0.6))
}
