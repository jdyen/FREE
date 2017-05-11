plot.FREEfitCV <-
function(x, ...){
  par(mfrow=c(1,1))
  plot(x$observed, x$predicted, bty='l', xlab="Observed", ylab="Predicted")
}
