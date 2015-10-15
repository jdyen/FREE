residuals.FREEfitCV <-
function(object, ...){
  resid <- object$predicted - object$observed
  resid
}
