##' @method residuals FREEfit
##' @export
residuals.FREEfit <-
function(object, ...){
  resid <- object$fitted - object$observed
  resid
}
