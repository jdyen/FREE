coef.FREEfit <-
function(object, ...){
  if (!is.null(object$coefs.sd)) {
    upper <- object$coefs.mean + 1.96 * object$coefs.sd
    lower <- object$coefs.mean - 1.96 * object$coefs.sd
  } else {
    upper <- NULL
    lower <- NULL
  }
  coef <- list(mean=object$coefs.mean, upper=upper, lower=lower)
  coef
}
