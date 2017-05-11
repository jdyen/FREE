summary.FREEfit <-
function(object, ...){
  if (!is.null(object$call)) {
    call.print <- object$call
  } else {
    call.print <- "Formula interface must be used to return call..."
  }
  if (!is.null(object$xIC)) {
    xIC.print <- object$xIC
  } else {
    xIC.print <- paste("xIC not computed for method ", object$method, sep="")
  }
  res <- list(call=call.print, method=object$method, coefficients=object$coefs.mean,
              coefficient.sd=object$coefs.sd, r2=object$r2, xIC=xIC.print)
  res
}
