summary.FREEfit <-
function(object, ...){
  if (!is.null(object$call)) {
    call.print <- object$call
  } else {
    call.print <- "formula interface must be used to return call"
  }
  if (!is.null(object$xIC)) {
    DIC.print <- object$xIC
  } else {
    DIC.print <- paste("DIC not computed for method ", object$method, sep="")
  }
  if (any(unlist(object$rhats) > 1.1)) {
    warn.print <- "some rhat values were greater than 1.1; consider increasing n.iters"
  }
  res <- list(call=call.print, method=object$method, coefficients=object$coefs.mean,
              coefficient.sd=object$coefs.sd, r2=object$r2, DIC=DIC.print,
              warning=warn.print)
  res
}
