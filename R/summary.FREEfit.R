##' @method summary FREEfit
##' @export
summary.FREEfit <-
function(object, ...){
  if (!is.null(object$call)) {
    call.print <- object$call
  } else {
    call.print <- "formula interface must be used to return call"
  }
  if (!is.null(object$DIC)) {
    DIC.print <- object$DIC
  } else {
    DIC.print <- paste("DIC not computed for method ", object$method, sep="")
  }
  res <- list(call=call.print, method=object$method, coefficients=object$coefs.mean,
              coefficient.sd=object$coefs.sd, r2=object$r2, DIC=DIC.print)
  if (!is.null(object$delta.mean)) {
    res$delta.mean <- object$delta.mean
  }
  if ((!is.null(object$gamma.mean)) & (object$method == "scalar")) {
    res$gamma.mean <- object$gamma.mean
  }
  if ((!is.null(object$rand.coefs.mean))) {
    res$rand.coefs <- object$rand.coefs.mean
  }
  if (any(unlist(object$rhats) > 1.1, na.rm=TRUE)) {
    warn.print <- "some rhat values were greater than 1.1; consider increasing n.iters"
    res$warning <- warn.print
  }
  res
}
