summary.FREEfitCV <-
function(object, ...){
  if (!is.null(object$call)) {
    call.print <- object$call
  } else {
    call.print <- "Formula interface must be used to return call..."
  }
  res <- list(call=call.print, observed=object$observed, predicted=object$predicted,
              cv.r2=object$cv.r2)
  res
}
