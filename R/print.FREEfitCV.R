print.FREEfitCV <-
function(x, ...){
  if (!is.null(x$call)) {
    call.print <- x$call
  } else {
    call.print <- "Formula interface must be used to return call..."
  }
  cat("Call:\n")
  print(call.print)
  cat("\nCross validation r2:\n")
  print(x$cv.r2)
}
