print.FREEfit <-
function(x, ...){
  if (!is.null(x$call)) {
    call.print <- x$call
  } else {
    call.print <- "Formula interface must be used to return call..."
  }
  cat("Call:\n")
  print(call.print)
  cat("\nMethod:\n")
  print(x$method)
  cat("\nCoefficients:\n")
  print(x$coefs.mean)
  cat("\nr2: ")
  cat(x$r2)
  if (!is.null(x$xIC)) {
    cat("\nxIC: ")
    cat(x$xIC)
  } else {
    cat("\nxIC: xIC not computed for method ", x$method, ".....", sep="")
  }
}
