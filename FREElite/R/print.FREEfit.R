print.FREEfit <-
function(x, ...){
  if (!is.null(x$call)) {
    call.print <- x$call
  } else {
    call.print <- "formula interface must be used to return call"
  }
  cat("Call:\n")
  print(call.print)
  cat("\nMethod:\n")
  print(x$method)
  cat("\nCoefficients:\n")
  print(x$coefs.mean)
  cat("\nr2: ")
  cat(x$r2)
  if (!is.null(x$DIC)) {
    cat("\nDIC: ")
    cat(x$DIC)
  } else {
    cat("\nDIC: DIC not computed for method ", x$method, ".....", sep="")
  }
  if (any(unlist(x$rhats) > 1.1)) {
    cat("\nWarning: some rhats were greater than 1.1; consider increasing n.iters")
  }
}
