FREEbsInternal <-
function(i, n.bs, y, x, bins, verbose=TRUE, ...){
  y.store <- y
  x.store <- x
  sites.to.bs <- sample(1:nrow(y), size=nrow(y), replace=TRUE)
  y <- y[sites.to.bs, ]
  x <- x[sites.to.bs, ]
  model <- FREEgamboost(y, x, bins, ...)
  class(model) <- "FREEfit"
  coefs <- coef(model)$mean
  if (verbose) {
    cat(paste(100 * {i / n.bs}, "% complete.....", sep=""), "\n")
    flush.console()
  }
  return(coefs)
}