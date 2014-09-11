FREEcvInternalScalar <-
function(i, n.obs, n.cv, y, x, bins, verbose, n.vars=n.vars, n.bins=n.bins, ...){
  y.store <- y
  x.store <- x
  inc <- floor(n.obs / n.cv)
  sites.to.cv <- {{i - 1} * inc + 1}:{i * inc}
  y <- y[-sites.to.cv]
  x <- x[-sites.to.cv, ]
  model <- FREEscalar(y, x, bins, n.vars=n.vars, n.bins=n.bins, ...)
  model$method <- "scalar"
  class(model) <- "FREEfit"
  observed <- y.store[sites.to.cv]
  predicted <- predict(model, newdata=x.store[sites.to.cv, ])
  if (verbose) {
    cat(paste(100 * {i / n.cv}, "% complete.....", sep=""), "\n")
    flush.console()
  }
  return(list(observed=observed, predicted=predicted))
}
