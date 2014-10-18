FREEcvInternal <-
function(i, n.obs, n.cv, y, x, bins, method, verbose, stan.model=NA, ...){
  y.store <- y
  x.store <- x
  inc <- floor(n.obs / n.cv)
  sites.to.cv <- {{i - 1} * inc + 1}:{i * inc}
  y <- y[-sites.to.cv, ]
  x <- x[-sites.to.cv, ]
  if (method == "fda") {
    model <- FREEfda(y, x, bins, ...)
    class(model) <- "FREEfit"
    model$method <- method
  }
  if (method == "gamboost") {
    model <- FREEgamboost(y, x, bins, ...)
    class(model) <- "FREEfit"
    model$method <- method
  }
  if (method == "INLA") {
    model <- FREEinla(y, x, bins, ...)
    class(model) <- "FREEfit"
    model$method <- method
  }
  if (method == "stan") {
    if (i == 1) {
      model <- FREEstan(y, x, bins, refresh=0, ...)
      stan.model <- model$stan.model
    } else {
      model <- FREEstan(y, x, bins, stan.model=stan.model, refresh=0, ...)
    }
    class(model) <- "FREEfit"
    model$method <- method
  }
  if (method == "BUGSspline") {
    model <- FREEbugs(y, x, bins, ...)
    class(model) <- "FREEfit"
    model$method <- method
  }
  if (method == "BUGSjump") {
    model <- FREEbugsJump(y, x, bins, ...)
    class(model) <- "FREEfit"
    model$method <- method
  }
  if (method == "FREE") {
    model <- FREEfree(y, x, bins, ...)
    class(model) <- "FREEfit"
    model$method <- method
  }
  observed <- y.store[sites.to.cv, ]
  predicted <- predict(model, newdata=x.store[sites.to.cv, ])
  if (verbose) {
    cat(paste(100 * {i / n.cv}, "% complete.....", sep=""), "\n")
    flush.console()
  }
  if ({method == "stan"} & {i == 1}) {
    return(list(observed=observed, predicted=predicted, stan.model=stan.model))
  } else {
    return(list(observed=observed, predicted=predicted))
  }
}
