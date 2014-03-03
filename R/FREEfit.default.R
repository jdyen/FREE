FREEfit.default <-
function(y, x, bins=1:ncol(y), method=c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"), ...){
  if (ncol(y) != length(bins)) {
    stop("Response data y should have the same number of columns as there are bins...",
         call.=FALSE)
  }
  if (!is.null(nrow(x))){
    if (nrow(y) != nrow(x)) {
      stop("Response and predictor data should have one row for each site", call.=FALSE)
    }
  }
  if (!any(method == c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"))) {
    warning("Method must be one of fda, gamboost, INLA, stan, BUGSspline, BUGSjump or FREE. fda method selected by default...", call.=FALSE)
    method <- "fda"
  }
  if (method == "fda") {
    model <- FREEfda(y=y, x=x, bins=bins, ...)
  }
  if (method == "gamboost") {
    model <- FREEgamboost(y=y, x=x, bins=bins, ...)
  }
  if (method == "INLA") {
    model <- FREEinla(y=y, x=x, bins=bins, ...)
  }
  if (method == "stan") {
    model <- FREEstan(y=y, x=x, bins=bins, ...)
  }
  if (method == "BUGSspline") {
    model <- FREEbugs(y=y, x=x, bins=bins, ...)
  }
  if (method == "BUGSjump") {
    model <- FREEbugsJump(y=y, x=x, bins=bins, ...)
  }
  if (method == "FREE") {
    model <- FREEfree(y=y, x=x, bins=bins, ...)
  }
  model$method <- method
  class(model) <- "FREEfit"
  model
}
