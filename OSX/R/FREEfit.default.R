FREEfit.default <-
function(y, x, bins=1:ncol(y), method=c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"), bootstrap.cis=FALSE, n.bs=50, ...){
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
    if (bootstrap.cis) {
      coefs <- sapply(1:n.bs, FREEbsInternal, n.bs=n.bs, y=y, x=x, bins=bins,
                      ..., simplify="array")
      coefs.sd <- apply(coefs, c(1,2), sd)
      model$coefs.sd <- coefs.sd
    }
  }
  if (method == "INLA") {
    model <- FREEinla(y=y, x=x, bins=bins, ...)
  }
  if (method == "stan") {
    stop("Method stan not supported for OSX version of FREE...", call.=FALSE)
  }
  if (method == "BUGSspline") {
    stop("Method BUGSspline not supported for OSX version of FREE...", call.=FALSE)
  }
  if (method == "BUGSjump") {
    stop("Method BUGSjump not supported for OSX version of FREE...", call.=FALSE)
  }
  if (method == "FREE") {
    model <- FREEfree(y=y, x=x, bins=bins, ...)
  }
  model$method <- method
  class(model) <- "FREEfit"
  model
}
