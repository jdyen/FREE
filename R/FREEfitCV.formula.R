FREEfitCV.formula <-
function(formula, data=list(), bins=NULL, method=c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"), n.cv=10, verbose=TRUE, ...){
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  if (length(grep("Intercept", colnames(x))) == 1) {
    x <- x[, -1]
  } else {
    x <- x
  }
  y <- model.response(mf)
  if (is.null(bins)) {
    bins <- 1:ncol(y)
  }
  if (length(method) > 1) {
    method <- "gamboost"
  }
  model <- FREEfitCV.default(y=y, x=x, bins=bins, method=method, verbose=verbose, n.cv=n.cv, ...)
  model$call <- match.call()
  model$formula <- formula
  model
}
