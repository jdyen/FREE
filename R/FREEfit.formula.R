FREEfit.formula <-
function(formula, data=list(), bins=NULL, method=c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"), ...){
  mf <- model.frame(formula=formula, data=data)
  if (is.null(model.response(mf))) {
    stop("Model response must be specified in formula.....", call.=FALSE)
  }
  x <- model.matrix(attr(mf, "terms"), data=mf)
  if (length(grep("Intercept", colnames(x))) == 1) {
     x <- x[, -1]
  } else {
    x <- x
  }
  x <- as.matrix(x)
  if (length(variable.names(mf)) >= 2) {
    if (ncol(x) == {length(variable.names(mf)) - 1}) {
      colnames(x) <- variable.names(mf)[2:{length(variable.names(mf))}]
      colnames(x) <- gsub(":", "", colnames(x))
    }
  }
  y <- model.response(mf)
  if (is.null(bins)) {
    bins <- 1:ncol(y)
  }
  if (length(method) > 1) {
    method <- "gamboost"
  }  
  model <- FREEfit.default(y=y, x=x, bins=bins, method=method, ...)
  model$call <- match.call()
  model$formula <- formula
  model
}
