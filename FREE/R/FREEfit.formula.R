FREEfit.formula <-
function(formula, data=list(), bins=NULL, groups=NULL, z=NULL,
         method=c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"), ...){
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
  if (is.matrix(model.response(mf))) {
    if (length(variable.names(mf)) >= 2) {
      if (ncol(x) == {length(variable.names(mf)) - 1}) {
        colnames(x) <- variable.names(mf)[2:{length(variable.names(mf))}]
        colnames(x) <- gsub(":", "", colnames(x))
      }
    }
    y <- model.response(mf)
    if (length(method) > 1) {
      if (is.null(groups)) {
        warning("Only one method should be supplied; fda method used by default...", call.=FALSE)
        method <- "fda"
      } else {
        warning("Only one method should be supplied; FREE method used as default for model with clustering variables...", call.=FALSE)
        method <- "FREE"
      }
    }
    model <- FREEfit.default(y=y, x=x, bins=bins, groups=groups, z=z, method=method, ...)
  } else {
    if (is.numeric(model.response(mf))) {
      y <- model.response(mf)
      if (is.null(bins)) {
        bins <- 1:ncol(x)
      }
      model <- FREEfit.default(y=y, x=x, bins=bins, groups=groups, z=z,
                               method=method, ...)
    } else {
      stop("Model response must be a numeric vector (scalar response model) or a matrix (function response model).", call.=FALSE)
    }
  }  
  model$call <- match.call()
  model$formula <- formula
  model
}
