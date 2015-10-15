FREEfitCV.formula <-
function(formula, data=list(), bins=NULL, groups=NULL, z=NULL, n.cv=10, verbose=TRUE, ...){
  mf <- model.frame(formula=formula, data=data)
  if (is.null(model.response(mf))) {
    stop("model response must be specified in formula", call.=FALSE)
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
  	if (is.null(bins)) {
  	  bins <- 1:ncol(y)
  	}
  	model <- FREEfitCV.default(y=y, x=x, bins=bins, groups=groups, z=z,
  	                           verbose=verbose, n.cv=n.cv, ...)
  } else {
  	if (is.numeric(model.response(mf))) {
  	  y <- model.response(mf)
  	  n.bins <- ncol(x) / (length(variable.names(mf)) - 1)
  	  n.vars <- length(variable.names(mf)) - 1
  	  if (is.null(bins)) {
  	    bins <- 1:n.bins
  	  }
  	  model <- FREEfitCV.default(y=y, x=x, bins=bins, groups=groups, z=z,
  	                             verbose=verbose, n.cv=n.cv, ...)
  	} else {
  	  stop("model response must be a numeric vector (scalar response model) or a matrix (function response model)", call.=FALSE)
  	}
  }
  model$call <- match.call()
  model$formula <- formula
  model
}
