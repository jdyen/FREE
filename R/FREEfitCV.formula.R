FREEfitCV.formula <-
function(formula, data=list(), bins=NULL, groups=NULL, z=NULL, n.cv=10,
         verbose=TRUE, n.iters=100, n.burnin=round(n.iters / 5),
         n.chains=3, par.run=FALSE, ...){
  mf <- model.frame(formula=formula, data=data)
  if (is.null(model.response(mf))) {
    stop("model response must be specified in formula", call.=FALSE)
  }
  x <- model.matrix(attr(mf, "terms"), data=mf)
  if ((length(grep("Intercept", colnames(x))) == 1) & (ncol(x) > 1)) {
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
  	                           verbose=verbose, n.cv=n.cv, n.iters=n.iters,
  	                           n.burnin=n.burnin, n.chains=n.chains,
  	                           par.run=par.run, ...)
  } else {
    if (is.numeric(model.response(mf))) {
      y <- model.response(mf)
      if (length(names(mf)) > 1) {
        x.use <- vector('list', length = (length(names(mf)) - 1))
        for (i in 2:length(names(mf))) {
          x.use[[i - 1]] <- as.matrix(mf[names(mf)[i]])
        }
      } else {
        x.use <- NULL
      }
      x <- x.use
      if (is.null(bins) & !is.null(x)) {
        bins <- vector('list', length = length(x))
        for (i in seq(along = x)) {
          bins[[i]] <- 1:ncol(x[[i]])
        }
      }
      model <- FREEfitCV.default(y=y, x=x, bins=bins, groups=groups, z=z,
                                 verbose=verbose, n.cv=n.cv, n.iters=n.iters,
                                 n.burnin=n.burnin, n.chains=n.chains,
                                 par.run=par.run, ...)
  	} else {
  	  stop("model response must be a numeric vector (scalar response model) or a matrix (function response model)", call.=FALSE)
  	}
  }
  model$call <- match.call()
  model$formula <- formula
  model
}
