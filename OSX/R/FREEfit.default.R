FREEfit.default <-
function(y, x, bins=1:ncol(y), method=c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"), n.vars=NULL, n.bins=NULL, bootstrap.cis=FALSE, n.bs=50, ...){
  if (is.matrix(y)) {
    if (ncol(y) != length(bins)) {
      stop("Response data y should have the same number of columns as there are bins...",
           call.=FALSE)
    }
    if (!is.null(nrow(x))){
      if (nrow(y) != nrow(x)) {
        stop("Response and predictor data should have one row for each site", call.=FALSE)
      }
    }
    if (!any(method == c("fda", "gamboost", "INLA", "stan",
                         "BUGSspline", "BUGSjump", "FREE"))) {
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
  } else {
  	if (is.numeric(y)) {
  	  if (is.null(n.vars) | is.null(n.bins)) {
  	    stop("The number of variables and number of bins must be specified.", call.=FALSE)
  	  }
  	  if (n.bins != length(bins)) {
  	    stop("Predictor data x should have the same number of columns as there are bins...",
  	         call.=FALSE)
  	  }
  	  if (!is.null(nrow(x))) {
  	    if (length(y) != nrow(x)) {
  	      stop("Response and predictor data should have one row/observation for each site.",
  	           call.=FALSE)
  	    }
  	  }
  	  model <- FREEscalar(y=y, x=x, bins=bins, n.vars=n.vars, n.bins=n.bins, ...)
  	  model$method <- "scalar"
  	} else {
  	  stop("Model response must be a numeric vector (scalar response model) or a matrix (function response model).", call.=FALSE)
  	}
  }
    class(model) <- "FREEfit"
    model
}
