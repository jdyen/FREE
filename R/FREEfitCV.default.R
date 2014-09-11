FREEfitCV.default <-
function(y, x, bins=1:ncol(y), method=c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"), n.cv=10, verbose=TRUE, n.vars=NULL, n.bins=NULL, ...){
  if (is.matrix(y)) {
  	if (ncol(y) != length(bins)) {
  	  stop("Response data y should have the same number of columns as there are bins...",
  	       call.=FALSE)
  	}
  	if (!is.null(nrow(x))) {
  	  if (nrow(y) != nrow(x)) {
  	    stop("Response and predictor data should have one row for each site", call.=FALSE)
  	  }
  	}
    n.obs <- nrow(y)
    observed <- NULL
    predicted <- NULL
    if (!any(method == c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"))) {
      warning("Method must be one of fda, gamboost, INLA, stan, BUGSspline, BUGSjump or FREE. fda method selected by default...", call.=FALSE)
      method <- "fda"
    }
    if (method == "stan") {
      cv.temp <- list()
      cv.temp.stan <- FREEcvInternal(i=1, n.obs=n.obs, n.cv=n.cv, y=y, x=x, bins=bins,
                                     method=method, verbose=verbose, ...)
      cv.temp[[1]] <- list(observed=cv.temp.stan$observed, predicted=cv.temp.stan$predicted)
      stan.model <- cv.temp.stan$stan.model
      cv.temp2 <- lapply(c(2:n.cv), FREEcvInternal, n.obs=n.obs, n.cv=n.cv, y=y, x=x, bins=bins,
                         method=method, verbose=verbose, stan.model=stan.model, ...)
      cv.temp <- c(cv.temp, cv.temp2)
    } else {
      cv.temp <- lapply(c(1:n.cv), FREEcvInternal, n.obs=n.obs, n.cv=n.cv, y=y, x=x, bins=bins,
                       method=method, verbose=verbose, ...)
    }
    cv.out <- FREEcvOutput(cv.temp)
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
      n.obs <- length(y)
      observed <- NULL
      predicted <- NULL
      cv.temp <- lapply(c(1:n.cv), FREEcvInternalScalar, n.obs=n.obs, n.cv=n.cv, y=y,
                        x=x, bins=bins, verbose=verbose, n.vars=n.vars, n.bins=n.bins,
                        ...)
      predicted <- sapply(cv.temp, function(x) x$predicted)
      observed <- sapply(cv.temp, function(x) x$observed)
      cv.out <- list(predicted=predicted, observed=observed)
    } else {
      stop("Model response must be a numeric vector (scalar response model) or a matrix (function response model).", call.=FALSE)
    }
  }
  observed <- cv.out$observed
  predicted <- cv.out$predicted
  cv.cor <- cor(c(observed), c(predicted))
  cv.r2 <- cv.cor * cv.cor
  CVout <- list(observed=observed, predicted=predicted, cv.r2=cv.r2)
  class(CVout) <- "FREEfitCV"
  CVout
}
