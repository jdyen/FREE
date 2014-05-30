FREEfitCV.default <-
function(y, x, bins=1:ncol(y), method=c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"), n.cv=10, verbose=TRUE, ...){
  n.obs <- nrow(y)
  observed <- NULL
  predicted <- NULL
  if (!any(method == c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"))) {
    warning("Method must be one of fda, gamboost, INLA, stan, BUGSspline, BUGSjump or FREE. fda method selected by default...", call.=FALSE)
    method <- "fda"
  }
  if (method == "stan") {
    stop("Method stan not supported for OSX version of FREE...", call.=FALSE)
  } else {
  	if (method == "BUGSjump" | method == "BUGSspline") {
	  stop("BUGS methods not supported for OSX version of FREE...", call.=FALSE)
  	} else {
    cv.temp <- lapply(c(1:n.cv), FREEcvInternal, n.obs=n.obs, n.cv=n.cv, y=y, x=x, bins=bins,
                     method=method, verbose=verbose, ...)
    }
  }
  cv.out <- FREEcvOutput(cv.temp)
  observed <- cv.out$observed
  predicted <- cv.out$predicted
  cv.cor <- cor(c(observed), c(predicted))
  cv.r2 <- cv.cor * cv.cor
  CVout <- list(observed=observed, predicted=predicted, cv.r2=cv.r2)
  class(CVout) <- "FREEfitCV"
  CVout
}
