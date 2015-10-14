FREEfitCV.default <-
function(y, x, bins=NULL, groups=NULL, z=NULL,
         method=c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"), n.cv=10, verbose=TRUE,
         n.vars=NULL, n.bins=NULL, ...){
  if (is.matrix(y)) {
    if (is.null(bins)) {
      bins <- 1:ncol(y)
    }
  	if (ncol(y) != length(bins)) {
  	  stop("Response data y should have the same number of columns as there are bins...",
  	       call.=FALSE)
  	}
  	if (!is.null(nrow(x))) {
  	  if (nrow(y) != nrow(x)) {
  	    stop("Response and predictor data should have one row for each site...", call.=FALSE)
  	  }
  	}
    if (!is.null(z)) {
      warning("Matrix z is only used for scalar models...", call.=FALSE)
    }
    n.obs <- nrow(y)
    observed <- NULL
    predicted <- NULL
  	if (is.null(groups)) {
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
  	  if (method != "FREE") {
  	    warning("Models with clustering variables are implemented using a Gibbs sampler; fda, gamboost, INLA, stan, BUGSspline and BUGSjump are not available...", call.=FALSE)
  	  }
  	  if (is.matrix(groups) | is.integer(groups) | is.numeric(groups)) {
  	    if (is.integer(groups)) {
  	      groups <- as.matrix(groups)
  	    }
  	    if (is.numeric(groups)) {
  	      warning("A numeric vector has been provided for groups; make sure that this vector is assigning each row of the response to a cluster...", call.=FALSE)
  	      groups <- as.matrix(groups)
  	    }
  	    if (nrow(groups) != nrow(y)) {
  	      stop("Response and groups should have the same number of rows...", call.=FALSE)
  	    }
  	    groups <- apply(groups, 2, function(x) as.integer(as.factor(x)))
  	    cv.out <- FREEsplineCV(n.cv=n.cv, y=y, x=x, groups=groups, w=bins, ...)
  	  } else {
  	    stop("The groups variable needs to be a matrix or integer/numeric vector...", call.=FALSE)
  	  }
  	}
  } else {
    if (is.numeric(y)) {
      if (is.null(bins)) {
        bins <- 1:ncol(x)
      }
      if (!is.matrix(x)) {
        stop("The predictor data should be a matrix with one row for each site...", call.=FALSE)
      }
      if (length(y) != nrow(x)) {
        stop("Response and predictor data should have one row/observation for each site...",
             call.=FALSE)
      }
      if (is.matrix(z) | is.numeric(z) | is.integer(z)) {
        if (is.numeric(z) | is.integer(z)) {
          z <- matrix(z, ncol=1)
        }
        if (nrow(z) != nrow(x)) {
          stop("Response and predictor data should have one row/observation for each site...", call.=FALSE)
        }
      } else {
          stop("Scalar predictors z must be either a matrix or a numeric or integer vector...", call.=FALSE)
      }
      if (is.matrix(groups) | is.integer(groups) | is.numeric(groups)) {
        if (is.integer(groups)) {
          groups <- as.matrix(groups)
        }
        if (is.numeric(groups)) {
          warning("A numeric vector has been provided for groups; make sure that this vector is assigning each row of the response to a cluster...", call.=FALSE)
          groups <- as.matrix(groups)
        }
        if (nrow(groups) != length(y)) {
          stop("Response and groups should have the same number of rows...", call.=FALSE)
        }
        groups <- apply(groups, 2, function(x) as.integer(as.factor(x)))
        cv.out <- FREEscalarCV(n.cv=n.cv, y=y, x=x, z=z, groups=groups, bins=bins, ...)
      } else {
        stop("The groups variable needs to be a matrix or integer/numeric vector...", call.=FALSE)
      }
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
