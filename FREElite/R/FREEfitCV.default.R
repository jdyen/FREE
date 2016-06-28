FREEfitCV.default <-
function(y, x, bins=NULL, groups=NULL, z=NULL, n.cv=10, verbose=TRUE, n.iters=100,
         n.burnin=round(n.iters / 5), n.chains=3, par.run=FALSE, ...){
  if (is.matrix(y)) {
    if (is.null(bins)) {
      bins <- 1:ncol(y)
    }
  	if (ncol(y) != length(bins)) {
  	  stop("response data y should have the same number of columns as there are bins",
  	       call.=FALSE)
  	}
  	if (!is.null(nrow(x))) {
  	  if (nrow(y) != nrow(x)) {
  	    stop("response and predictor data should have one row for each site", call.=FALSE)
  	  }
  	}
    if (!is.null(z)) {
      warning("matrix z is only used for scalar models", call.=FALSE)
    }
    n.obs <- nrow(y)
    observed <- NULL
    predicted <- NULL
  	if (is.null(groups)) {
      cv.out <- FREEsplineCV(n.cv=n.cv, y=y, x=x, groups=groups, bins=bins, ...)
  	} else {
  	  if (is.matrix(groups) | is.integer(groups) | is.numeric(groups)) {
  	    if (is.integer(groups)) {
  	      groups <- as.matrix(groups)
  	    }
  	    if (is.numeric(groups)) {
  	      warning("a numeric vector has been provided for groups; make sure that this vector is assigning each row of the response to a cluster", call.=FALSE)
  	      groups <- as.matrix(groups)
  	    }
  	    if (nrow(groups) != nrow(y)) {
  	      stop("response and groups should have the same number of rows", call.=FALSE)
  	    }
  	    groups <- apply(groups, 2, function(x) as.integer(as.factor(x)))
  	    cv.out <- FREEsplineCV(n.cv=n.cv, y=y, x=x, groups=groups, bins=bins,
                               n.iters=n.iters, n.burnin=n.burnin, n.chains=n.chains,
                               par.run=par.run, verbose=verbose, ...)
  	  } else {
  	    stop("the groups variable needs to be a matrix or integer/numeric vector", call.=FALSE)
  	  }
  	}
  } else {
    if (is.numeric(y)) {
      if (is.null(bins)) {
       bins <- vector('list', length = length(x))
        for (i in seq(along = x)) {
          bins[[i]] <- 1:ncol(x[[i]])
        }
      }
      if (!all(sapply(x, is.matrix))) {
        stop("the predictor data should be a matrix with one row for each site", call.=FALSE)
      }
      if (!all(sapply(x, nrow) == length(y))) {
        stop("response and predictor data should have one row/observation for each site",
             call.=FALSE)
      }
      if (is.null(z)) {
        z <- matrix(0, nrow = length(y), ncol = 2)
      }
      if (is.matrix(z) | is.numeric(z) | is.integer(z)) {
        if (!is.matrix(z)) {
          z <- matrix(z, ncol=1)
        }
        if (!all(nrow(z) == sapply(x, nrow))) {
          stop("response and predictor data should have one row/observation for each site", call.=FALSE)
        }
      } else {
          stop("scalar predictors z must be either a matrix or a numeric or integer vector", call.=FALSE)
      }
      if (is.null(groups)) {
        groups <- matrix(rep(c(1, 2), times = length(y)), ncol = 2)
        for (i in 1:ncol(groups)) {
          groups[, i] <- as.integer(groups[, i])
        }
      }
      if (is.matrix(groups) | is.integer(groups) | is.numeric(groups)) {
        if (!is.matrix(groups)) {
          groups <- as.matrix(groups)
        }
        if (!all(apply(groups, 2, function(x) all((x %% 1) == 0)))) {
          warning("groups are not integers; make sure that this vector is assigning each row of the response to a cluster",
                  call. = FALSE)
          groups <- as.matrix(groups)
        }
        if (nrow(groups) != length(y)) {
          stop("response and groups should have the same number of rows", call.=FALSE)
        }
        groups <- apply(groups, 2, function(x) as.integer(as.factor(x)))
        cv.out <- FREEscalarCV(n.cv=n.cv, y=y, x=x, z=z, groups=groups, bins=bins,
                               n.iters=n.iters, n.burnin=n.burnin, n.chains=n.chains,
                               par.run=par.run, verbose=verbose, ...)
      } else {
        stop("the groups variable needs to be a matrix or integer/numeric vector", call.=FALSE)
      }
    } else {
      stop("model response must be a numeric vector (scalar response model) or a matrix (function response model)", call.=FALSE)
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