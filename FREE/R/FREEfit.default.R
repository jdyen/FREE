FREEfit.default <-
function(y, x, bins=NULL, groups=NULL, z=NULL,
         method=c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"),
         bootstrap.cis=FALSE, n.bs=50, ...)
{
  if (is.matrix(y)) {
    if (is.null(bins)) {
      bins <- 1:ncol(y)
    }
    if (ncol(y) != length(bins)) {
      stop("Response data y should have the same number of columns as there are bins...",
           call.=FALSE)
    }
    if (!is.null(nrow(x))){
      if (nrow(y) != nrow(x)) {
        stop("Response and predictor data should have one row for each site", call.=FALSE)
      }
    }
    if (!is.null(z)) {
      warning("Matrix z is only used for scalar models...", call.=FALSE)
    }
    if (is.null(groups)) {
      if (!any(method == c("fda", "gamboost", "INLA", "stan",
                           "BUGSspline", "BUGSjump", "FREE"))) {
        warning("Method must be one of fda, gamboost, INLA, stan, BUGSspline, BUGSjump or FREE. fda method selected by default...", call.=FALSE)
        method <- "fda"
      }
      if (method == "fda") {
        model <- FREEfda(y=y, x=x, bins=bins, ...)
      }
      if (method == "gamboost") {
        warning("The FREE package no longer supports the gamboost method; consider using FREE, INLA or BUGSjump methods...", call.=FALSE)
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
        warning("The FREE package no longer supports the stan method; consider using FREE, INLA or BUGSjump methods...", call.=FALSE)
        model <- FREEstan(y=y, x=x, bins=bins, ...)
      }
      if (method == "BUGSspline") {
        if (Sys.info()["sysname"] == "Windows") {
          warning("The FREE package no longer supports the BUGSspline method; consider using FREE, INLA or BUGSjump methods...", call.=FALSE)
          model <- FREEbugs(y=y, x=x, bins=bins, ...)
        } else {
          stop("The BUGSspline method requires the R2WinBUGS package; install this package or consider using the FREE or INLA methods...", call.=FALSE)
        }
      }
      if (method == "BUGSjump") {
        if (Sys.info()["sysname"] == "Windows") {
          model <- FREEbugsJump(y=y, x=x, bins=bins, ...)
        } else {
          stop("The BUGSspline method requires the R2WinBUGS package; install this package or consider using the FREE or INLA methods...", call.=FALSE)
        }
      }
      if (method == "FREE") {
        model <- FREEfree(y=y, x=x, bins=bins, ...)
      }
      model$method <- method
    } else {
      if (method != "FREE") {
        warning("Models with clustering variables are implemented using a Gibbs sampler; fda, gamboost, INLA, stan, BUGSspline and BUGSjump are not available...", call.=FALSE)
      }
      if (is.matrix(groups) | is.integer(groups) | is.numeric(groups)) {
        if (is.integer(groups)) {
          groups <- as.matrix(groups)
        } else {
          if (is.numeric(groups)) {
            warning("A numeric vector has been provided for groups; make sure that this vector is assigning each row of the response to a cluster...", call.=FALSE)
            groups <- as.matrix(groups)
          }
        }
        if (nrow(groups) != nrow(y)) {
          stop("Response and groups should have the same number of rows...", call.=FALSE)
        }
        groups <- apply(groups, 2, function(x) as.integer(as.factor(x)))
        model <- FREEspline(y=y, x=x, groups=groups, w=bins, ...)
        model$bins <- bins
        model$method <- "FREEmixed"
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
          }
          if (is.matrix(groups) | is.integer(groups) | is.numeric(groups)) {
            if (is.integer(groups)) {
              groups <- as.matrix(groups)
            }
            if (is.numeric(groups)) {
              warning("A numeric vector has been provided for groups; make sure that this vector is assigning each observation of the response to a cluster...", call.=FALSE)
              groups <- as.matrix(groups)
            }
            if (nrow(groups) != length(y)) {
              stop("Response and groups should have the same number of rows/observations...", call.=FALSE)
            }
            groups <- apply(groups, 2, function(x) as.integer(as.factor(x)))
          }
          model <- FREEscalar(y=y, x=x, z=z, groups=groups, bins=bins, ...)
          model$bins <- bins
    	  model$method <- "scalar"
    	} else {
    	  stop("Model response must be a numeric vector (scalar response model) or a matrix (function response model)...", call.=FALSE)
    	}
    }
  class(model) <- "FREEfit"
  model
}
