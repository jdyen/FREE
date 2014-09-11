FREEdataSim <-
function(function.class=c("response", "predictor"), n=100, n.vars=2, bins=NULL, y.mean=NULL, x.val=NULL, beta.mean=NULL, error.mean=0, error.sd=1, cov.mean=NULL, cov.sd=NULL, struc.err.mean=0, struc.err.sd=1){
  if (length(function.class) == 2) {
    warning("function.class has not been specified and a function response model is used as default.", call.=FALSE)
    function.class <- "response"
  }
  if (function.class == "response") {
    if (!is.null(x.val)) {
      if (n.vars != ncol(x.val)) {
        warning("x.val does not define n.vars variables: some variables set to default values.....", call.=FALSE)
      }
    }
    if (!is.null(beta.mean)) {
      if (ncol(as.matrix(beta.mean)) == 1) {
  	    beta.mean <- t(matrix(beta.mean))
  	  } else {
  	    beta.mean <- matrix(beta.mean)
      }
    }
    if (!is.null(beta.mean)) {
      if (n.vars != nrow(beta.mean)) {
        warning("beta.mean does not define n.vars variables: some beta set to default values.....", call.=FALSE)
      }
    }
    if (!is.null(beta.mean)) {
      if (nrow(beta.mean) > n.vars) {
        stop("beta.mean should have no more than n.vars rows.....", call.=FALSE)
      }
    }
    if (!is.null(x.val)) {
      if (is.null(n)) {
  	    n <- length(x.val)
  	  } else {
        if (length(x.val) != n){
          stop("x.val must be the same length as n.....", call.=FALSE)
        }
      }
    }
    if ({!is.null(cov.mean) & length(cov.mean) != n.vars} | 
        {!is.null(cov.sd) & length(cov.sd) != n.vars}) {
      stop("cov.mean and cov.sd must have n.vars values.....", call.=FALSE)
    }
    if (!is.null(y.mean) & !is.null(beta.mean) & is.null(bins)) {
      if (length(y.mean) != length(beta.mean)) {
        stop("y.mean and beta.mean must have the same length.....", call.=FALSE)
      } else {
        bins <- length(y.mean)
      }
    }
    if (!is.null(y.mean) & is.null(beta.mean) & is.null(bins)) {
      bins <- 1:length(y.mean)
    }
    if (is.null(y.mean) & !is.null(beta.mean) & is.null(bins)) {
      bins <- 1:ncol(beta.mean)
    }
    if (!is.null(bins)) {
      if (!is.null(y.mean)) {
        if (length(y.mean) != length(bins)) {
          stop("y.mean and bins must specify the same length vector.....", call.=FALSE)
        }
      }
      if (!is.null(beta.mean)) {
        if (ncol(beta.mean) != length(bins)) {
          stop("beta.mean and bins must specify the same length vector.....", call.=FALSE)
        }
      }
      if (length(bins) == 1) {
        x <- seq(-10, 10, length=bins)
      } else {
        x <- bins
      }
    } else {
      x <- seq(-10, 10, by=1)
    }
    n.bins <- length(x)
    if (!is.null(y.mean)) {
      y.mean <- y.mean
    } else {
      y.mean <- sin(0.25 * x + 1.6)
    }
    y.noise.iid <- matrix(rnorm(length(x) * n, mean=error.mean, sd=error.sd), nrow=n)
    y.noise.ar1 <- matrix(unlist(lapply(1:n, function(x, n.sim=n.bins) arima.sim(list(ar=runif(1, 0.5, 0.85)), n=n.sim)[1:n.sim])), nrow=n, byrow=TRUE)
    y.noise.mvn <- mvrnorm(n=n, mu=rep(0, n.bins), Sigma=MakeCovMat(n=n.bins))
    y.noise.iid <- sweep(y.noise.iid, c(1, 2), apply(y.noise.iid, 1,
                     function(x) max(abs(x))), "/")
    y.noise.ar1 <- sweep(y.noise.ar1, c(1, 2), apply(y.noise.ar1, 1,
                     function(x) max(abs(x))), "/")
    y.noise.mvn <- sweep(y.noise.mvn, c(1, 2), apply(y.noise.mvn, 1,
                     function(x) max(abs(x))), "/")
    cov.val <- matrix(NA, nrow=n, ncol=n.vars)
    if (is.null(cov.mean)) {
      cov.mean <- sample(1:10, size=n.vars, replace=TRUE)
    }
    if (is.null(cov.sd)) {
      cov.sd <- sample(c(0.5, 1, 1.5, 2), size=n.vars, replace=TRUE)
    }
    if (!is.null(x.val)) {
      x.val <- matrix(x.val)
  	  if (ncol(x.val) > 1) {
        for (i in 1:ncol(x.val)) {
          cov.val[, i] <- x.val[, i]
        }
      } else {
        cov.val[, 1] <- x.val
      }
      if (ncol(x.val) != n.vars) {
        n.cov.missing <- n.vars - ncol(x.val)
        for (i in 1:n.cov.missing) {
          cov.val[, {ncol(x.val) + i}] <- rnorm(n, mean=cov.mean[i], sd=cov.sd[i])
        }
      }
    } else {
      for (i in 1:n.vars) {
        cov.val[, i] <- rnorm(n, mean=cov.mean[i], sd=cov.sd[i])
      }
    }
    cov.val <- apply(cov.val, 2, function(x) x / max(x))
    y.cov <- matrix(NA, ncol=length(x), nrow=n.vars)
    if (!is.null(beta.mean)) {
      for (i in 1:nrow(beta.mean)) {
        y.cov[i, ] <- beta.mean[i, ]
      }
      if (nrow(beta.mean) != n.vars) {
        n.beta.missing <- n.vars - nrow(beta.mean)
        for (i in 1:n.beta.missing) {
          y.cov[{nrow(beta.mean) + i}, ] <- runif(1, 0.1, 2) * sin(runif(1, 0.1, 1) * x + i * {pi / 3})
        }
      }
    } else {
      for (i in 1:n.vars) {
        y.cov[i, ] <- runif(1, 0.1, 2) * sin(runif(1, 0.1, 1) * x + i * {pi / 3})
      }
    }
    y.obs.iid <- sweep(y.noise.iid, 2, y.mean, "+")
    y.obs.ar1 <- sweep(y.noise.ar1, 2, y.mean, "+")
    y.obs.mvn <- sweep(y.noise.mvn, 2, y.mean, "+")
    struc.err <- rnorm(nrow(y.obs.iid), mean=struc.err.mean, sd=struc.err.sd)
    y.obs.iid <- sweep(y.obs.iid, c(1, 2), struc.err, "+")
    y.obs.ar1 <- sweep(y.obs.ar1, c(1, 2), struc.err, "+")
    y.obs.mvn <- sweep(y.obs.mvn, c(1, 2), struc.err, "+")
    y.beta <- cov.val %*% y.cov
    y.obs.iid <- y.obs.iid + y.beta
    y.obs.ar1 <- y.obs.ar1 + y.beta
    y.obs.mvn <- y.obs.mvn + y.beta
    cov.val <- data.frame(cov.val)
    for (i in 1:ncol(cov.val)) {
      colnames(cov.val)[i] <- paste("VAR", i, sep="")
    }
  } else {
    if (function.class == "predictor") {
      if (!is.null(x.val)) {
        if ((n.vars > 1) & (n.vars != dim(x.val)[3])) {
          warning("x.val does not define n.vars variables: some variables set to default values.....", call.=FALSE)
        }
      }
      if (!is.null(beta.mean)) {
        if (ncol(as.matrix(beta.mean)) == 1) {
    	  beta.mean <- t(matrix(beta.mean))
  	    } else {
  	      beta.mean <- as.matrix(beta.mean)
        }
      }
      if (!is.null(beta.mean)) {
        if (n.vars != nrow(beta.mean)) {
          warning("beta.mean does not define n.vars variables: some beta set to default values.....", call.=FALSE)
        }
      }
      if (!is.null(beta.mean)) {
        if (nrow(beta.mean) > n.vars) {
          stop("beta.mean should have no more than n.vars rows.....", call.=FALSE)
        }
      }
      if (!is.null(x.val)) {
        if (is.null(n)) {
    	  n <- nrow(x.val)
  	    } else {
          if (nrow(x.val) != n){
            stop("x.val must have n rows.....", call.=FALSE)
          }
        }
      }
      if ({!is.null(cov.mean) & length(cov.mean) != n.vars} | 
          {!is.null(cov.sd) & length(cov.sd) != n.vars}) {
        stop("cov.mean and cov.sd must have n.vars values.....", call.=FALSE)
      }
      if (!is.null(x.val) & !is.null(beta.mean) & is.null(bins)) {
      	if (ncol(x.val) != ncol(beta.mean)) {
      	  stop("x.val and beta.mean must have the same number of columns/bins.....",
      	       call.=FALSE)
      	} else {
          bins <- 1:ncol(beta.mean)
        }
      }
      if (!is.null(x.val) & is.null(beta.mean) & is.null(bins)) {
        bins <- 1:ncol(x.val)
      }
      if (is.null(x.val) & !is.null(beta.mean) & is.null(bins)) {
        bins <- 1:ncol(beta.mean)
      }
      if (!is.null(bins)) {
        if (length(bins) == 1) {
          bins <- seq(-10, 10, length=bins)
        }
        if (!is.null(beta.mean)) {
          if (ncol(beta.mean) != length(bins)) {
            stop("beta.mean and bins must specify the same length vector.....", call.=FALSE)
          }
        }
        x <- bins
      } else {
        x <- seq(-10, 10, by=1)
      }
      n.bins <- length(x)
      if (!is.null(y.mean)) {
      	if (length(y.mean) > 1) {
      	  warning("y.mean has length > 1; only the first value will be used.....", call.=FALSE)
      	  y.mean <- y.mean[1]
      	} else {
          y.mean <- y.mean
        }
      } else {
        y.mean <- rnorm(1, sd=10)
      }
      y.noise.iid <- rnorm(n, mean=error.mean, sd=error.sd)
      y.noise.ar1 <- NULL
      y.noise.mvn <- NULL
      cov.val <- array(NA, dim=c(n, n.bins, n.vars))
      if (is.null(cov.mean)) {
        cov.mean <- matrix(sample(1:10, size=(n.vars * n.bins), replace=TRUE),
                           ncol=n.bins)
      }
      if (is.null(cov.sd)) {
        cov.sd <- matrix(sample(c(0.5, 1, 1.5, 2), size=(n.vars * n.bins), replace=TRUE),
                         ncol=n.bins)
      }
      if (!is.null(x.val)) {
        if ((n.vars > 1) & (!is.na(dim(x.val)[3]))) {
          for (i in 1:dim(x.val)[3]) {
            cov.val[, , i] <- x.val[, , i]
          }
        } else {
          cov.val[, , 1] <- x.val
        }
       	if (is.na(dim(x.val)[3])) {
          n.vars.given <- 1
        } else {
          n.vars.given <- dim(x.val)[3]
        }
        if (n.vars.given != n.vars) {
          n.cov.missing <- n.vars - n.vars.given
          for (i in 1:n.cov.missing) {
            cov.val[, , {n.vars.given + i}] <- matrix(rnorm(n * n.bins,
                       mean=cov.mean[i, ], sd=cov.sd[i, ]), ncol=n.bins)
          }
        }
      } else {
        for (i in 1:n.vars) {
          cov.val[, , i] <- matrix(rnorm(n * n.bins, mean=cov.mean[i, ], sd=cov.sd[i, ]),
                                   ncol=n.bins)
        }
      }
      y.cov <- matrix(NA, ncol=length(x), nrow=n.vars)
      if (!is.null(beta.mean)) {
        for (i in 1:nrow(beta.mean)) {
          y.cov[i, ] <- beta.mean[i, ]
        }
        if (nrow(beta.mean) != n.vars) {
          n.beta.missing <- n.vars - nrow(beta.mean)
          for (i in 1:n.beta.missing) {
            y.cov[{nrow(beta.mean) + i}, ] <- runif(1, 0.1, 2) * sin(runif(1, 0.1, 1) * x + i * {pi / 3})
          }
        }
      } else {
        for (i in 1:n.vars) {
          y.cov[i, ] <- runif(1, 0.1, 2) * sin(runif(1, 0.1, 1) * x + i * {pi / 3})
        }
      }
      y.obs.iid <- y.mean + y.noise.iid
      y.obs.ar1 <- NULL
      y.obs.mvn <- NULL
      struc.err <- rep(0, n)
      y.obs.iid <- y.obs.iid + struc.err
      y.obs.ar1 <- NULL
      y.obs.mvn <- NULL
      y.beta <- rep(0, n)
      for (i in 1:n.vars) {
        y.beta <- y.beta + apply(sweep(cov.val[, , i], 2, y.cov[i, ], "*"), 1, sum)
      }
      y.obs.iid <- y.obs.iid + y.beta
      y.obs.ar1 <- NULL
      y.obs.mvn <- NULL
    } else {
      stop("function.class must be set to either predictor or response.", call.=FALSE)
    }
  }
  return(list(y.iid=y.obs.iid, y.ar1=y.obs.ar1, y.mvn=y.obs.mvn, x=cov.val, mean.real=y.mean, beta.real=y.cov, error.iid=y.noise.iid, error.ar1=y.noise.ar1, error.mvn=y.noise.mvn, bins=x))
}
