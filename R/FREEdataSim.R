FREEdataSim <-
function(n=100, n.vars=2, bins=NULL, y.mean=NULL, x.val=NULL, beta.mean=NULL, error.mean=0, error.sd=1, cov.mean=NULL, cov.sd=NULL){
  library(MASS)
  if (!is.null(x.val)) {
    if (n.vars != ncol(x.val)) {
      warning("x.val does not define n.vars variables: some variables set to default values.....", call.=FALSE)
    }
  }
  if (!is.null(beta.mean)) {
  	if (ncol(matrix(beta.mean)) == 1) {
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
  max.noise <- max(y.mean) / 4
  y.noise.ar1 <- matrix(unlist(lapply(1:n, function(x, n.sim=n.bins) arima.sim(list(ar=runif(1, 0.5, 0.85)), n=n.sim)[1:n.sim])), nrow=n, byrow=TRUE)
  y.noise.mvn <- mvrnorm(n=n, mu=rep(0, n.bins), Sigma=MakeCovMat(n=n.bins))
  y.noise.iid <- max.noise * sweep(y.noise.iid, c(1, 2), apply(y.noise.iid, 1,
                   function(x) max(abs(x))), "/")
  y.noise.ar1 <- max.noise * sweep(y.noise.ar1, c(1, 2), apply(y.noise.ar1, 1,
                   function(x) max(abs(x))), "/")
  y.noise.mvn <- max.noise * sweep(y.noise.mvn, c(1, 2), apply(y.noise.mvn, 1,
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
  y.beta <- cov.val %*% y.cov
  y.obs.iid <- y.obs.iid + y.beta
  y.obs.ar1 <- y.obs.ar1 + y.beta
  y.obs.mvn <- y.obs.mvn + y.beta
  cov.val <- data.frame(cov.val)
  for (i in 1:ncol(cov.val)) {
    colnames(cov.val)[i] <- paste("VAR", i, sep="")
  }
  return(list(y.iid=y.obs.iid, y.ar1=y.obs.ar1, y.mvn=y.obs.mvn, x=cov.val, mean.real=y.mean, beta.real=y.cov, error.iid=y.noise.iid, error.ar1=y.noise.ar1, error.mvn=y.noise.mvn, bins=x))
}
