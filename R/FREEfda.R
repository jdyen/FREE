FREEfda <-
function(y, x, bins, y.basis="bspline", nbasis.y=12, norder.y=6, beta.basis="bspline", nbasis.beta=12, norder.beta=6, loglam=NULL){
  # Fits a basis function model for regression where the response variable
  #   is a function and there are one or more scalar predictor variables
  #
  # Dependencies:
  #   fda
  #   Matrix
  #   lattice
  #
  # Args:
  #   blah
  #
  # Returns
  #   blah
  if (y.basis == "fourier") {
    y.basis <- create.fourier.basis(c(min(bins), max(bins)), nbasis=nbasis.y)
  } else {
    y.basis <- create.bspline.basis(c(min(bins), max(bins)), nbasis=nbasis.y,
                                    norder=norder.y)
  }
  basismat.y <- eval.basis(bins, y.basis)
  y.fd <- Data2fd(argvals=bins, y=t(y), basisobj=y.basis)
  if (beta.basis == "fourier") {
    beta.basis <- create.fourier.basis(c(min(bins), max(bins)), nbasis=nbasis.beta)
  } else {
    beta.basis <- create.bspline.basis(c(min(bins), max(bins)), nbasis=nbasis.beta,
                                       norder=norder.beta)
  }
  basismat.beta <- eval.basis(bins, beta.basis)
  beta.fdPar <- fdPar(beta.basis)
  beta.list <- vector(mode="list", length={ncol(x) + 1})
  for (i in seq(along=beta.list)) {
    beta.list[[i]] <- beta.fdPar
    if (!is.null(loglam)) {
      beta.list[[i]]$lambda <- 10 ^ loglam
    }
  }
  data.list <- list()
  data.list[[1]] <- rep(1, nrow(x))
  for (i in 1:ncol(x)) {
    data.list[[{i + 1}]] <- x[, i]
  }
  names(data.list) <- c("int", colnames(x))
  fRegress.fit <- fRegress(y.fd, xfdlist=data.list, betalist=beta.list)
  fitted.coefs <- fRegress.fit$yhatfd$fd$coefs
  y.coefs <- y.fd$coefs
  beta.coefs <- lapply(fRegress.fit$betaestlist, function(x) return(x$fd))
  fitted.betas <- basismat.beta %*% sapply(beta.coefs, function(x) x$coef)
  fitted.vals <- basismat.y %*% fRegress.fit$yhatfd$fd$coefs
  yhatmat <- eval.fd(bins, fRegress.fit$yhatfdobj$fd)
  rmatb <- t(y) - yhatmat
  SigmaEb <- var(t(rmatb))
  cmap.temp <- smooth.basis(bins, t(y), y.fd)
  y2cMap <- cmap.temp$y2cMap
  fda.sd <- fRegress.stderr(fRegress.fit, y2cMap, SigmaEb)
  beta.sd.est <- fda.sd$betastderrlist
  fitted.sd <- basismat.beta %*% sapply(beta.sd.est, function(x) x$coef)
  colnames(fitted.betas) <- c("INTERCEPT", colnames(x))
  colnames(fitted.sd) <- colnames(fitted.betas)
  r <- cor(c(t(fitted.vals)), c(y))
  r2 <- r * r
  family <- "Gaussian"
  xIC <- NULL
  return(list(fitted=t(fitted.vals), observed=y, coefs.mean=t(fitted.betas), 
              coefs.sd=t(fitted.sd), r2=r2, family=family, bins=bins, xIC=xIC))
}
