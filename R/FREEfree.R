FREEfree <-
function(y, x, bins, family="gaussian", errors="iid", model=c("iid", "rw2"), n.chains=3, n.iters=10000, n.burnin=n.iters/5, n.thin=n.iters/100, hypers=list(psi=1, phi=1, sigma2_beta=10), inits=list(vari=10, betas=NULL)) {
  # add intercept to predictors
  x <- cbind(rep(1, nrow(y)), x)
  # extract indices
  n <- nrow(y)
  np <- ncol(x)
  nj <- ncol(y)
  # prepare output objects
  beta.store <- vector(mode="list", length=np)
  for (i in 1:np) {
    beta.store[[i]] <- matrix(NA, nrow={{n.iters - n.burnin} / n.thin}, ncol=nj)
    names(beta.store)[[i]] <- paste("beta", i, sep="")
  }
  sigma2.store <- vector(mode="numeric", length={{n.iters - n.burnin} / n.thin})
  # reformat data
  response <- as.vector(t(y))
  preds <- as.vector(t(x))
  # set hyperparameters
  psis <- hypers$psi
  phis <- hypers$phi
  var_beta <- hypers$sigma2_beta
  # initialise parameters
  if (!is.null(inits$betas)) {
    betas <- inits$betas
  } else {
    betas <- rnorm(np * nj, mean = 0, sd = 1)
  }
  vari <- inits$vari
  # run MCMC loop
  for (i in 1:n.iters) {
    updates <- update_beta(model=model, response=response, preds=preds, n=n, np=np,
                           nj=nj, betas=betas, vari=vari, psis=psis, phis=phis,
                           var_beta=var_beta)
    betas <- updates$betas
    vari <- updates$vari
    # save every n.thin iterations after n.burnin
    if ({i > n.burnin} & {{i %% n.thin} == 0}) {
      for (j in 1:np) {
        beta.store[[j]][{{i - n.burnin} / n.thin}, 1:nj] <- updates$betas[{{j - 1}*nj + 1}:{{j - 1}*nj + nj}]
      }
      sigma2.store[{i - n.burnin} / n.thin] <- vari
    }
  }
  coefs.mean <- matrix(sapply(beta.store, function(x) { return(apply(x, 2, mean)) }), nrow=np, byrow=TRUE)
  coefs.sd <- matrix(sapply(beta.store, function(x) {return(apply(x, 2, sd))}), nrow=np, byrow=TRUE)
  fitted <- matrix(preds, ncol=np, byrow=TRUE) %*% coefs.mean
  xIC <- NULL
  family <- "gaussian"
  y <- matrix(response, ncol=nj, byrow=TRUE)
  r <- cor(c(fitted), c(y))
  r2 <- r * r
  bins <- 1:nj
  return(list(fitted=fitted, observed=y, coefs.mean=coefs.mean, coefs.sd=coefs.sd,
              r2=r2, family=family, bins=bins, xIC=xIC))
}