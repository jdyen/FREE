FREEscalar <-
function(y, x, bins, n.vars, n.bins, family="gaussian", errors="iid",
         n.chains=3, n.iters=10000, n.burnin=n.iters/5, n.thin=n.iters/100,
         hypers=list(psi=1, phi=100, sigma2_alpha=10, sigma2_beta=10),
         inits=list(vari=10, alpha=NULL, betas=NULL)) {
  # extract indices
  n <- length(y)
  np <- n.vars
  nj <- n.bins
  # reformat data
  response <- y
  preds <- c(t(x))
  # set hyperparameters
  psis <- hypers$psi
  phis <- hypers$phi
  var_alpha <- hypers$sigma2_alpha
  var_beta <- hypers$sigma2_beta
  # initialise chains
  alpha.mean.temp <- vector(mode="numeric", length=n.chains)
  alpha.sd.temp <- vector(mode="numeric", length=n.chains)
  coefs.mean.temp <- array(0, dim=c(np, nj, n.chains))
  coefs.sd.temp <- array(0, dim=c(np, nj, n.chains))
  fitted.temp <- array(0, dim=c(n, n.chains))
  for (w in 1:n.chains) {
    # prepare output objects
    alpha.store <- vector(mode="numeric", length={{n.iters - n.burnin} / n.thin})
    beta.store <- vector(mode="list", length=np)
    for (i in 1:np) {
      beta.store[[i]] <- matrix(NA, nrow={{n.iters - n.burnin} / n.thin}, ncol=nj)
      names(beta.store)[[i]] <- paste("beta", i, sep="")
    }
    sigma2.store <- vector(mode="numeric", length={{n.iters - n.burnin} / n.thin})
    # initialise parameters
    if (!is.null(inits$alpha)) {
      alpha <- inits$alpha
    } else {
      alpha <- rnorm(1, mean=0, sd=1)
    }
    if (!is.null(inits$betas)) {
      betas <- inits$betas
    } else {
      betas <- rnorm(np * nj, mean = 0, sd = 1)
    }
    vari <- inits$vari
    # run MCMC loop
    for (i in 1:n.iters) {
      updates <- update_beta_scalar(response=response, preds=preds, n=n, np=np,
                             nj=nj, alpha=alpha, betas=betas, vari=vari, psis=psis, phis=phis,
                             var_beta=var_beta, var_alpha=var_alpha)
      alpha <- updates$alpha
      betas <- updates$betas
      vari <- updates$vari
      # save every n.thin iterations after n.burnin
      if ({i > n.burnin} & {{i %% n.thin} == 0}) {
      	alpha.store[{i - n.burnin} / n.thin] <- alpha
        for (j in 1:np) {
          beta.store[[j]][{{i - n.burnin} / n.thin}, 1:nj] <-
                    updates$betas[{{j - 1}*nj + 1}:{{j - 1}*nj + nj}]
        }
        sigma2.store[{i - n.burnin} / n.thin] <- vari
      }
    }
    alpha.mean.temp[w] <- mean(alpha.store)
    alpha.sd.temp[w] <- sd(alpha.store)
    coefs.mean.temp[, , w] <- matrix(sapply(beta.store, function(x) {
    	            return(apply(x, 2, mean)) }), nrow=np, byrow=TRUE)
    coefs.sd.temp[, , w] <- matrix(sapply(beta.store, function(x) {
    	            return(apply(x, 2, sd))}), nrow=np, byrow=TRUE)
    for (j in 1:np) {
      fitted.temp[, w] <- fitted.temp[, w] +
                          x[, ((j - 1) * nj + 1):(j * nj)] %*% coefs.mean.temp[j, , w]
    }
    fitted.temp[, w] <- fitted.temp[, w] + alpha.mean.temp[w]
  }
  alpha.mean <- mean(alpha.mean.temp)
  alpha.sd <- mean(alpha.sd.temp)
  coefs.mean <- apply(coefs.mean.temp, c(1, 2), mean)
  coefs.sd <- apply(coefs.sd.temp, c(1, 2), mean)
  fitted <- apply(fitted.temp, 1, mean)
  xIC <- NULL
  family <- "gaussian"
  y <- c(response)
  r <- cor(c(fitted), c(y))
  r2 <- r * r
  bins <- 1:nj
  return(list(fitted=fitted, observed=y, intercept.mean=alpha.mean,
              intercept.sd=alpha.sd, coefs.mean=coefs.mean,
              coefs.sd=coefs.sd, r2=r2, family=family, bins=bins, xIC=xIC))
}