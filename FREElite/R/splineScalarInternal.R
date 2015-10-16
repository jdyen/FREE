splineScalarInternal <- function(chain, y, x, z, groups, degree, n_knots, n.iter,
                                   n.burnin, n.thin, inits, n, n_q, n_G_q, n_k,
                                   grid, endpoints, phi1, psi1, phi2, psi2, s2_alpha,
                                   s2_beta, s2_delta)
{
  # prepare output objects
  n_keep <- (n.iter - n.burnin) / n.thin
  fitted.store <- array(dim=c(n, n_keep))
  loglik.store <- vector("numeric", length=n_keep)
  sigma2.store <- vector("numeric", length=n_keep)
  sigma2_gamma.store <-  array(dim=c(n_q, n_keep))
  coefs.store <-  array(dim=c(length(grid), n_keep))
  beta.store <- array(dim=c(n_knots, n_keep))
  alpha.store <- vector("numeric", length=n_keep)
  delta.store <- array(dim=c(n_k, n_keep))
  fp.sd.store <- array(dim=c(n_q, n_keep))
  gamma.store <- vector("list", length=n_q)
  for (q in 1:n_q) {
    gamma.store[[q]] <- array(dim=c(n_G_q[q], n_keep))
  }  
  llik_all <- vector("numeric", length=n.iter)
  # initialise chain
  if (!is.null(inits$alpha)) {
    alpha <- inits$alpha
  } else {
    alpha <- rnorm(1)
  }
  if (!is.null(inits$beta)) {
    beta <- inits$beta
  } else {
    beta <- rnorm(n_knots)
  }
  if (!is.null(inits$gamma)) {
      gamma <- inits$gamma
  } else {
    gamma <- vector("list", length=n_q)
    for (q in 1:n_q) {
      gamma[[q]] <- rnorm(n_G_q[q] - 1)
      gamma[[q]] <- c(gamma[[q]], -sum(gamma[[q]][1:(n_G_q[q] - 1)]))
    }
  }
  if (!is.null(inits$delta)) {
    delta <- inits$delta
  } else {
    delta <- rnorm(n_k)
  }
  theta <- seq(5, ncol(x) - 4, length=n_knots - degree)
  bs_beta <- calc_bs_scalar(1:ncol(x), theta, degree, c(0, ncol(x) + 1))
  sigma2 <- inits$sigma2
  if (length(inits$sigma2_gamma) == n_q) {
    sigma2_gamma <- inits$sigma2_gamma
  } else {
    sigma2_gamma <- rep(inits$sigma2_gamma, n_q)
  }
    
  llik_init <- lnL_scalar(y=y, x=x, groups=as.matrix(groups), beta=beta, gamma=gamma, delta=delta,
                   z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta)
  # run through each iteration
  for (i in 1:n.iter) {
    tmp <- updateBetasScalar(y=y, x=x, z=z, groups=groups, alpha=alpha,
                               beta=beta, gamma=gamma, delta=delta,
                               sigma2=sigma2, sigma2_gamma=sigma2_gamma,
                               bs_beta=bs_beta, phi1=phi1, psi1=psi1,
                               phi2=phi2, psi2=psi2,
                               s2_alpha, s2_beta=s2_beta, s2_delta,
                               n_q, n_G_q)
    alpha <- tmp$alpha
    beta <- tmp$beta
    gamma <- tmp$gamma
    delta <- tmp$delta
    sigma2 <- tmp$sigma2
    sigma2_gamma <- tmp$sigma2_gamma

    llik_all[i] <- tmp$lnL

    # save outputs
    if ({i > n.burnin} & {{i %% n.thin} == 0}) {
      current_iter <- {{i - n.burnin} / n.thin}
      fitted.store[, current_iter] <- fitted_scalar(x=x, z=z, groups=groups,
                                                    beta=beta, gamma=gamma,
                                                    delta=delta, bs_beta=bs_beta)
      loglik.store[current_iter] <- llik_all[i]
      coefs.store[, current_iter] <- coefs_calc_scalar(beta=beta, theta=theta, degree=degree,
                                                grid=grid, endpoints=endpoints)
      sigma2.store[current_iter] <- sigma2
      sigma2_gamma.store[, current_iter] <- sigma2_gamma
      alpha.store[current_iter] <- alpha
      beta.store[, current_iter] <- beta
      delta.store[, current_iter] <- delta
      for (q in 1:n_q) {
        gamma.store[[q]][, current_iter] <- gamma[[q]]
        fp.sd.store[q, current_iter] <- sd(gamma[[q]])
      }
    }
  }
  allout <- list(fitted=fitted.store, coefs=coefs.store, loglik=loglik.store,
                 gamma=gamma.store, sigma2=sigma2.store, beta=beta.store,
                 sigma2_gamma=sigma2_gamma.store, llik_all=llik_all,
                 bs_beta=bs_beta, alpha=alpha.store, delta=delta.store,
                 fp.sd.store=fp.sd.store)
  return(allout)
}

splineScalarInternal2 <- function(chain, y, x, z, groups, degree, n_knots, n.iter,
                                    n.burnin, n.thin, inits, n, n_q, n_G_q, n_k,
                                    grid, endpoints, phi1, psi1, phi2, psi2, s2_alpha,
                                    s2_beta, s2_delta)
{
  # prepare output objects
  n_keep <- (n.iter - n.burnin) / n.thin
  fitted.store <- array(dim=c(n, n_keep))
  loglik.store <- vector("numeric", length=n_keep)
  sigma2.store <- vector("numeric", length=n_keep)
  sigma2_gamma.store <-  array(dim=c(n_q, n_keep))
  coefs.store <-  array(dim=c(length(grid), n_keep))
  beta.store <- array(dim=c(n_knots, n_keep))
  alpha.store <- vector("numeric", length=n_keep)
  delta.store <- array(dim=c(n_k, n_keep))
  fp.sd.store <- array(dim=c(n_q, n_keep))
  gamma.store <- vector("list", length=n_q)
  for (q in 1:n_q) {
    gamma.store[[q]] <- array(dim=c(n_G_q[q], n_keep))
  }  
  llik_all <- vector("numeric", length=n.iter)
  # initialise chain
  if (!is.null(inits$alpha)) {
    alpha <- inits$alpha
  } else {
    alpha <- rnorm(1)
  }
  if (!is.null(inits$beta)) {
    beta <- inits$beta
  } else {
    beta <- rnorm(n_knots)
  }
  gamma <- vector("list", length=n_q)
  for (q in 1:n_q) {
    gamma[[q]] <- rep(0, n_G_q[q])
  }
  if (!is.null(inits$delta)) {
    delta <- inits$delta
  } else {
    delta <- rnorm(n_k)
  }
  theta <- seq(5, ncol(x) - 4, length=n_knots - degree)
  bs_beta <- calc_bs_scalar(1:ncol(x), theta, degree, c(0, ncol(x) + 1))
  sigma2 <- inits$sigma2
  if (length(inits$sigma2_gamma) == n_q) {
    sigma2_gamma <- inits$sigma2_gamma
  } else {
    sigma2_gamma <- rep(inits$sigma2_gamma, n_q)
  }
  
  llik_init <- lnL_scalar(y=y, x=x, groups=as.matrix(groups), beta=beta, gamma=gamma, delta=delta,
                   z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta)
  # run through each iteration
  for (i in 1:n.iter) {
    tmp <- updateBetasScalar2(y=y, x=x, z=z, groups=groups, alpha=alpha,
                                beta=beta, gamma=gamma, delta=delta,
                                sigma2=sigma2, sigma2_gamma=sigma2_gamma,
                                bs_beta=bs_beta, phi1=phi1, psi1=psi1,
                                phi2=phi2, psi2=psi2,
                                s2_alpha, s2_beta=s2_beta, s2_delta,
                                n_q, n_G_q)
    alpha <- tmp$alpha
    beta <- tmp$beta
    gamma <- tmp$gamma
    delta <- tmp$delta
    sigma2 <- tmp$sigma2
    sigma2_gamma <- tmp$sigma2_gamma
    
    llik_all[i] <- tmp$lnL
    
    # save outputs
    if ({i > n.burnin} & {{i %% n.thin} == 0}) {
      current_iter <- {{i - n.burnin} / n.thin}
      fitted.store[, current_iter] <- fitted_scalar(x=x, z=z, groups=groups,
                                                    beta=beta, gamma=gamma,
                                                    delta=delta, bs_beta=bs_beta)
      loglik.store[current_iter] <- llik_all[i]
      coefs.store[, current_iter] <- coefs_calc_scalar(beta=beta, theta=theta, degree=degree,
                                                grid=grid, endpoints=endpoints)
      sigma2.store[current_iter] <- sigma2
      sigma2_gamma.store[, current_iter] <- sigma2_gamma
      alpha.store[current_iter] <- alpha
      beta.store[, current_iter] <- beta
      delta.store[, current_iter] <- delta
      for (q in 1:n_q) {
        gamma.store[[q]][, current_iter] <- gamma[[q]]
        fp.sd.store[q, current_iter] <- sd(gamma[[q]])
      }
    }
  }
  allout <- list(fitted=fitted.store, coefs=coefs.store, loglik=loglik.store,
                 gamma=gamma.store, sigma2=sigma2.store, beta=beta.store,
                 sigma2_gamma=sigma2_gamma.store, llik_all=llik_all,
                 bs_beta=bs_beta, alpha=alpha.store, delta=delta.store,
                 fp.sd.store=fp.sd.store)
  return(allout)
}

splineScalarInternal3 <- function(chain, y, x, z, groups, degree, n_knots, n.iter,
                                    n.burnin, n.thin, inits, n, n_q, n_G_q, n_k,
                                    grid, endpoints, phi1, psi1, phi2, psi2, s2_alpha,
                                    s2_beta, s2_delta)
{
  # prepare output objects
  n_keep <- (n.iter - n.burnin) / n.thin
  fitted.store <- array(dim=c(n, n_keep))
  loglik.store <- vector("numeric", length=n_keep)
  sigma2.store <- vector("numeric", length=n_keep)
  sigma2_gamma.store <-  array(dim=c(n_q, n_keep))
  coefs.store <-  array(dim=c(length(grid), n_keep))
  beta.store <- array(dim=c(n_knots, n_keep))
  alpha.store <- vector("numeric", length=n_keep)
  delta.store <- array(dim=c(n_k, n_keep))
  fp.sd.store <- array(dim=c(n_q, n_keep))
  gamma.store <- vector("list", length=n_q)
  for (q in 1:n_q) {
    gamma.store[[q]] <- array(dim=c(n_G_q[q], n_keep))
  }  
  llik_all <- vector("numeric", length=n.iter)
  # initialise chain
  if (!is.null(inits$alpha)) {
    alpha <- inits$alpha
  } else {
    alpha <- rnorm(1)
  }
  if (!is.null(inits$beta)) {
    beta <- inits$beta
  } else {
    beta <- rnorm(n_knots)
  }
  if (!is.null(inits$gamma)) {
    gamma <- inits$gamma
  } else {
    gamma <- vector("list", length=n_q)
    for (q in 1:n_q) {
      gamma[[q]] <- rnorm(n_G_q[q] - 1)
      gamma[[q]] <- c(gamma[[q]], -sum(gamma[[q]][1:(n_G_q[q] - 1)]))
      fp.sd.store[q, current_iter] <- sd(gamma[[q]])
    }
  }
  delta <- rep(0, n_k)
  theta <- seq(5, ncol(x) - 4, length=n_knots - degree)
  bs_beta <- calc_bs_scalar(1:ncol(x), theta, degree, c(0, ncol(x) + 1))
  sigma2 <- inits$sigma2
  if (length(inits$sigma2_gamma) == n_q) {
    sigma2_gamma <- inits$sigma2_gamma
  } else {
    sigma2_gamma <- rep(inits$sigma2_gamma, n_q)
  }
  
  llik_init <- lnL_scalar(y=y, x=x, groups=as.matrix(groups), beta=beta, gamma=gamma, delta=delta,
                   z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta)
  # run through each iteration
  for (i in 1:n.iter) {
    tmp <- updateBetasScalar3(y=y, x=x, z=z, groups=groups, alpha=alpha,
                                beta=beta, gamma=gamma, delta=delta,
                                sigma2=sigma2, sigma2_gamma=sigma2_gamma,
                                bs_beta=bs_beta, phi1=phi1, psi1=psi1,
                                phi2=phi2, psi2=psi2,
                                s2_alpha, s2_beta=s2_beta, s2_delta,
                                n_q, n_G_q)
    alpha <- tmp$alpha
    beta <- tmp$beta
    gamma <- tmp$gamma
    delta <- tmp$delta
    sigma2 <- tmp$sigma2
    sigma2_gamma <- tmp$sigma2_gamma
    
    llik_all[i] <- tmp$lnL
    
    # save outputs
    if ({i > n.burnin} & {{i %% n.thin} == 0}) {
      current_iter <- {{i - n.burnin} / n.thin}
      fitted.store[, current_iter] <- fitted_scalar(x=x, z=z, groups=groups,
                                                    beta=beta, gamma=gamma,
                                                    delta=delta, bs_beta=bs_beta)
      loglik.store[current_iter] <- llik_all[i]
      coefs.store[, current_iter] <- coefs_calc_scalar(beta=beta, theta=theta, degree=degree,
                                                grid=grid, endpoints=endpoints)
      sigma2.store[current_iter] <- sigma2
      sigma2_gamma.store[, current_iter] <- sigma2_gamma
      alpha.store[current_iter] <- alpha
      beta.store[, current_iter] <- beta
      delta.store[, current_iter] <- delta
      for (q in 1:n_q) {
        gamma.store[[q]][, current_iter] <- gamma[[q]]
        fp.sd.store[q, current_iter] <- sd(gamma[[q]])
      }
    }
  }
  allout <- list(fitted=fitted.store, coefs=coefs.store, loglik=loglik.store,
                 gamma=gamma.store, sigma2=sigma2.store, beta=beta.store,
                 sigma2_gamma=sigma2_gamma.store, llik_all=llik_all,
                 bs_beta=bs_beta, alpha=alpha.store, delta=delta.store,
                 fp.sd.store=fp.sd.store)
  return(allout)
}

splineScalarInternal4 <- function(chain, y, x, z, groups, degree, n_knots, n.iter,
                                    n.burnin, n.thin, inits, n, n_q, n_G_q, n_k,
                                    grid, endpoints, phi1, psi1, phi2, psi2, s2_alpha,
                                    s2_beta, s2_delta)
{
  # prepare output objects
  n_keep <- (n.iter - n.burnin) / n.thin
  fitted.store <- array(dim=c(n, n_keep))
  loglik.store <- vector("numeric", length=n_keep)
  sigma2.store <- vector("numeric", length=n_keep)
  sigma2_gamma.store <-  array(dim=c(n_q, n_keep))
  coefs.store <-  array(dim=c(length(grid), n_keep))
  beta.store <- array(dim=c(n_knots, n_keep))
  alpha.store <- vector("numeric", length=n_keep)
  delta.store <- array(dim=c(n_k, n_keep))
  fp.sd.store <- array(dim=c(n_q, n_keep))
  gamma.store <- vector("list", length=n_q)
  for (q in 1:n_q) {
    gamma.store[[q]] <- array(dim=c(n_G_q[q], n_keep))
  }  
  llik_all <- vector("numeric", length=n.iter)
  # initialise chain
  if (!is.null(inits$alpha)) {
    alpha <- inits$alpha
  } else {
    alpha <- rnorm(1)
  }
  if (!is.null(inits$beta)) {
    beta <- inits$beta
  } else {
    beta <- rnorm(n_knots)
  }
  gamma <- vector("list", length=n_q)
  for (q in 1:n_q) {
    gamma[[q]] <- rep(0, n_G_q[q])
  }
  delta <- rep(0, n_k)
  theta <- seq(5, ncol(x) - 4, length=n_knots - degree)
  bs_beta <- calc_bs_scalar(1:ncol(x), theta, degree, c(0, ncol(x) + 1))
  sigma2 <- inits$sigma2
  if (length(inits$sigma2_gamma) == n_q) {
    sigma2_gamma <- inits$sigma2_gamma
  } else {
    sigma2_gamma <- rep(inits$sigma2_gamma, n_q)
  }
  
  llik_init <- lnL_scalar(y=y, x=x, groups=as.matrix(groups), beta=beta, gamma=gamma, delta=delta,
                   z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta)
  # run through each iteration
  for (i in 1:n.iter) {
    tmp <- updateBetasScalar4(y=y, x=x, z=z, groups=groups, alpha=alpha,
                                beta=beta, gamma=gamma, delta=delta,
                                sigma2=sigma2, sigma2_gamma=sigma2_gamma,
                                bs_beta=bs_beta, phi1=phi1, psi1=psi1,
                                phi2=phi2, psi2=psi2,
                                s2_alpha, s2_beta=s2_beta, s2_delta,
                                n_q, n_G_q)
    alpha <- tmp$alpha
    beta <- tmp$beta
    gamma <- tmp$gamma
    delta <- tmp$delta
    sigma2 <- tmp$sigma2
    sigma2_gamma <- tmp$sigma2_gamma
    
    llik_all[i] <- tmp$lnL
    
    # save outputs
    if ({i > n.burnin} & {{i %% n.thin} == 0}) {
      current_iter <- {{i - n.burnin} / n.thin}
      fitted.store[, current_iter] <- fitted_scalar(x=x, z=z, groups=groups,
                                                    beta=beta, gamma=gamma,
                                                    delta=delta, bs_beta=bs_beta)
      loglik.store[current_iter] <- llik_all[i]
      coefs.store[, current_iter] <- coefs_calc_scalar(beta=beta, theta=theta, degree=degree,
                                                grid=grid, endpoints=endpoints)
      sigma2.store[current_iter] <- sigma2
      sigma2_gamma.store[, current_iter] <- sigma2_gamma
      alpha.store[current_iter] <- alpha
      beta.store[, current_iter] <- beta
      delta.store[, current_iter] <- delta
      for (q in 1:n_q) {
        gamma.store[[q]][, current_iter] <- gamma[[q]]
        fp.sd.store[q, current_iter] <- sd(gamma[[q]])
      }
    }
  }
  allout <- list(fitted=fitted.store, coefs=coefs.store, loglik=loglik.store,
                 gamma=gamma.store, sigma2=sigma2.store, beta=beta.store,
                 sigma2_gamma=sigma2_gamma.store, llik_all=llik_all,
                 bs_beta=bs_beta, alpha=alpha.store, delta=delta.store,
                 fp.sd.store=fp.sd.store)
  return(allout)
}