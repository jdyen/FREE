splineSliceInternal <- function(chain, y, x, groups, w, degree, n_knots, n.iter,
n.burnin, n.thin, inits, n, n_j, n_k, n_p, n_q, n_G_q,
n_t, grid, endpoints, sigma2_hyper_a, sigma2_hyper_b,
sigma2_gamma_hyper_a, sigma2_gamma_hyper_b, beta_hyper,
bin_id)
{
  # prepare output objects
  n_keep <- (n.iter - n.burnin) / n.thin
  fitted.store <- array(dim=c(n, max(n_j, na.rm=TRUE), n_keep))
  loglik.store <- vector("numeric", length=n_keep)
  rho.store <- vector("numeric", length=n_keep)
  sigma2.store <- vector("numeric", length=n_keep)
  sigma2_gamma.store <-  array(dim=c(n_q, n_keep))
  coefs.store <-  array(dim=c(n_k, length(grid), n_keep))
  beta.store <- array(dim=c(n_k, n_p, n_keep))
  theta1.store <- array(dim=c(n_knots, n_keep))
  rand.coefs.store <- vector("list", length=n_q)
  gamma.store <- vector("list", length=n_q)
  for (q in 1:n_q) {
    rand.coefs.store[[q]] <- array(dim=c(n_G_q[q], length(grid), n_keep))
    gamma.store[[q]] <- array(dim=c(n_G_q[q], n_t, n_keep))
  }
  llik_all <- vector("numeric", length=n.iter)
  # initialise chain
  if (!is.null(inits$beta)) {
    beta <- inits$beta
  } else {
    beta <- matrix(rnorm((n_k * n_p), mean=0, sd=0.1), ncol=n_p)
  }
  if (!is.null(inits$gamma)) {
    if (is.vector(gamma, mode="list") && all(sapply(gamma, is.vector, mode="list"))) {
      gamma <- inits$gamma
    } else {
      stop("Error in FREEspline: initial values for gamma have the wrong dimensions...",
      call.=FALSE)
    }
  } else {
    gamma <- vector("list", length=n_q)
    for (q in 1:n_q) {
      gamma[[q]] <- matrix(rnorm((n_G_q[q] * n_t), mean=0, sd=0.1), ncol=n_t)
    }
  }
  theta1 <- seq(min(w, na.rm=TRUE) - 0.5, max(w, na.rm=TRUE) + 0.5,
  length=n_knots)
  b_splines_mat <- calc_bs(w, theta1, degree,
  c(min(w, na.rm=TRUE) - 1, max(w, na.rm=TRUE) + 1))
  rho <- inits$rho
  sigma2 <- inits$sigma2
  if (length(inits$sigma2_gamma) == n_q) {
    sigma2_gamma <- inits$sigma2_gamma
  } else {
    sigma2_gamma <- rep(inits$sigma2_gamma, n_q)
  }
  
  llik_init <- lnL(y, x, groups, beta, gamma,
  sigma2, rho, b_splines_mat,
  beta_hyper, n, n_j, n_k, n_q, bin_id)
  # run through each iteration
  for (i in 1:n.iter) {
    tmp <- updateBetas(y=y, x=x, groups=groups, beta=beta, gamma=gamma, rho=rho,
    sigma2=sigma2, sigma2_gamma=sigma2_gamma, b_splines_mat=b_splines_mat,
    sigma2_hyper_a=sigma2_hyper_a, sigma2_hyper_b=sigma2_hyper_b,
    sigma2_gamma_hyper_a=sigma2_gamma_hyper_a,
    sigma2_gamma_hyper_b=sigma2_gamma_hyper_b, beta_hyper=beta_hyper,
    n=n, n_j=n_j, n_k=n_k, n_q=n_q, n_G_q=n_G_q, n_p=n_p, n_t=n_t,
    bin_id=bin_id)
    beta <- tmp$beta
    gamma <- tmp$gamma
    rho <- tmp$rho
    sigma2 <- tmp$sigma2
    sigma2_gamma <- tmp$sigma2_gamma
    
    llik_all[i] <- lnL(y, x, groups, beta, gamma,
    sigma2, rho, b_splines_mat,
    beta_hyper, n, n_j, n_k, n_q, bin_id)
    
    # save outputs
    if ({i > n.burnin} & {{i %% n.thin} == 0}) {
      current_iter <- {{i - n.burnin} / n.thin}
      fitted.store[, , current_iter] <- fitted_calc(x=x, groups=groups, beta=beta,
      gamma=gamma, b_splines_mat=b_splines_mat,
      n=n, n_j=n_j, n_q=n_q, bin_id)
      loglik.store[current_iter] <- lnL(y, x, groups, beta=beta, gamma=gamma,
      sigma2, rho, b_splines_mat,
      beta_hyper, n, n_j, n_k, n_q, bin_id)
      coefs.store[, , current_iter] <- coefs_calc(beta=beta, theta=theta1,
      degree=degree, grid=grid,
      endpoints=endpoints)
      rho.store[current_iter] <- rho
      sigma2.store[current_iter] <- sigma2
      sigma2_gamma.store[, current_iter] <- sigma2_gamma
      beta.store[, , current_iter] <- beta
      for (q in 1:n_q) {
        rand.coefs.store[[q]][, , current_iter] <- coefs_calc(beta=gamma[[q]],
        theta=theta1, degree=degree,
        grid=grid, endpoints=endpoints)
        gamma.store[[q]][, , current_iter] <- gamma[[q]]
      }
    }
  }
  allout <- list(fitted.store, coefs.store, loglik.store, rand.coefs.store, gamma.store,
  rho.store, sigma2.store, beta.store, sigma2_gamma.store, llik_all)
  return(allout)
}