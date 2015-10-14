spline_slice_internal <- function(chain, y, x, groups, w, degree, n_knots_beta, n_knots_gamma, n.iter,
                                  n.burnin, n.thin, inits, n, n_j, n_k, n_p, n_q, n_G_q,
                                  n_t, grid, endpoints, sigma2_hyper_a, sigma2_hyper_b,
                                  sigma2_gamma_hyper_a, sigma2_gamma_hyper_b, beta_hyper)
{
  # prepare output objects
  n_keep <- (n.iter - n.burnin) / n.thin
  fitted.store <- array(dim=c(n, max(n_j, na.rm=TRUE), n_keep))
  loglik.store <- vector("numeric", length=n_keep)
  rho.store <- vector("numeric", length=n_keep)
  sigma2.store <- vector("numeric", length=n_keep)
  sigma2_gamma.store <-  array(dim=c(n_q, n_keep))
  coefs.store <-  array(dim=c(n_k, length(grid), n_keep))
  beta.store <- vector("list", length=n_k)
  theta1.store <- vector("list", length=n_k)
  for (k in 1:n_k) {
    beta.store[[k]] <- array(dim=c(n_p[k], n_keep))
    theta1.store[[k]] <- array(dim=c(n_knots_beta[k], n_keep))
  }
  rand.coefs.store <- vector("list", length=n_q)
  gamma.store <- vector("list", length=n_q)
  theta2.store <- vector("list", length=n_q)
  for (q in 1:n_q) {
    rand.coefs.store[[q]] <- array(dim=c(n_G_q[q], length(grid), n_keep))
    gamma.store[[q]] <- vector("list", length=n_G_q[q])
    theta2.store[[q]] <- vector("list", length=n_G_q[q])
    for (qq in 1:n_G_q[q]) {
      gamma.store[[q]][[qq]] <- array(dim=c(n_t[[q]][qq], n_keep))
      theta2.store[[q]][[qq]] <- array(dim=c(n_knots_gamma[[q]][qq], n_keep))
    }
  }  
  llik_all <- vector("numeric", length=n.iter)
  beta_all <- vector("list", length=n_k)
  gamma_all <- vector("list", length=n_q)
  for (k in 1:n_k) {
    beta_all[[k]] <- array(dim=c(n_p[k], n.iter))
  }
  for (q in 1:n_q) {
    gamma_all[[q]] <- vector("list", length=n_G_q[q])
    for (qq in 1:n_G_q[q]) {
      gamma_all[[q]][[qq]] <- array(dim=c(n_t[[q]][qq], n.iter))
    }
  }
  # initialise chain
  if (!is.null(inits$beta)) {
    beta <- inits$beta
  } else {
    beta <- vector("list", length=n_k)
    for (k in 1:n_k) {
      beta[[k]] <- rnorm(n_p[k])
    }
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
      gamma[[q]] <- vector("list", length=n_G_q[q])
      for (j in 1:n_G_q[q]) {
        gamma[[q]][[j]] <- rnorm(n_t[[q]][j])
      }
    }
  }
  if (is.null(inits$theta1)) {
    theta1 <- vector("list", length=n_k)
    for (k in 1:n_k) {
      theta1[[k]] <- seq(min(w, na.rm=TRUE) - 0.5, max(w, na.rm=TRUE) + 0.5,
                         length=n_knots_beta[k])
    }
  }
  if (is.null(inits$theta2)) {
    theta2 <- vector("list", length=n_q)
    for (q in 1:n_q) {
      theta2[[q]] <- vector("list", length=n_G_q[q])
      for (qq in 1:n_G_q[q]) {
        theta2[[q]][[qq]] <- seq(min(w, na.rm=TRUE) - 0.5, max(w, na.rm=TRUE) + 0.5,
                                 length=n_knots_gamma[[q]][qq])
      }
    }
  }
  b_splines_beta <- vector("list", length=n_k)
  for (i in 1:n_k) {
    b_splines_beta[[i]] <- calc_bs(w, theta1[[i]], degree,
                                   c(min(w, na.rm=TRUE) - 1, max(w, na.rm=TRUE) + 1))
  }
  b_splines_gamma <- vector("list", length=n_q)
  for (i in 1:n_q) {
    b_splines_gamma[[i]] <- vector("list", length=n_G_q[i])
    for (q in 1:n_G_q[i]) {
      b_splines_gamma[[i]][[q]] <- calc_bs(w, theta2[[i]][[q]], degree,
                                           c(min(w, na.rm=TRUE) - 1, max(w, na.rm=TRUE) + 1))
    }
  }
  rho <- inits$rho
  sigma2 <- inits$sigma2
  if (length(inits$sigma2_gamma) == n_q) {
    sigma2_gamma <- inits$sigma2_gamma
  } else {
    sigma2_gamma <- rep(inits$sigma2_gamma, n_q)
  }
    
  llik_init <- lnL(y, x, w, groups, beta, gamma,
                   theta1, theta2, sigma2, rho,
                   b_splines_beta, b_splines_gamma,
                   sigma2_hyper_a, sigma2_hyper_b,
                   sigma2_gamma_hyper_a, sigma2_gamma_hyper_b,
                   beta_hyper, n, n_j, n_k, n_q)
  # run through each iteration
  for (i in 1:n.iter) {
    tmp <- update_betas(y=y, x=x, groups=groups, w=w, beta=beta, gamma=gamma, rho=rho,
                        theta1=theta1, theta2=theta2, sigma2=sigma2,
                        sigma2_gamma=sigma2_gamma, b_splines_beta=b_splines_beta,
                        b_splines_gamma=b_splines_gamma,
                        sigma2_hyper_a=sigma2_hyper_a, sigma2_hyper_b=sigma2_hyper_b,
                        sigma2_gamma_hyper_a=sigma2_gamma_hyper_a,
                        sigma2_gamma_hyper_b=sigma2_gamma_hyper_b,
                        beta_hyper=beta_hyper, n=n, n_j=n_j, n_k=n_k, n_q=n_q,
                        n_G_q=n_G_q, n_p=n_p, n_t=n_t, degree=degree)
    beta <- tmp$beta
    gamma <- tmp$gamma
    rho <- tmp$rho
    theta1 <- tmp$theta1
    theta2 <- tmp$theta2
    sigma2 <- tmp$sigma2
    sigma2_gamma <- tmp$sigma2_gamma
    b_splines_beta <- tmp$b_splines_beta
    b_splines_gamma <- tmp$b_splines_gamma

    llik_all[i] <- lnL(y, x, w, groups, beta, gamma,
                       theta1, theta2, sigma2, rho,
                       b_splines_beta, b_splines_gamma,
                       sigma2_hyper_a, sigma2_hyper_b,
                       sigma2_gamma_hyper_a, sigma2_gamma_hyper_b,
                       beta_hyper, n, n_j, n_k, n_q)
    for (k in 1:n_k) {
      beta_all[[k]][, i] <- beta[[k]]
    }
    for (q in 1:n_q) {
      for (qq in 1:n_G_q[q]) {
        gamma_all[[q]][[qq]][, i] <- gamma[[q]][[qq]]
      }
    }

    # save outputs
    if ({i > n.burnin} & {{i %% n.thin} == 0}) {
      current_iter <- {{i - n.burnin} / n.thin}
      fitted.store[, , current_iter] <- fitted_calc(w=w, x=x, groups=groups,
                         beta=beta, gamma=gamma, theta1=theta1, theta2=theta2, n=n,
                         n_j=n_j, n_k=n_k, n_q=n_q, n_G_q=n_G_q, n_p=n_p, n_t=n_t,
                         degree=degree)
       loglik.store[current_iter] <- lnL(y, x, w, groups, beta=beta, gamma=gamma,
                   theta1, theta2, sigma2, rho,
                   b_splines_beta, b_splines_gamma,
                   sigma2_hyper_a, sigma2_hyper_b,
                   sigma2_gamma_hyper_a, sigma2_gamma_hyper_b,
                   beta_hyper, n, n_j, n_k, n_q)
      coefs.store[, , current_iter] <- coefs_calc(beta=beta, theta=theta1,
    	                                                 degree=degree, grid=grid,
    	                                                 endpoints=endpoints)
      rho.store[current_iter] <- rho
      sigma2.store[current_iter] <- sigma2
      sigma2_gamma.store[, current_iter] <- sigma2_gamma
      for (k in 1:n_k) {
        beta.store[[k]][, current_iter] <- beta[[k]]
        theta1.store[[k]][, current_iter] <- theta1[[k]]
      }
      for (q in 1:n_q) {
        rand.coefs.store[[q]][, , current_iter] <- coefs_calc(beta=gamma[[q]],
                                                   theta=theta2[[q]], degree=degree,
                                                   grid=grid, endpoints=endpoints)
        for (qq in 1:n_G_q[q]) {
          gamma.store[[q]][[qq]][, current_iter] <- gamma[[q]][[qq]]
          theta2.store[[q]][[qq]][, current_iter] <- theta2[[q]][[qq]]
        }
      }
    }
  }
  allout <- list(fitted.store, coefs.store, loglik.store, rand.coefs.store, gamma.store,
                 theta1.store, theta2.store, rho.store, sigma2.store, beta.store,
                 sigma2_gamma.store, llik_all, gamma_all, beta_all)
  return(allout)
}