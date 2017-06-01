# FREEspline
# AN R FUNCTION FOR FITTING A HIERARCHICAL FUNCTION REGRESSION MODEL USING A GIBBS SAMPLER
# 
# Jian D. L. Yen, 11 November 2014
# 
# This code is intended for the R package FREE. It fits a function regression model with a
#   hierarchical structure on some parameters to account for grouping of subjects
#   (e.g., in space or in time). Model parameters are estimated using a Gibbs sampler
#   with slice sampling from non-standard conditional densities.
# 
# A general model description is in
# Yen J.D.L. et al. (2015) Function regression in ecology and evolution: FREE.
#   Methods in Ecology and Evolution. 6:17
#
# Arguments:
#
# y                        the values of the response function in a n x max(n_j^(i)) matrix
# x                        the values of the predictor variables in a n x n_k matrix
# groups                   the group IDs for each grouping variable q in a n x n_q matrix
# w                        the observation points for each y_ij in a n x max(n_j^(i)) matrix
# degree                   the degree of the spline functions
# n_knots_beta             number of knots in the B-splines for the fixed effects
# n_knots_gamma            number of knots in the B-splines for the quasi-random effects
# n_iter                   number of MCMC iterations
# n_burnin                 number of MCMC iterations to discard as a burn-in period
# n_thin                   thinning rate: every n_thin iterations after n_burnin will be stored
# n_chains                 number of MCMC chains
# hypers                   list with hyperparameter values for the model variance
#                            (phi_main and psi_main), the variance of the quasi-random
#                            effects (phi_gamma and psi_gamma) and the
#                            fixed effects (sigma2_beta)
# inits                    list with initial values for beta, gamma, rho, theta1, theta2,
#                            sigma2 and sigma2_gamma. These values are NOT checked prior to
#                            model fitting.

FREEscalar <- function(y,
                       x,
                       z,
                       groups,
                       bins,
                       degree = 3,
                       n_knots = 8,
                       n.iters = 1000,
                       n.burnin = round(n.iters / 5),
                       n.thin = 1,
                       n.chains = 3,
                       hypers = list(phi1 = 0.1, psi1 = 0.1,
                                     phi2 = 0.1, psi2 = 0.1,
                                     s2_alpha = 10,
                                     s2_beta = 10,
                                     s2_delta = 10),
                       inits = list(alpha = NULL, beta = NULL,
                                    gamma = NULL, delta = NULL,
                                    sigma2 = 1, sigma2_gamma = 1),
                       par.run = FALSE)
{
  # data prep
  n <- length(y)
  mod.type <- 1
  if (is.null(z)) {
    if (is.null(x)) {
      stop("at least one of x or z must be provided to FREEscalar",
           call. = FALSE)
    }
    mod.type <- 3
    z <- matrix(rnorm(2 * n), ncol = 2)
  }
  if (is.null(x)) {
    x <- vector('list', length = 1)
    x[[1]] <- matrix(rnorm((5 * n)), ncol = 5)
    bins <- vector('list', length = length(x))
    for (i in seq(along = x)) {
      bins[[i]] <- 1:ncol(x[[i]])
    }
    mod.type <- 5
  }
  if (is.null(groups)) {
    mod.type <- mod.type + 1
    groups <- matrix(sample(1:4, size = (2 * n), replace = TRUE),
                     ncol = 2)
  }
  groups <- groups - 1
  if (is.null(ncol(z))) {
    z <- matrix(z, ncol = 1)
  }
  n_k <- ncol(z)
  n_q <- ncol(groups)
  n_G_q <- apply(groups, 2, function(x) length(unique(x)))
    
  # set up grid for output of coefficients
  grid <- sort(unique(unlist(bins)))
  endpoints <- c(grid[1] - 1, grid[length(grid)] + 1)
  
  # set up priors and hypers
  phi1 <- hypers$phi1
  psi1 <- hypers$psi1
  phi2 <- hypers$phi2
  psi2 <- hypers$psi2
  s2_alpha <- hypers$s2_alpha
  s2_beta <- hypers$s2_beta
  s2_delta <- hypers$s2_delta

  # run through each chain
  spline_mod <- switch(mod.type,
                       splineScalarInternal,
                       splineScalarInternal2,
                       splineScalarInternal3,
                       splineScalarInternal4,
                       splineScalarInternal5,
                       splineScalarInternal6)
  if ((Sys.info()["sysname"] == "Darwin") & (par.run)) {
    mod <- mclapply(1:n.chains,
                    spline_mod,
                    y, x, z,
                    groups,
                    degree, n_knots,
                    n.iters, n.burnin, n.thin,
                    inits,
                    n, n_q, n_G_q, n_k,
                    grid, endpoints,
                    phi1, psi1, phi2, psi2,
                    s2_alpha, s2_beta, s2_delta,
                    mc.cores = ifelse(n.chains < 4, n.chains, 4))
  } else {
    mod <- lapply(1:n.chains,
                  spline_mod,
                  y, x, z,
                  groups,
                  degree, n_knots,
                  n.iters, n.burnin, n.thin,
                  inits,
                  n, n_q, n_G_q, n_k,
                  grid, endpoints,
                  phi1, psi1, phi2, psi2,
                  s2_alpha, s2_beta, s2_delta)
  }
  
  # summarise all outputs from all chains
  n.keep <- ncol(mod[[1]][[1]])
  fitted.mean <- apply(array(unlist(lapply(mod, function(x) x[[1]])),
                             dim = c(nrow(mod[[1]][[1]]), n.keep, n.chains)),
                       1, mean)
  fitted.sd <- apply(array(unlist(lapply(mod, function(x) x[[1]])),
                           dim = c(nrow(mod[[1]][[1]]), n.keep, n.chains)),
                     1, sd)
  coefs.mean <- apply(array(unlist(lapply(mod, function(x) x[[2]])),
                            dim = c(nrow(mod[[1]][[2]]), ncol(mod[[1]][[2]]),
                                    n.keep, n.chains)),
                      c(1, 2), mean)
  coefs.sd <- apply(array(unlist(lapply(mod, function(x) x[[2]])),
                          dim = c(nrow(mod[[1]][[2]]), ncol(mod[[1]][[2]]),
                                  n.keep, n.chains)),
                    c(1, 2), sd)
  gamma.mean <- vector("list", length = n_q)
  gamma.sd <- vector("list", length = n_q)
  for (q in 1:n_q) {
    gamma.mean[[q]] <- apply(array(unlist(lapply(mod, function(x) x[[4]][[q]])),
                                   dim = c(nrow(mod[[1]][[4]][[q]]),
                                           n.keep, n.chains)),
                             1, mean)
    gamma.sd[[q]] <- apply(array(unlist(lapply(mod, function(x) x[[4]][[q]])),
                                 dim = c(nrow(mod[[1]][[4]][[q]]),
                                         n.keep, n.chains)),
                           1, sd)
  }
  sigma2.mean <- mean(unlist(lapply(mod, function(x) x$sigma2)))
  sigma2.sd <- sd(unlist(lapply(mod, function(x) x$sigma2)))
  alpha.mean <- mean(unlist(lapply(mod, function(x) x$alpha)))
  alpha.sd <- sd(unlist(lapply(mod, function(x) x$alpha)))
  beta.mean <- apply(array(unlist(lapply(mod, function(x) x$beta)),
                           dim = c(nrow(mod[[1]]$beta), ncol(mod[[1]]$beta),
                                   n.keep, n.chains)),
                     c(1, 2), mean)
  beta.sd <- apply(array(unlist(lapply(mod, function(x) x$beta)),
                         dim = c(nrow(mod[[1]]$beta), ncol(mod[[1]]$beta),
                                 n.keep, n.chains)),
                   c(1, 2), sd)
  delta.mean <- apply(array(unlist(lapply(mod, function(x) x$delta)),
                            dim = c(nrow(mod[[1]]$delta),
                                    n.keep, n.chains)),
                      1, mean)
  delta.sd <- apply(array(unlist(lapply(mod, function(x) x$delta)),
                          dim = c(nrow(mod[[1]]$delta), n.keep, n.chains)),
                    1, sd)
  fp.sd.mean <- apply(array(unlist(lapply(mod, function(x) x$fp.sd.store)),
                            dim = c(nrow(mod[[1]]$fp.sd.store),
                                    n.keep, n.chains)),
                      1, mean)
  fp.sd.sd <- apply(array(unlist(lapply(mod, function(x) x$fp.sd.store)),
                          dim = c(nrow(mod[[1]]$fp.sd.store),
                                  n.keep, n.chains)),
                    1, sd)
  sigma2_gamma.mean <- apply(array(unlist(lapply(mod, function(x) x$sigma2_gamma)),
                                   dim = c(nrow(mod[[1]]$sigma2_gamma),
                                           n.keep, n.chains)),
                             1, mean)
  sigma2_gamma.sd <- apply(array(unlist(lapply(mod, function(x) x$sigma2_gamma)),
                                 dim = c(nrow(mod[[1]]$sigma2_gamma),
                                         n.keep, n.chains)),
                           1, sd)
  llik_all <- lapply(mod, function(x) x$llik_all)
  gamma_all <- lapply(mod, function(x) x$gamma_all)
  beta_all <- lapply(mod, function(x) x$beta_all)
  delta_all <- lapply(mod, function(x) x$delta_all)
  alpha_all <- lapply(mod, function(x) x$alpha_all)
  
  # summarise all outputs from each chain
  if (n.chains > 1) {
    coefs.chain.var <- apply(array(unlist(lapply(mod,
                                                 function(x) apply(x$coefs,
                                                                   c(1, 2),
                                                                   var))),
                                   dim = c(nrow(mod[[1]]$coefs),
                                           ncol(mod[[1]]$coefs),
                                           n.chains)),
                             c(1, 2), mean)
    gamma.chain.var <- vector("list", length = n_q)
    coefs.chain.var2 <- n.iters * apply(array(unlist(lapply(mod,
                                                            function(x) apply(x$coefs,
                                                                              c(1, 2),
                                                                              mean))),
                                              dim = c(nrow(mod[[1]]$coefs),
                                                      ncol(mod[[1]]$coefs),
                                                      n.chains)),
                                        c(1, 2), var)
    gamma.chain.var2 <- vector("list", length = n_q)
    for (q in 1:n_q) {
      gamma.chain.var[[q]] <- apply(array(unlist(lapply(mod, function(x) apply(x$gamma[[q]], 1, var))),
                                          dim=c(nrow(mod[[1]]$gamma[[q]]), n.chains)),
                                    1, mean)
      gamma.chain.var2[[q]] <- n.iters * apply(array(unlist(lapply(mod, function(x) apply(x$gamma[[q]], 1, mean))),
                                                    dim=c(nrow(mod[[1]]$gamma[[q]]), n.chains)),
                                              1, var)
    }
    sigma2.chain.var <- mean(array(unlist(lapply(mod, function(x) var(x$sigma2))),
                                   dim=c(n.chains)))
    beta.chain.var <- apply(array(unlist(lapply(mod, function(x) apply(x$beta, c(1, 2), var))),
                                  dim=c(nrow(mod[[1]]$beta), ncol(mod[[1]]$beta), n.chains)),
                            c(1, 2), mean)
    alpha.chain.var <- mean(array(unlist(lapply(mod, function(x) var(x$alpha))),
                                  dim=c(n.chains)))
    delta.chain.var <- apply(array(unlist(lapply(mod, function(x) apply(x$delta, 1, var))),
                                   dim=c(nrow(mod[[1]]$delta), n.chains)),
                             1, mean)
    sigma2_gamma.chain.var <- apply(array(unlist(lapply(mod, function(x) apply(x$sigma2_gamma, 1, var))),
                                          dim=c(nrow(mod[[1]]$sigma2_gamma), n.chains)),
                                    1, mean)
    sigma2.chain.var2 <- n.iters * var(array(unlist(lapply(mod, function(x) mean(x$sigma2))),
                                            dim=c(n.chains)))
    beta.chain.var2 <- n.iters * apply(array(unlist(lapply(mod, function(x) apply(x$beta, c(1, 2), mean))),
                                            dim=c(nrow(mod[[1]]$beta), ncol(mod[[1]]$beta), n.chains)),
                                      c(1, 2), var)
    alpha.chain.var2 <- n.iters * var(array(unlist(lapply(mod, function(x) mean(x$alpha))),
                                           dim=c(n.chains)))
    delta.chain.var2 <- n.iters * apply(array(unlist(lapply(mod, function(x) apply(x$delta, 1, mean))),
                                             dim=c(nrow(mod[[1]]$delta), n.chains)),
                                       1, var)
    sigma2_gamma.chain.var2 <- n.iters * apply(array(unlist(lapply(mod, function(x) apply(x$sigma2_gamma, 1, mean))),
                                                    dim=c(nrow(mod[[1]]$sigma2_gamma), n.chains)),
                                              1, var)
    
    # calculate Rhats
    rhat.calc <- function(chain.var, chain.var2, n.its) {
      rhat <- sqrt(((1 - (1 / n.its)) * chain.var + (1 / n.its) * chain.var2) / chain.var)
      return(rhat)
    }
    coefs.rhat <- rhat.calc(coefs.chain.var, coefs.chain.var2, n.iters)
    gamma.rhat <- vector("list", length=n_q)
    for (q in 1:n_q) {
      gamma.rhat[[q]] <- rhat.calc(gamma.chain.var[[q]], gamma.chain.var2[[q]], n.iters)
    }
    sigma2.rhat <- rhat.calc(sigma2.chain.var, sigma2.chain.var2, n.iters)
    beta.rhat <- rhat.calc(beta.chain.var, beta.chain.var2, n.iters)
    alpha.rhat <- rhat.calc(alpha.chain.var, alpha.chain.var2, n.iters)
    delta.rhat <- rhat.calc(delta.chain.var, delta.chain.var2, n.iters)
    sigma2_gamma.rhat <- rhat.calc(sigma2_gamma.chain.var, sigma2_gamma.chain.var2, n.iters)
  } else {
    coefs.rhat <- NULL
    rand.coefs.rhat <- NULL
    theta2.rhat <- NULL
    gamma.rhat <- NULL
    sigma2.rhat <- NULL
    theta1.rhat <- NULL
    beta.rhat <- NULL
    alpha.rhat <- NULL
    delta.rhat <- NULL
    rho.rhat <- NULL
    sigma2_gamma.rhat <- NULL
  }
  rhats <- list(coefs=coefs.rhat, alpha=alpha.rhat, beta=beta.rhat, gamma=gamma.rhat,
                delta=delta.rhat, sigma2=sigma2.rhat,
                sigma2_gamma=sigma2_gamma.rhat)
  
  # calculate DIC from mean log likelihood and log likelihood at mean(parameters)
  loglik.mean <- mean(unlist(lapply(mod, function(x) x$loglik)))
  theta <- seq(5, max(sapply(x, ncol)) - 4, length=n_knots - degree)
  bs_beta <- vector('list', length = length(x))
  for (i in seq(along = bs_beta)) {
    bs_beta[[i]] <- calc_bs(1:ncol(x[[i]]), theta, degree, c(0, ncol(x[[i]]) + 1))
  }
  loglik.param.bar <- lnL_scalar(y=y, x=x, groups=groups, beta=beta.mean, gamma=gamma.mean,
                          delta=delta.mean, z=z, alpha=alpha.mean, sigma2=sigma2.mean,
                          bs_beta=bs_beta)
  DIC <- -2 * loglik.mean + 2 * loglik.param.bar
      
  # calculate all returned values
  family <- "gaussian"
  if (mod.type == 2) {
    # z but no groups
    gamma.mean <- NULL
    gamma.sd <- NULL
    sigma2_gamma.mean <- NULL
    sigma2_gamma.sd <- NULL
    fp.sd.mean <- NULL
    fp.sd.sd <- NULL
  }
  if (mod.type == 3) {
    # groups but no z
    delta.mean <- NULL
    delta.sd <- NULL
  }
  if (mod.type == 4) {
    # no z and no groups
    gamma.mean <- NULL
    gamma.sd <- NULL
    sigma2_gamma.mean <- NULL
    sigma2_gamma.sd <- NULL
    fp.sd.mean <- NULL
    fp.sd.sd <- NULL
    delta.mean <- NULL
    delta.sd <- NULL
  }
  if (mod.type == 5) {
    # groups and z but no x
    beta.mean <- NULL
    beta.sd <- NULL
    coefs.mean <- NULL
    coefs.sd <- NULL
  }
  if (mod.type == 6) {
    # no x and no groups
    gamma.mean <- NULL
    gamma.sd <- NULL
    sigma2_gamma.mean <- NULL
    sigma2_gamma.sd <- NULL
    fp.sd.mean <- NULL
    fp.sd.sd <- NULL
    beta.mean <- NULL
    beta.sd <- NULL
    coefs.mean <- NULL
    coefs.sd <- NULL
  }
  r2 <- cor(c(fitted.mean), c(y)) * cor(c(fitted.mean), c(y), use = "complete")
  return(list(fitted = fitted.mean, fitted.sd = fitted.sd,
              observed = y,
              coefs.mean = coefs.mean,
              coefs.sd = coefs.sd,
              r2 = r2, family = family,
              DIC = DIC, rhats = rhats,
              sigma2.mean = sigma2.mean,
              sigma2.sd = sigma2.sd,
              sigma2_gamma.mean = sigma2_gamma.mean, 
              sigma2_gamma.sd = sigma2_gamma.sd,
              alpha.mean = alpha.mean,
              alpha.sd = alpha.sd,
              beta.mean = beta.mean,
              beta.sd = beta.sd,
              gamma.mean = gamma.mean,
              gamma.sd = gamma.sd,
              delta.mean = delta.mean,
              delta.sd = delta.sd,
              llik_all = llik_all, fp.sd.mean = fp.sd.mean,
              fp.sd.sd = fp.sd.sd))
}