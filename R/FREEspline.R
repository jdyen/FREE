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
# degree                   the degree of the spline functions
# n_knots                  number of knots in the B-splines
# n_iter                   number of MCMC iterations
# n_burnin                 number of MCMC iterations to discard as a burn-in period
# n_thin                   thinning rate: every n_thin iterations after n_burnin will be stored
# n_chains                 number of MCMC chains
# hypers                   list with hyperparameter values for the model variance
#                            (phi_main and psi_main), the variance of the quasi-random
#                            effects (phi_gamma and psi_gamma) and the
#                            fixed effects (sigma2_beta)
# inits                    list with initial values for beta, gamma, rho,
#                            sigma2 and sigma2_gamma. These values are NOT checked prior to
#                            model fitting.
# par.run

FREEspline <- function(y, x, groups, bins=NULL, degree=3, n_knots=5, n.iters=1000,
                       n.burnin=round(n.iters / 5), n.thin=1, n.chains=3,
                       hypers=list(psi_main=0.1, phi_main=0.1, psi_gamma=0.1,
                                   phi_gamma=0.1, sigma2_beta=10),
                       inits=list(beta=NULL, gamma=NULL, rho=0.5,
                                  sigma2=1, sigma2_gamma=1),
                       par.run=FALSE)
{
  # data prep
  if (is.null(x)) {
    x <- matrix(rep(1, nrow(y)), ncol=1)
  } else {
    x <- as.matrix(cbind(rep(1, nrow(y)), x))
  }
  n <- nrow(y)
  n_j <- apply(y, 1, function(x) sum(!is.na(x)))
  n_k <- ncol(x)
  n_p <- n_knots + degree
  n_t <- n_p
  n_q <- ncol(groups)
  n_G_q <- as.integer(apply(groups, 2, function(x) length(unique(x))))
  
  # set up grid for output of coefficients
  if (is.null(bins)) {
    bin_vals <- sapply(unique(n_j), function(x) seq(0, 1, length=x))
  } else {
    bin_vals <- bins
  }
  all_bins <- sort(unique(round(unlist(bin_vals), 2)))
  bin_id <- vector("list", length=n)
  for (i in 1:n) {
    bin_vals_temp <- round(seq(0, 1, length=n_j[i]), 2)
    bin_id[[i]] <- match(bin_vals_temp, all_bins)
  }
  grid <- seq(min(all_bins, na.rm=TRUE), max(all_bins, na.rm=TRUE), length=ncol(y))
  endpoints <- c(min(all_bins, na.rm=TRUE) - 1, max(all_bins, na.rm=TRUE) + 1)
  
  # set up priors and hypers
  sigma2_hyper_a <- hypers$psi_main
  sigma2_hyper_b <- hypers$phi_main
  sigma2_gamma_hyper_a <- hypers$psi_gamma
  sigma2_gamma_hyper_b <- hypers$phi_gamma
  beta_hyper <- hypers$sigma2_beta
  
  y <- ifelse(is.na(y), 0, y)
  
  print(dim(y))
  print(dim(x))
  print(dim(groups))
  print(c(n, n_j, n_k, n_p, n_q, n_G_q, n_t))
  
  # run through each chain
  if ((Sys.info()["sysname"] != "Windows") & (par.run)) {
    mod <- mclapply(1:n.chains, splineSliceInternal, y, x, groups, all_bins, degree, n_knots,
                    n.iters, n.burnin, n.thin, inits, n, n_j, n_k, n_p, n_q,
                    n_G_q, n_t, grid, endpoints, sigma2_hyper_a, sigma2_hyper_b,
                    sigma2_gamma_hyper_a, sigma2_gamma_hyper_b, beta_hyper, bin_id,
                    mc.cores=ifelse(n.chains < 4, n.chains, 4))
  } else {
    mod <- lapply(1:n.chains, splineSliceInternal, y, x, groups, all_bins, degree, n_knots,
                  n.iters, n.burnin, n.thin, inits, n, n_j, n_k, n_p, n_q,
                  n_G_q, n_t, grid, endpoints, sigma2_hyper_a, sigma2_hyper_b,
                  sigma2_gamma_hyper_a, sigma2_gamma_hyper_b, beta_hyper, bin_id)
  }
  
  # summarise all outputs from all chains
  n.keep <- dim(mod[[1]][[2]])[3]
  fitted.mean <- apply(array(unlist(lapply(mod, function(x) x[[1]])),
                             dim=c(nrow(mod[[1]][[1]]), ncol(mod[[1]][[1]]), n.keep, n.chains)),
                       c(1, 2), mean)
  fitted.sd <- apply(array(unlist(lapply(mod, function(x) x[[1]])),
                           dim=c(nrow(mod[[1]][[1]]), ncol(mod[[1]][[1]]), n.keep, n.chains)),
                     c(1, 2), sd)
  coefs.mean <- apply(array(unlist(lapply(mod, function(x) x[[2]])),
                            dim=c(nrow(mod[[1]][[2]]), ncol(mod[[1]][[2]]), n.keep, n.chains)),
                      c(1, 2), mean)
  coefs.sd <- apply(array(unlist(lapply(mod, function(x) x[[2]])),
                          dim=c(nrow(mod[[1]][[2]]), ncol(mod[[1]][[2]]), n.keep, n.chains)),
                    c(1, 2), sd)
  rand.coefs.mean <- vector("list", length=n_q)
  rand.coefs.sd <- vector("list", length=n_q)
  gamma.mean <- vector("list", length=n_q)
  gamma.sd <- vector("list", length=n_q)
  for (q in 1:n_q) {
    rand.coefs.mean[[q]] <- apply(array(unlist(lapply(mod, function(x) x[[4]][[q]])),
                                        dim=c(nrow(mod[[1]][[4]][[q]]), ncol(mod[[1]][[4]][[q]]), n.keep, n.chains)),
                                  c(1, 2), mean)
    rand.coefs.sd[[q]] <- apply(array(unlist(lapply(mod, function(x) x[[4]][[q]])),
                                      dim=c(nrow(mod[[1]][[4]][[q]]), ncol(mod[[1]][[4]][[q]]), n.keep, n.chains)),
                                c(1, 2), sd)
    gamma.mean[[q]] <- vector("list", length=n_G_q[q])
    gamma.sd[[q]] <- vector("list", length=n_G_q[q])
    for (qq in 1:n_G_q[q]) {
      gamma.mean[[q]][[qq]] <- apply(array(unlist(lapply(mod, function(x) x[[5]][[q]][qq, , ])),
                                           dim=c(nrow(mod[[1]][[5]][[q]][qq, , ]), n.keep, n.chains)),
                                     1, mean)
      gamma.sd[[q]][[qq]] <- apply(array(unlist(lapply(mod, function(x) x[[5]][[q]][qq, , ])),
                                         dim=c(nrow(mod[[1]][[5]][[q]][qq, , ]), n.keep, n.chains)),
                                   1, sd)
    }
  }
  sigma2.mean <- mean(unlist(lapply(mod, function(x) x[[7]])))
  sigma2.sd <- sd(unlist(lapply(mod, function(x) x[[7]])))
  beta.mean <- vector("list", length=n_k)
  beta.sd <- vector("list", length=n_k)
  for (k in 1:n_k) {
    beta.mean[[k]] <- apply(array(unlist(lapply(mod, function(x) x[[8]][k, , ])),
                                  dim=c(nrow(mod[[1]][[8]][k, , ]), n.keep, n.chains)),
                            1, mean)
    beta.sd[[k]] <- apply(array(unlist(lapply(mod, function(x) x[[8]][k, , ])),
                                dim=c(nrow(mod[[1]][[8]][k, , ]), n.keep, n.chains)),
                          1, sd)
  }
  rho.mean <- mean(unlist(lapply(mod, function(x) x[[6]])))
  rho.sd <- sd(unlist(lapply(mod, function(x) x[[6]])))
  sigma2_gamma.mean <- apply(array(unlist(lapply(mod, function(x) x[[9]])),
                                   dim=c(nrow(mod[[1]][[9]]), n.keep, n.chains)),
                             1, mean)
  sigma2_gamma.sd <- apply(array(unlist(lapply(mod, function(x) x[[9]])),
                                 dim=c(nrow(mod[[1]][[9]]), n.keep, n.chains)),
                           1, sd)
  llik_all <- lapply(mod, function(x) x[[10]])
  
  # summarise all outputs from each chain
  if (n.chains > 1) {
    coefs.chain.var <- apply(array(unlist(lapply(mod, function(x) apply(x[[2]], c(1, 2), var))),
                                   dim=c(nrow(mod[[1]][[2]]), ncol(mod[[1]][[2]]), n.chains)),
                             c(1, 2), mean)
    coefs.chain.var2 <- n.iters * apply(array(unlist(lapply(mod, function(x) apply(x[[2]], c(1, 2), mean))),
                                              dim=c(nrow(mod[[1]][[2]]), ncol(mod[[1]][[2]]), n.chains)),
                                        c(1, 2), var)
    rand.coefs.chain.var <- vector("list", length=n_q)
    gamma.chain.var <- vector("list", length=n_q)
    rand.coefs.chain.var2 <- vector("list", length=n_q)
    gamma.chain.var2 <- vector("list", length=n_q)
    for (q in 1:n_q) {
      rand.coefs.chain.var[[q]] <- apply(array(unlist(lapply(mod, function(x) apply(x[[4]][[q]], c(1, 2), var))),
                                               dim=c(nrow(mod[[1]][[4]][[q]]), ncol(mod[[1]][[4]][[q]]), n.chains)),
                                         c(1, 2), mean)
      gamma.chain.var[[q]] <- vector("list", length=n_G_q[q])
      rand.coefs.chain.var2[[q]] <- n.iters * apply(array(unlist(lapply(mod, function(x) apply(x[[4]][[q]], c(1, 2), mean))),
                                                          dim=c(nrow(mod[[1]][[4]][[q]]), ncol(mod[[1]][[4]][[q]]), n.chains)),
                                                    c(1, 2), var)
      gamma.chain.var2[[q]] <- vector("list", length=n_G_q[q])
      for (qq in 1:n_G_q[q]) {
        gamma.chain.var[[q]][[qq]] <- apply(array(unlist(lapply(mod, function(x) apply(x[[5]][[q]][qq, ,], 1, var))),
                                                  dim=c(nrow(mod[[1]][[5]][[q]][qq, ,]), n.chains)),
                                            1, mean)
        gamma.chain.var2[[q]][[qq]] <- n.iters * apply(array(unlist(lapply(mod, function(x) apply(x[[5]][[q]][qq, , ], 1, mean))),
                                                             dim=c(nrow(mod[[1]][[5]][[q]][qq, , ]), n.chains)),
                                                       1, var)
      }
    }
    sigma2.chain.var <- mean(array(unlist(lapply(mod, function(x) var(x[[7]]))),
                                   dim=c(n.chains)))
    beta.chain.var <- vector("list", length=n_k)
    sigma2.chain.var2 <- n.iters * var(array(unlist(lapply(mod, function(x) mean(x[[7]]))),
                                             dim=c(n.chains)))
    beta.chain.var2 <- vector("list", length=n_k)
    for (k in 1:n_k) {
      beta.chain.var[[k]] <- apply(array(unlist(lapply(mod, function(x) apply(x[[8]][k, ,], 1, var))),
                                         dim=c(nrow(mod[[1]][[8]][k, ,]), n.chains)),
                                   1, mean)
      beta.chain.var2[[k]] <- n.iters * apply(array(unlist(lapply(mod, function(x) apply(x[[8]][k, , ], 1, mean))),
                                                    dim=c(nrow(mod[[1]][[8]][k, , ]), n.chains)),
                                              1, var)
    }
    rho.chain.var <- mean(array(unlist(lapply(mod, function(x) var(x[[6]]))),
                                dim=c(n.chains)))
    rho.chain.var2 <- n.iters * var(array(unlist(lapply(mod, function(x) mean(x[[6]]))),
                                          dim=c(n.chains)))
    sigma2_gamma.chain.var <- apply(array(unlist(lapply(mod, function(x) apply(x[[9]], 1, var))),
                                          dim=c(nrow(mod[[1]][[9]]), n.chains)),
                                    1, mean)
    sigma2_gamma.chain.var2 <- n.iters * apply(array(unlist(lapply(mod, function(x) apply(x[[9]], 1, mean))),
                                                     dim=c(nrow(mod[[1]][[9]]), n.chains)),
                                               1, var)
    
    # calculate Rhats
    rhat.calc <- function(chain.var, chain.var2, n.its) {
      rhat <- sqrt(((1 - (1 / n.its)) * chain.var + (1 / n.its) * chain.var2) / chain.var)
      return(rhat)
    }
    coefs.rhat <- rhat.calc(coefs.chain.var, coefs.chain.var2, n.iters)
    rand.coefs.rhat <- vector("list", length=n_q)
    gamma.rhat <- vector("list", length=n_q)
    for (q in 1:n_q) {
      rand.coefs.rhat[[q]] <- rhat.calc(rand.coefs.chain.var[[q]], rand.coefs.chain.var2[[q]], n.iters)
      gamma.rhat[[q]] <- vector("list", length=n_G_q[q])
      for (qq in 1:n_G_q[q]) {
        gamma.rhat[[q]][[qq]] <- rhat.calc(gamma.chain.var[[q]][[qq]], gamma.chain.var2[[q]][[qq]], n.iters)
      }
    }
    sigma2.rhat <- rhat.calc(sigma2.chain.var, sigma2.chain.var2, n.iters)
    beta.rhat <- vector("list", length=n_k)
    for (k in 1:n_k) {
      beta.rhat[[k]] <- rhat.calc(beta.chain.var[[k]], beta.chain.var2[[k]], n.iters)
    }
    rho.rhat <- rhat.calc(rho.chain.var, rho.chain.var2, n.iters)
    sigma2_gamma.rhat <- rhat.calc(sigma2_gamma.chain.var, sigma2_gamma.chain.var2, n.iters)
  } else {
    coefs.rhat <- NULL
    rand.coefs.rhat <- NULL
    gamma.rhat <- NULL
    sigma2.rhat <- NULL
    beta.rhat <- NULL
    rho.rhat <- NULL
    sigma2_gamma.rhat <- NULL
  }
  rhats <- list(coefs=coefs.rhat, rand.coefs=rand.coefs.rhat,
                gamma=gamma.rhat, sigma2=sigma2.rhat, beta=beta.rhat,
                rho=rho.rhat, sigma2_gamma=sigma2_gamma.rhat)
  
  # calculate DIC from mean log likelihood and log likelihood at mean(parameters)
  loglik.mean <- mean(unlist(lapply(mod, function(x) x[[3]])))
  theta1 <- seq(min(all_bins, na.rm=TRUE) - 0.5, max(all_bins, na.rm=TRUE) + 0.5,
                length=n_knots)
  b_splines_mat <- calc_bs(all_bins, theta1, degree,
                           c(min(all_bins, na.rm=TRUE) - 1, max(all_bins, na.rm=TRUE) + 1))
  beta.mean2 <- matrix(unlist(beta.mean), ncol=n_p, byrow=TRUE)
  gamma.mean2 <- vector("list", length=n_q)
  for (i in seq(along=gamma.mean2)) {
    gamma.mean2[[i]] <- matrix(unlist(gamma.mean[[i]]), ncol=n_t, byrow=TRUE)
  }
  loglik.param.bar <- lnL(y=y, x=x, groups=groups, beta=beta.mean2, gamma=gamma.mean2,
                          sigma2=sigma2.mean, rho=rho.mean, b_splines_mat=b_splines_mat,
                          beta_hyper, n, n_j, n_k, n_q, bin_id)
  DIC <- -2 * loglik.mean + 2 * loglik.param.bar
  
  # calculate all returned values
  family <- "gaussian"
  y.tmp <- NULL
  fitted.tmp <- NULL
  for (i in 1:n) {
    y.tmp <- c(y.tmp, y[i, 1:n_j[i]])
    fitted.tmp <- c(fitted.tmp, fitted.mean[i, 1:n_j[i]])
  }
  r2 <- cor(y.tmp, fitted.tmp)
  return(list(fitted=fitted.mean, fitted.sd=fitted.sd, observed=y,
              coefs.mean=coefs.mean, coefs.sd=coefs.sd,
              rand.coefs.mean=rand.coefs.mean, rand.coefs.sd=rand.coefs.sd,
              r2=r2, family=family, DIC=DIC, rhats=rhats,
              sigma2.mean=sigma2.mean, sigma2.sd=sigma2.sd,
              sigma2_gamma.mean=sigma2_gamma.mean,
              sigma2_gamma.sd=sigma2_gamma.sd,
              beta.mean=beta.mean, gamma.mean=gamma.mean,
              rho.mean=rho.mean, rho.sd=rho.sd,
              llik_all=llik_all))
}
