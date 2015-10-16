FREEscalarCV <- function(n.cv=10, y, x, z, groups, bins, 
                         degree=3, n_knots=8, n.iters=100,
                         n.burnin=round(n.iters / 5), n.thin=1, n.chains=3,
                         hypers=list(phi1=0.1, psi1=0.1, phi2=0.1, psi2=0.1, s2_alpha=10,
                                     s2_beta=10, s2_delta=10),
                         inits=list(alpha=NULL, beta=NULL, gamma=NULL, delta=NULL,
                                    sigma2=1, sigma2_gamma=1),
                         par.run=FALSE, verbose=TRUE, ...)
{
  pred <- NULL
  n.out <- floor(length(y) / n.cv)
  theta <- seq(5, ncol(x) - 4, length=n_knots - degree)
  bs_beta <- calc_bs(1:ncol(x), theta, degree, c(0, ncol(x) + 1))
  for (i in 1:n.cv) {
    if (verbose) {
      print(paste("starting fold ", i, " of ", n.cv, "...", sep=""))
      flush.console()
    }
    if (i < n.cv) {
      subset <- seq((i - 1) * n.out + 1, i * n.out, by=1)
    } else {
      subset <- seq((i - 1) * n.out + 1, length(y), by=1)
    }
    if (is.matrix(groups[-subset, ])) {
      groups.use <- matrix(0, nrow=nrow(groups[-subset, ]), ncol=ncol(groups))
    } else {
      groups.use <- matrix(0, nrow=nrow(as.matrix(groups[-subset, ])), ncol=ncol(groups))
    }
    for (j in 1:ncol(groups)) {
      groups.use[, j] <- as.integer(as.factor(groups[-subset, j]))
    }
    mod <- FREEscalar(y=y[-subset], x=x[-subset, ], z=z[-subset, ],
                      groups=groups.use, bins=bins,
                      degree=degree, n_knots=n_knots, n.iters=n.iters,
                      n.burnin=n.burnin, n.thin=n.thin, n.chains=n.chains,
                      hypers=hypers, inits=inits, par.run=par.run, ...)
    pred <- c(pred, fitted_scalar_cv(x=x[subset, ], z=z[subset, ],
                                     beta=mod$beta.mean, delta=mod$delta.mean,
                                     bs_beta=bs_beta))
    if (verbose) {
      print(paste("finished fold ", i, " of ", n.cv, "...", sep=""))
      flush.console()
    }
  }
  r2 <- cor(y, pred) ** 2
  return(list(predicted=pred, observed=y))
}