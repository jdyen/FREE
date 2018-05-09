FREEsplineCV <- function(n.cv, y, x, groups, bins, degree=3, n_knots=5,
                         n.iters=100, n.chains=3, n.thin=1,
                         n.burnin=round(n.iters / 5), par.run=FALSE,
                         hypers=list(psi_main=0.1, phi_main=0.1, psi_gamma=0.1,
                                      phi_gamma=0.1, sigma2_beta=10),
                         verbose=TRUE, ...) {
  mod <- vector("list", length=n.cv)
  for (i in 1:n.cv) {
    if (verbose) {
      print(paste("starting fold ", i, " of ", n.cv, "...", sep=""))
      flush.console()
    }
    mod[[i]] <- cv_inner(i, n.cv=n.cv, y=y, x=x,
                         groups=groups, bins=bins, degree=degree,
                         n_knots=n_knots, 
                         n.iters=n.iters, n.chains=n.chains, n.thin=n.thin,
                         n.burnin=n.burnin, par.run=par.run,
                         hypers=hypers, ...)
    if (verbose) {
      print(paste("finished fold ", i, " of ", n.cv, "...", sep=""))
      flush.console()
    }
  }
  pred <- unlist(lapply(mod, function(x) x$pred))
  obs <- unlist(lapply(mod, function(x) x$y))
  return(list(predicted=pred, observed=obs))
}

cv_inner <- function(i, n.cv, y, x, groups, bins, degree=3, n_knots=5,
                     n.iters, n.chains, n.thin,
                     n.burnin, par.run,
                     hypers=list(psi_main=0.1, phi_main=0.1, psi_gamma=0.1,
                                 phi_gamma=0.1, sigma2_beta=10), ...)
{
  l.out <- floor(nrow(y) / n.cv)
  if (i < n.cv) {
    subset <- seq((i - 1) * l.out + 1, i * l.out, by=1)
  } else {
  	subset <- seq((i - 1) * l.out + 1, nrow(y), by=1)
  }
  y.use <- y[-subset, ]
  x.use <- x[-subset, ]
  if (!is.null(groups)) {
    groups.use <- matrix(as.integer(as.factor(groups[-subset, ])), ncol=ncol(groups))
    groups.use <- apply(groups.use, 2, function(x) as.integer(as.factor(x)))
  } else {
    groups.use <- NULL
  }
  mod <- FREEspline(y=y.use, x=x.use, groups=groups.use, bins=NULL, degree=degree,
                    n_knots=n_knots,
                    n.iters=n.iters, n.burnin=n.burnin, n.thin=n.thin,
                    n.chains=n.chains,
                    hypers=hypers, par.run=par.run, ...)
  x.tmp <- cbind(rep(1, length(subset)), x[subset, ])
  mod.pred <- as.matrix(x.tmp) %*% mod$coefs.mean
  return(list(pred=mod.pred, y=y[subset, ]))
}
