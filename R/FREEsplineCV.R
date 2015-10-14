FREEsplineCV <- function(n.cv, y, x, groups, w, degree=3, n_knots_beta=11,
                         n_knots_gamma=11, n.iters, n.chains=1, n.thin=1,
                          hypers=list(psi_main=0.1, phi_main=0.1, psi_gamma=0.1,
                                      phi_gamma=0.1, sigma2_beta=10), ...) {
  mod <- vector("list", length=n.cv)
  for (i in 1:n.cv) {
    mod[[i]] <- cv_inner(i, n.cv=n.cv, y=y, x=x,
                         groups=groups, w=w, degree=degree,
                         n_knots_beta=n_knots_beta, n_knots_gamma=n_knots_gamma,
                         n.iters=n.iters, n.chains=n.chains, n.thin=n.thin,
                         hypers=hypers, ...)
  }
  pred <- unlist(lapply(mod, function(x) x$pred))
  obs <- unlist(lapply(mod, function(x) x$y))
  return(list(predicted=pred, observed=obs))
}

cv_inner <- function(i, n.cv, y, x, groups, w, degree=3, n_knots_beta=11,
                     n_knots_gamma=11, n.iters, n.chains=1, n.thin=1,
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
  groups.use <- matrix(as.integer(as.factor(groups[-subset, ])), ncol=ncol(groups))
  groups.use <- apply(groups.use, 2, function(x) as.integer(as.factor(x)))
  mod <- FREEspline(y=y.use, x=x.use, groups=groups.use, w=w, degree=degree,
                    n_knots_beta=n_knots_beta, n_knots_gamma=n_knots_gamma,
                    n.iters=n.iters, n.chains=n.chains, n.thin=n.thin,
                    hypers=hypers, ...)
  x.tmp <- cbind(rep(1, length(subset)), x[subset, ])
  mod.pred <- as.matrix(x.tmp) %*% mod$coefs.mean
  return(list(pred=mod.pred, y=y[subset, ]))
}
