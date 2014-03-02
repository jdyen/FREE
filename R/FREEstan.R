FREEstan <-
function(y, x, bins, stan.file=NULL, stan.model=NA, Kt=12, iid.er=FALSE, n.chains=3, n.iters=2000, n.burnin=n.iters/2, n.thin=1, verbose=FALSE, refresh=max(n.iters/10, 1), ...){
  if (get_cppo()$mode != "fast") {
    set_cppo("fast")
  }
  if (is.null(stan.file) & is.na(stan.model)) {
    stan.file <- stanCodeDefault()
  }
  D <- ncol(y)
  Y <- y
  N <- nrow(Y)
  data.list <- NULL
  data.list <- rep(1,N)
  for (i in 1:ncol(x)) {
    data.list <- cbind(data.list, x[1:N, i])
  }
  colnames(data.list) <- c("int", colnames(x))
  k <- ncol(data.list)
  Kt <- Kt
  grid <- seq(min(bins), max(bins), length=D)
  BS <- bs(1:D, df=Kt, intercept=TRUE, degree=3)
  if (!iid.er) {
    alpha <- 0.1
    diff0 <- diag(1, D, D)
    diff2 <- matrix(rep(c(1, -2, 1, rep(0, D - 2)), D - 2)[1:{{D - 2} * D}], D - 2, D, byrow=T)
    P0 <- t(BS) %*% t(diff0) %*% diff0 %*% BS
    P2 <- t(BS) %*% t(diff2) %*% diff2 %*% BS
    CovMat <- alpha * P0 * {1 - alpha} * P2
    CovMat <- round(CovMat, 4)
  } else {
    CovMat <- round(diag(1, Kt), 4)
  }
  dat <- list(Y=Y, X=data.list, N=N, D=D, k=k, Kt=Kt, BS=BS[1:D,1:Kt], CovMat=CovMat)
  if (is.na(stan.model)) {
    model <- stan(model_code=stan.file, fit=stan.model, data=dat, chains=n.chains, iter=n.iters,
                  warmup=n.burnin, verbose=verbose, thin=n.thin, refresh=refresh, ...)
  } else {
    model <- stan(fit=stan.model, data=dat, chains=n.chains, iter=n.iters, warmup=n.burnin,
                  verbose=verbose, thin=n.thin, refresh=refresh, ...)
  }
  fitted.betas.bs <- matrix(monitor(as.array(model), print=FALSE)
                            [1:{ncol(data.list) * Kt}, "mean"], ncol=ncol(data.list),
                            byrow=TRUE)
  beta.upper.bs <- matrix(monitor(as.array(model), print=FALSE)
                          [1:{ncol(data.list) * Kt}, "97.5%"], ncol=ncol(data.list),
                          byrow=TRUE)
  beta.lower.bs <- matrix(monitor(as.array(model), print=FALSE)
                          [1:{ncol(data.list) * Kt}, "2.5%"], ncol=ncol(data.list), byrow=TRUE)
  fitted.betas <- BS %*% fitted.betas.bs
  beta.upper <- BS %*% beta.upper.bs
  beta.lower <- BS %*% beta.lower.bs
  beta.sd.TEMP <- beta.upper - fitted.betas
  colnames(fitted.betas) <- c("INTERCEPT", colnames(x))
  colnames(beta.sd.TEMP) <- colnames(fitted.betas)
  Rhats <- monitor(as.array(model), print=FALSE)[,"Rhat"]
  if (any(Rhats > 1.1)) {
    warning("Some Rhats are greater than 1.1, so convergence is not guaranteed...", 
            call.=FALSE)
  }
  fitted.y.mean <- as.matrix(data.list) %*% t(fitted.betas)
  fitted.y.upper <- as.matrix(data.list) %*% t(beta.upper)
  fitted.y.lower <- as.matrix(data.list) %*% t(beta.lower)
  r <- cor(c(fitted.y.mean), c(y))
  r2 <- r * r
  xIC <- NULL
  family <- "Gaussian"
  return(list(fitted=fitted.y.mean, observed=y, coefs.mean=t(fitted.betas),
              coefs.sd=t(beta.sd.TEMP), r2=r2, family=family, bins=bins, xIC=xIC,
              stan.model=model))
}