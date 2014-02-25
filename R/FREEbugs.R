FREEbugs <-
function(y, x, bins, bugs.file=NULL, Kt=12, iid.er=FALSE, n.chains=3, n.iters=2000, n.burnin=n.iters/2, n.thin=1, debug=FALSE, bugs.dir=NULL){
  if (is.null(bugs.file)) {
    bugs.file <- "FREEbugsSplineTemp.txt"
    MakeBUGSsplineFile(filename="FREEbugsSplineTemp.txt", ARmod=!iid.er)
  }
  if (is.null(bugs.dir)) {
    bugs.dir <- "C:/Users/jdyen/Documents/WinBUGS14"
  }
  N <- nrow(y)
  p <- ncol(y)
  X <- t(x)
  X <- rbind(rep(1, N), X)
  v <- nrow(X)
  BS <- bs(x=bins, df=Kt, degree=3, intercept=TRUE)
  BS <- BS[1:p, 1:Kt]
  if (iid.er == FALSE) {
    alpha <- 0.1
    diff0 <- diag(1, p, p)
    diff2 <- matrix(rep(c(1, -2, 1, rep(0,p-2)), p-2)[1:{{p-2}*p}], p-2, p, byrow=T)
    P0 <- t(BS) %*% t(diff0) %*% diff0 %*% BS
    P2 <- t(BS) %*% t(diff2) %*% diff2 %*% BS
    CovMat <- alpha * P0 * {1-alpha} * P2
    CovMat <- round(CovMat, 4)
  } else {
    CovMat <- round(diag(1, Kt), 4)
  }
  save.params <- c("mu", "beta")
  bugdata <- list("y", "X", "N", "p", "v", "BS", "Kt", "CovMat")
  if (iid.er) {
    initials <- function(){
      list(tau_y=1, tau_beta=rep(1, v), beta=matrix(rnorm(Kt * v), ncol=Kt), tau_site=1,
           site.e=rep(0, N))
    }
  } else {
    initials <- function(){
      list(tau_y=1, tau_beta=rep(1, v), beta=matrix(rnorm(Kt * v), ncol=Kt), tau_site=1,
           site.e=rep(0, N), rho.m=0.01, rho=rep(0.01, N))
    }
  } 
  model <- bugs(data=bugdata, inits=initials, model.file=bugs.file,
                parameters.to.save=save.params, n.chains=n.chains,
                n.iter=n.iters, n.burnin=n.burnin, n.thin=n.thin, debug=debug,
                bugs.directory=bugs.dir)
  fitted.y.mean <- model$mean$mu
  fitted.beta.bs <- t(model$mean$beta)
  fitted.beta.sd <- t(model$sd$beta)
  fitted.betas <- BS %*% fitted.beta.bs
  beta.sd.TEMP <- BS %*% fitted.beta.sd
  colnames(fitted.betas) <- c("INTERCEPT", colnames(x))
  colnames(beta.sd.TEMP) <- colnames(fitted.betas)
  if (any(colnames(model$summary) == "Rhat")) {
    Rhats <- model$summary[,"Rhat"]
    if (any(Rhats > 1.1)) {
      warning("Some Rhats are greater than 1.1, so convergence is not guaranteed...", call.=FALSE)
    }
  }
  r <- cor(c(fitted.y.mean), c(y))
  r2 <- r * r
  xIC <- model$DIC
  family <- "Gaussian"
  return(list(fitted=fitted.y.mean, observed=y, coefs.mean=t(fitted.betas), coefs.sd=t(beta.sd.TEMP), r2=r2, family=family, bins=bins, xIC=xIC))
}
