FREEbugsJump <-
function(y, x, bins, family="gaussian", errors="ar1", order=3, cont=order, bugs.file=NULL, n.chains=3, n.iters=5000, n.burnin=n.iters/2, n.thin=1, debug=FALSE, bugs.dir=NULL, ...){
  if (is.null(bugs.file)) {
    bugs.file <- "FREEbugsJumpTemp.txt"
    bugs.switch <- TRUE
  } else {
    bugs.switch <- FALSE
  }
  if (is.null(bugs.dir)) {
    bugs.dir <- "c:/Program Files/WinBUGS14/"
  }
  Nsites <- nrow(y)
  Nclasses <- ncol(y)
  response <- y
  predictors <- x
  predictors <- cbind(rep(1, Nsites), predictors)
  Q <- ncol(predictors)
  if (bugs.switch) {
    MakeBUGSfile(Q=Q, filename="FREEbugsJumpTemp.txt", order=order, cont=cont, ARmod={errors == "ar1"})
  }
  size.class <- bins
  maxsd <- c(5, rep(10, Q), 5)
  bugdata <- list("Nsites", "Nclasses", "size.class", "response", "maxsd", "Q", "predictors")
  if (errors == "ar1") {
    inits <- function(){
      list(alpha0=0, rho=rep(0.01, Nsites), rho.m=0.01, k=rep(0, Q + 2),
           sd=maxsd, site.e=rnorm(Nsites), beta=rep(0, Q), sd.e=rep(1,3))
    }
  } else {
    inits <- function(){
      list(alpha0=0, k=rep(0, Q + 2), sd=maxsd, site.e=rnorm(Nsites), beta=rep(0, Q), sd.e=rep(1,3))
    }
  }
  parameters <- c("mu", "beta.fun")
  model <- bugs(data=bugdata, inits=inits, model.file=bugs.file,
                parameters.to.save=parameters, n.chains=n.chains, n.iter=n.iters,
                n.burnin=n.burnin, n.thin=n.thin, debug=debug, bugs.directory=bugs.dir,
                ...)
  fitted.y.mean <- model$mean$mu
  fitted.beta.mean <- model$mean$beta.fun
  fitted.beta.sd <- model$sd$beta.fun
  rownames(fitted.beta.mean) <- c("INTERCEPT", colnames(x))
  rownames(fitted.beta.sd) <- colnames(fitted.beta.mean)
  if (any(colnames(model$summary) == "Rhat")) {
    Rhats <- model$summary[,"Rhat"]
    if (any(Rhats > 1.1)) {
      warning("Some Rhats are greater than 1.1, so convergence is not guaranteed...",
              call.=FALSE)
    }
  }
  r <- cor(c(fitted.y.mean), c(y))
  r2 <- r * r
  xIC <- model$DIC
  return(list(fitted=fitted.y.mean, observed=y, coefs.mean=fitted.beta.mean, coefs.sd=fitted.beta.sd, r2=r2, family=family, bins=bins, xIC=xIC))
}
