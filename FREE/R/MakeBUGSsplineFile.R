MakeBUGSsplineFile <-
function(filename="FREEbugsSplineTemp.txt", ARmod=FALSE){
  sink(filename, append=FALSE)
  if (ARmod) {

cat("model{
  for (i in 1:N) {
    for (j in 1:p) {
      y[i, j] ~ dnorm(mu2[i, j], tau_y)
      mu2[i, j] <- mu[i, j] + cor.e[i, j]
      mu[i, j] <- sum(temp.var2[i, ]) + site.e[i]
      temp.var2[i, j] <- inprod(temp.var[j, ], X[, i])
      resid[i, j] <- y[i, j] - mu[i, j]
    }
    site.e[i] ~ dnorm(0, tau_site)
    cor.e[i, 1] <- 0
    rho[i] ~ dnorm(rho.m, 10)
    for (j in 2:p) {
      cor.e[i, j] <- rho[i] * resid[i, (j - 1)]
    }", fill=TRUE)
  } else {
    cat("model{
  for (i in 1:N) {
    for (j in 1:p) {
      y[i, j] ~ dnorm(mu[i, j], tau_y)
      mu[i, j] <- sum(temp.var2[i, ]) + site.e[i]
      temp.var2[i, j]  <- inprod(temp.var[j, ], X[, i])
    }
    site.e[i] ~ dnorm(0, tau_site)", fill=TRUE)
  }
  cat("
  }
  for (i in 1:p) {
    for (j in 1:v) {
      temp.var[i, j] <- inprod(BS[i, ], beta[j, ])
    }
  }
  for (i in 1:Kt) {
    mu_beta[i] <- 0
  }
  for (i in 1:v) {
    beta[i, 1:Kt] ~ dmnorm(mu_beta[], TauMat[, ])
  }
  tau_y ~ dgamma(0.01, 0.01)
  tau_site ~ dgamma(0.01, 0.01)
  TauMat[1:Kt, 1:Kt] <- inverse(CovMat[, ])
  for (i in 1:v) {
    tau_beta[i] ~ dgamma(0.01, 0.01)
  }
  ", fill=TRUE)
  if (ARmod) {
    cat("  rho.m ~ dunif(-1, 1)", "\n")
  }
  cat("}")
  sink(file=NULL)
}
