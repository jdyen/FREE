MakeBUGSfile <-
function(Q, filename="FREEbugsJumpTemp.txt", order=3, cont=order, ARmod=FALSE){
  sink(filename, append=FALSE)
  if (ARmod) {
    cat("
    model{
      for (i in 1:Nsites) {
        for (j in 1:Nclasses) {
          response[i, j] ~ dnorm(mu2[i, j], tau[1])

          mu2[i, j] <- mu[i, j] + cor.e[i, j]

          mu[i, j] <- alpha.fun[j] + sum(funs[i, j, ]) + site.e[i]

          for (v in 1:Q) {
            funs[i, j, v] <- beta.fun[v, j] * Xc[i, v]
          }

          resid[i, j] <- response[i, j] - mu[i, j]     
        }

        site.e[i] ~ dnorm(0, tau[2])
        cor.e[i, 1] <- 0
        rho[i] ~ dnorm(rho.m, 10)

        for (j in 2:Nclasses) {
          cor.e[i, j] <- rho[i] * resid[i, (j - 1)]
        }", fill=TRUE)
  } else {
    cat("

      model{

      for (i in 1:Nsites) {
        for (j in 1:Nclasses) {
          response[i, j] ~ dnorm(mu[i, j], tau[1])

          mu[i, j] <- alpha.fun[j] + sum(funs[i, j, ]) + site.e[i]

          for(v in 1:Q) {
            funs[i, j, v] <- beta.fun[v, j] * Xc[i, v]
          }

          resid[i, j] <- response[i, j] - mu[i, j]     
        }

        site.e[i]~dnorm(0,tau[2])", fill=T)
  }


  cat("
        for (v in 1:Q) {
          Xc[i, v] <- cut(predictors[i, v])
          predictors[i, v] ~ dnorm(0, 1)
        }
      }  

      alpha0 ~ dnorm(0, 0.001)

      tau.beta.int <- 1 / pow(maxsd[Q + 2], 2)

      for (v in 1:Q) {
        beta[v] ~ dnorm(0, tau.beta.int)
      }

      for (j in 1:Nclasses) {
        alpha.fun[j] <- alpha0 + alpha.fun0[j] - alpha.fun0[1]
", fill=TRUE)

  for (v in 1:Q) {
    cat(paste("        beta.fun[", v, ", j] <- beta[", v, "] + beta.fun", v, "[j] - beta.fun", v, "[1]", sep=""), "\n")
  }

  cat("
      }

      alpha.fun0[1:Nclasses] <- jump.pw.poly.df.cub(size.class[1:Nclasses], k[1], tau.beta[1])


      for (i in 1:(Q + 2)) {
        k[i] ~ dbin(0.5, kmax)
        tau.beta[i] <- 1 / pow(sd[i], 2)
        sd[i] ~ dunif(0, maxsd[i])
      }

      tau[1] <- 1 / pow(sd.e[1], 2)
      kmax <- Nclasses / 2

      for (t in 2:3) {
        tau[t] <- 1 / pow(sd.e[t], 2)
        sd.e[t] ~ dunif(0, maxsd[(t - 1)])
      }
      sd.e[1] ~ dunif(0, 5)
", fill=TRUE)

  for (v in 1:Q) {
    cat(paste("      beta.fun", v, "[1:Nclasses] <- jump.pw.poly.df.gen(size.class[1:Nclasses], k[", v + 1, "], tau.beta[", v + 1, "], ", order, ",", cont, ")", sep=""), "\n")
  }

  if (ARmod) {
    cat("      rho.m ~ dunif(-1, 1)", "\n")
  }

  cat("    }")
  sink(file=NULL)
}
