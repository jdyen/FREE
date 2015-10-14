installFREE <- function() {
  if(!require(fda)) {
    install.packages("fda")
  }
  if(!require(INLA)) {
    source("http://www.math.ntnu.no/inla/givemeINLA.R")
  }
  if(!require(MASS)) {
    install.packages("MASS")
  }
  if(!require(mboost)) {
    install.packages("mboost")
  }
  if(!require(parallel)) {
    install.packages("parallel")
  }
  if(!require(R2WinBUGS)) {
    install.packages("R2WinBUGS")
  }
  if(!require(Rcpp)) {
    install.packages("Rcpp")
  }
  if(!require(rstan)) {
    install.packages("rstan", type="source",
    repos = c(getOption("repos"), rstan = "http://wiki.rstan-repo.googlecode.com/git/"),
    INSTALL_opts = "--merge-multiarch")
  }
  install.packages("FREE_2.0.tar.gz", type="source", repos=NULL)
}