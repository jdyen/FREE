installFREE <- function(OSX.install=FALSE) {
  if(!require(fda)) {
    install.packages("fda")
  }
  if(!require(INLA)) {
    source("http://www.math.ntnu.no/inla/givemeINLA.R")
  }
  if(!require(mboost)) {
    install.packages("mboost")
  }
  if(!require(MASS)) {
    install.packages("MASS")
  }
  if(!require(Rcpp)) {
    install.packages("Rcpp")
  }
  if(!require(BayesX)) {
    install.packages("BayesX")
  }
  if(!require(boot)) {
    install.packages("boot")
  }
  if(!require(coda)) {
    install.packages("coda")
  }
  if(!require(inline)) {
    install.packages("inline")
  }
  if(!require(lattice)) {
    install.packages("lattice")
  }
  if(!require(maptools)) {
    install.packages("maptools", repos="http://R-Forge.R-project.org")
  }
  if(!require(Matrix)) {
    install.packages("Matrix")
  }
  if(!require(parallel)) {
    install.packages("parallel")
  }
  if(!require(sp)) {
    install.packages("sp")
  }
  if(!require(survival)) {
    install.packages("survival")
  }
  if (!OSX.install) {
    if(!require(R2WinBUGS)) {
      install.packages("R2WinBUGS")
    }
    if(!require(rstan)) {
      install.packages("rstan", type="source",
      repos = c(getOption("repos"), rstan = "http://wiki.rstan-repo.googlecode.com/git/"),
      INSTALL_opts = "--merge-multiarch")
    }
  }
  install.packages("FREE_1.0.tar.gz", type="source", repos=NULL)
}