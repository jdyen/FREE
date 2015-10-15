MakeCovMat <-
function(n=10){
  decay <- exp(c(-1:-n)/3)
  CovMat <- matrix(0, nrow=n, ncol=n)
  for (i in 1:{nrow(CovMat) - 1}) {
    CovMat[row(CovMat) == i][{i + 1}:nrow(CovMat)] <- decay[1:{length(decay) - i}]
  }
  diag(CovMat) <- 1
  CovMat[lower.tri(CovMat)] <- t(CovMat)[lower.tri(CovMat)]
  CovMat
}
