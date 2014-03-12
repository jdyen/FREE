MakeMassiveTridiag <-
function(tridiag, n.times){
  mat.temp <- lapply(1:n.times, 
                     function(i, tridiag, n.times){
                       mat.out <- matrix(0, nrow={nrow(tridiag) * n.times},
                                         ncol=ncol(tridiag))
                       mat.out[{{i - 1} * nrow(tridiag) + 1}:{i * nrow(tridiag)},
                                 1:ncol(tridiag)] <-
                                 tridiag
                       return(mat.out)
                     },
                     tridiag=tridiag, n.times=n.times)
  mat <- do.call("cbind", mat.temp)
  mat
}
