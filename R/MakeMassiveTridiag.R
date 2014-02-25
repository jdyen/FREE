MakeMassiveTridiag <-
function(tridiag, n.times){
  # Function to make a pseudo-diagonal matrix by repeating one smaller matrix along
  #   the main diagonal of a larger matrix
  #
  # Args:
  #   blah
  #
  # Returns:
  #   blah
  mat0 <- matrix(0, nrow=nrow(tridiag), ncol=ncol(tridiag))
  mat.temp <- lapply(1:n.times, 
                     function(i, mat0, tridiag){
                       mat.out <- rbind(do.call("rbind", rep(list(mat0), {i - 1})),
                       tridiag, do.call("rbind", rep(list(mat0), {n.times - i})))
                       return(mat.out)
                     },
                     mat0=mat0, tridiag=tridiag)
  mat <- do.call("cbind", mat.temp)
  mat
}
