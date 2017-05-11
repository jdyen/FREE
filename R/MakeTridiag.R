MakeTridiag <-
function(diagon, upper, lower, nrow=NULL){
  if (is.null(nrow) & {length(diagon) == 1}) {
    stop("diagon must have length greater than 1 if nrow is not set.....", call.=FALSE)
  }
  if ({length(diagon) == 1} & !is.null(nrow)){
    diagon <- rep(diagon, nrow)
  } else {
    nrow <- length(diagon)
  }
  if ( {length(upper) != {nrow - 1}} & length(upper) != 1){
    warning("upper should either be a single value or a vector of length nrow-1.....",
            call.=FALSE)
  }
  if ( {length(lower) != {nrow - 1}} & length(lower) != 1){
    warning("lower should either be a single value or a vector of length nrow-1.....",
            call.=FALSE)
  }
  ncol <- nrow
  mat <- diag(diagon, nrow, ncol)
  R <- row(mat)
  C <- col(mat)
  mat[C == {R + 1}] <- upper
  mat[C == {R - 1}] <- lower
  return(mat)
}
