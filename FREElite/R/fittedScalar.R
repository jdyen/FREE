fitted_scalar <- function(x, z, groups, beta, gamma, delta, bs_beta) {
  fitted <- (x %*% bs_beta) %*% beta
  if (length(delta) > 1) {
    fitted <- fitted + (z %*% delta)
  } else {
    fitted <- fitted + (z * delta)
  }
  for (q in seq(along=gamma)) {
    fitted <- fitted + gamma[[q]][groups[, q] + 1]
  }
  return(fitted)
}

coefs_calc_scalar <- function(beta, theta, degree, grid, endpoints) {
  bs_beta <- calc_bs_scalar(grid, theta, degree, endpoints)
  out <- c(bs_beta %*% beta)
  return(out)
}

fitted_scalar_cv <- function(x, z, beta, delta, bs_beta) {
  fitted <- (x %*% bs_beta) %*% beta
  if (length(delta) > 1) {
    fitted <- fitted + (z %*% delta)
  } else {
    fitted <- fitted + (z * delta)
  }
  return(fitted)
}