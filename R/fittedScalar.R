fitted_scalar <- function(x, z, groups, beta, gamma, delta, bs_beta) {
  fitted <- rep(0, nrow(x[[1]]))
  for (i in 1:length(x)) {
    fitted <- fitted + (x[[i]] %*% bs_beta[[i]]) %*% beta[i, ]
  }
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
  out <- matrix(NA, nrow = nrow(beta), ncol = length(grid))
  bs_beta <- vector('list', length = nrow(beta))
  for (k in 1:nrow(beta)) {
    bs_beta[[k]] <- calc_bs_scalar(grid, theta, degree, endpoints)
    out[k, ] <- c(bs_beta[[k]] %*% beta[k, ])
  }
  return(out)
}

fitted_scalar_cv <- function(x, z, beta, delta, bs_beta) {
  if (!is.null(bs_beta)) {
    fitted <- rep(0, nrow(x[[1]]))
    for (i in 1:length(x)) {
      fitted <- fitted + (x[[i]] %*% bs_beta[[i]]) %*% beta[i, ]
    }
  } else {
    fit_length <- ifelse(length(delta) > 1, nrow(z), length(z))
    print(fit_length)
    fitted <- rep(0, fit_length)
  }
  if (length(delta) > 1) {
    fitted <- fitted + (z %*% delta)
  } else {
    fitted <- fitted + (z * delta)
  }
  return(fitted)
}