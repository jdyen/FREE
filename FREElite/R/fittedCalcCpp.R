# fitted_calc
#
# Description: function to calculate fitted values for a function reegression model
#   based on B-spline basis functions
#
# Arguments
#
# x           predictor variables
# groups      group IDs for clustering variables
# beta        coefficients for main predictors
# gamma       random curves for group variables
# b_splines_mat
# n           number of subjects
# n_j         number of observations for each subject
# n_q         number of grouping variables
# bin_id

fitted_calc <- function(x, groups, beta, gamma, b_splines_mat,
n, n_j, n_q, bin_id) {
  out <- matrix(NA, nrow=n, ncol=max(n_j, na.rm=TRUE))
  for (i in 1:n) {
    out[i, 1:n_j[i]] <- x[i, ] %*% t(b_splines_mat[bin_id[[i]], ] %*% t(beta))
    for (q in 1:n_q) {
      out[i, 1:n_j[i]] <- t(b_splines_mat[bin_id[[i]], ] %*% gamma[[q]][groups[i, q], ])
    }
  }
  return(out)
}