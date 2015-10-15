# fitted_calc
#
# Description: function to calculate fitted values for a function reegression model
#   based on B-spline basis functions
#
# Arguments
#
# w           bin IDs for observations
# x           predictor variables
# groups      group IDs for clustering variables
# beta        coefficients for main predictors
# gamma       random curves for group variables
# theta1      knots for beta B-splines
# theta2      knots for gamma B-splines
# n           number of subjects
# n_j         number of observations for each subject
# n_k         number of predictor variables
# n_q         number of grouping variables
# n_G_q       number of groups for each grouping variable
# n_p         number of knots for each beta B-spline
# n_t         number of knots for each random gamma B-spline

fitted_calc <- function(w, x, groups, beta, gamma, theta1, theta2, n, n_j, n_k, n_q,
                        n_G_q, n_p, n_t, degree) {
  out <- matrix(0, nrow=n, ncol=max(n_j, na.rm=TRUE))
  if (length(w) == n_j) {
    w <- matrix(rep(w, n), ncol=n_j, byrow=TRUE)
  }
  for (k in 1:n_k) {
    b_sp <- calc_b_splines(theta1[[k]], w, degree=degree)
    out <- out + sweep(apply(b_sp, 3, function(x) x %*% beta[[k]]), 1, x[, k], "*")
  }
  for (q in 1:n_q) {
    for (qq in 1:n_G_q[q]) {
      b_sp <- calc_b_splines(theta2[[q]][[qq]], w, degree=degree)
      out <- out + ifelse(groups[, q] == c(unique(groups[, q])[qq]), 1, 0) *
             apply(b_sp, 3, function(x) x %*% gamma[[q]][[qq]])
    }
  }
  return(out)
}