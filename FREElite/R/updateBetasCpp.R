updateBetas <- function(y, x, groups, beta, gamma, rho, sigma2,
sigma2_gamma, b_splines_mat, sigma2_hyper_a, sigma2_hyper_b,
sigma2_gamma_hyper_a, sigma2_gamma_hyper_b, beta_hyper,
n, n_j, n_k, n_q, n_G_q, n_p, n_t, bin_id){
  
  for (i in 1:n_k) {
    for (j in 1:n_p) {
      beta[i, j] <- slice_sample_beta(beta[i, j], c(i, j), y, x, groups,
      beta, gamma, sigma2, rho, b_splines_mat,
      beta_hyper, n, n_j, n_k, n_q,
      sigma2_gamma, bin_id)
    }
  }
  
  for (q in 1:n_q) {
    gamma_tmp <- rep(0, n_t)
    for (j in 1:(n_G_q[q] - 1)) {
      for (i in 1:n_t) {
        gamma[[q]][j, i] <- slice_sample_gamma(gamma[[q]][j, i], c(q, i, j),
        y, x, groups, beta, gamma,
        sigma2, rho, b_splines_mat, beta_hyper,
        n, n_j, n_k, n_q, sigma2_gamma, bin_id)
      }
      gamma_tmp <- gamma_tmp + c((ginv(b_splines_mat) %*% b_splines_mat) %*% gamma[[q]][j, ])
      if (n_G_q[q] > 1) {
        gamma[[q]][n_G_q[q], ] <- -c(gamma_tmp)
      } else {
        gamma[[q]][n_G_q[q], ] <- rep(0, n_t)
      }
    }
  }
  
  sigma2 <- sample_sigma2(y, x, groups, beta, gamma, sigma2, rho, b_splines_mat,
  sigma2_hyper_a, sigma2_hyper_b,
  n, n_j, n_k, n_q, bin_id)
  
  rho <- slice_sample_rho(rho, c(1, 2), y, x, groups, beta, gamma,
  sigma2, rho, b_splines_mat, beta_hyper, n, n_j,
  n_k, n_q, sigma2_gamma, bin_id)
  
  for (q in 1:n_q) {
    sigma2_gamma[q] <- sample_sigma2_gamma(gamma, sigma2_gamma_hyper_a,
    sigma2_gamma_hyper_b, n_q, q)
  }
  
  return(list(beta=beta, gamma=gamma, rho=rho, sigma2=sigma2,
  sigma2_gamma=sigma2_gamma))
}