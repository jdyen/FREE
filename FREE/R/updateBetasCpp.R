update_betas <- function(y, x, groups, w, beta, gamma, rho, theta1, theta2, sigma2,
                         sigma2_gamma, b_splines_beta, b_splines_gamma,
                         sigma2_hyper_a, sigma2_hyper_b,
                         sigma2_gamma_hyper_a, sigma2_gamma_hyper_b, beta_hyper,
                         n, n_j, n_k, n_q, n_G_q, n_p, n_t, degree)
{

for (i in 1:n_k) {
  for (j in 1:n_p[i]) {
    beta[[i]][j] <- slice_sample_beta(beta[[i]][j], c(i, j), 
                                      y, x, w, groups, beta, gamma, theta1,
                                      theta2, sigma2, rho, b_splines_beta,
                                      b_splines_gamma, sigma2_hyper_a, sigma2_hyper_b,
                                      sigma2_gamma_hyper_a, sigma2_gamma_hyper_b,
                                      beta_hyper, n, n_j, n_k, n_q, sigma2_gamma)
  }
}

for (q in 1:n_q) {
  gamma_tmp <- rep(0, n_t[[q]][n_G_q[q]])
  for (j in 1:(n_G_q[q] - 1)) {
    for (i in 1:n_t[[q]][j]) {
      gamma[[q]][[j]][i] <- slice_sample_gamma(gamma[[q]][[j]][i],
                                               c(q, i, j), 
                                               y, x, w, groups, beta, gamma, theta1,
                                               theta2, sigma2, rho, b_splines_beta,
                                               b_splines_gamma, sigma2_hyper_a,
                                               sigma2_hyper_b, sigma2_gamma_hyper_a,
                                               sigma2_gamma_hyper_b, beta_hyper,
                                               n, n_j, n_k, n_q, sigma2_gamma)
    }
    b_sp <- calc_bs(w, theta2[[q]][[j]], degree, c(min(w) - 1, max(w) + 1))
    b_sp2 <- calc_bs(w, theta2[[q]][[n_G_q[q]]], degree, c(min(w) - 1, max(w) + 1))
    gamma_tmp <- gamma_tmp + c((ginv(b_sp2) %*% b_sp) %*% gamma[[q]][[j]])
    if (n_G_q[q] > 1) {
      gamma[[q]][[n_G_q[q]]] <- -c(gamma_tmp)
    } else {
      gamma[[q]][[n_G_q[q]]] <- rep(0, n_t[[q]][n_G_q[q]])
    }
  }
}

sigma2 <- sample_sigma2(y, x, w, groups, beta, gamma, theta1, theta2,
                        sigma2, rho, b_splines_beta, b_splines_gamma,
                        sigma2_hyper_a, sigma2_hyper_b,
                        sigma2_gamma_hyper_a, sigma2_gamma_hyper_b,
                        beta_hyper, n, n_j, n_k, n_q)

rho <- slice_sample_rho(rho, c(1, 2), y, x, w, groups, beta,
                        gamma, theta1,
                        theta2, sigma2, rho, b_splines_beta, b_splines_gamma,
                        sigma2_hyper_a, sigma2_hyper_b, sigma2_gamma_hyper_a,
                        sigma2_gamma_hyper_b, beta_hyper, n, n_j, n_k, n_q,
                        sigma2_gamma)

# calculate B-spline value
b_splines_beta <- vector("list", length=n_k)
for (i in 1:n_k) {
  b_splines_beta[[i]] <- calc_bs(w, theta1[[i]], degree,
                                 c(min(w, na.rm=TRUE) - 1, max(w, na.rm=TRUE) + 1))
}
b_splines_gamma <- vector("list", length=n_q)
for (i in 1:n_q) {
  b_splines_gamma[[i]] <- vector("list", length=n_G_q[i])
  for (q in 1:n_G_q[i]) {
    b_splines_gamma[[i]][[q]] <- calc_bs(w, theta2[[i]][[q]], degree,
                                         c(min(w, na.rm=TRUE) - 1, max(w, na.rm=TRUE) + 1))
  }
}

for (q in 1:n_q) {
  sigma2_gamma[q] <- sample_sigma2_gamma(gamma, sigma2_gamma_hyper_a,
                                         sigma2_gamma_hyper_b,
                                         n_q, q)
}

return(list(beta=beta, gamma=gamma, rho=rho, theta1=theta1, theta2=theta2,
            sigma2=sigma2, sigma2_gamma=sigma2_gamma, b_splines_beta=b_splines_beta,
            b_splines_gamma=b_splines_gamma))
}