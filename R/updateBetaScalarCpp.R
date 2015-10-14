update_betas_scalar <- function(y, x, z, groups, alpha, beta, gamma, delta, sigma2, sigma2_gamma,
                                bs_beta, phi1, psi1, phi2, psi2, s2_alpha, s2_beta, s2_delta,
                                n_q, n_G_q) {
  # update alpha
  alpha <- sample_alpha_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
                        z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta, s2_alpha=s2_alpha)

  # update betas
  for (p in 1:length(beta)) {
    beta[p] <- sample_beta_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
                           z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta,
                           s2_beta=s2_beta, p_id=(p - 1))
  }
  
  # update deltas
  for (k in 1:length(delta)) {
    delta[k] <- sample_delta_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma,
                             delta=delta, z=z, alpha=alpha, sigma2=sigma2,
                             bs_beta=bs_beta, s2_delta=s2_delta, k_id=(k - 1))
  }

  # update gammas
  for (q in 1:n_q) {
    for (g in 1:(n_G_q[q] - 1)) {
      gamma[[q]][g] <- sample_gamma_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma,
                                    delta=delta, z=z, alpha=alpha, sigma2=sigma2,
                                    bs_beta=bs_beta, sigma2_gamma=sigma2_gamma,
                                    q_id=(q - 1), G_id=(g - 1))  
    }
    gamma[[q]][n_G_q[q]] <- -sum(gamma[[q]][1:(n_G_q[q] - 1)])
  }
            
  # update sigma2
  sigma2 <- sample_sigma2_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
                          z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta,
                          phi1=phi1, psi1=psi1)

  # update sigma2_gamma
  for (q in 1:n_q) {
  	sigma2_gamma[q] <- sample_sigma2_gamma_scalar(gamma=gamma, phi2=phi2, psi2=psi2, q_id=(q - 1))
  }
                      
  # calculate lnL for checks
  lnL <- lnL_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
             z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta)

  return(list(alpha=alpha, beta=beta, gamma=gamma, delta=delta, sigma2=sigma2,
              sigma2_gamma=sigma2_gamma, lnL=lnL))
}

update_betas_scalar2 <- function(y, x, z, groups, alpha, beta, gamma, delta, sigma2, sigma2_gamma,
                                 bs_beta, phi1, psi1, phi2, psi2, s2_alpha, s2_beta, s2_delta,
                                 n_q, n_G_q) {
  # update alpha
  alpha <- sample_alpha_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
                        z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta, s2_alpha=s2_alpha)
  
  # update betas
  for (p in 1:length(beta)) {
    beta[p] <- sample_beta_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
                           z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta,
                           s2_beta=s2_beta, p_id=(p - 1))
  }
  
  # update deltas
  for (k in 1:length(delta)) {
    delta[k] <- sample_delta_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma,
                             delta=delta, z=z, alpha=alpha, sigma2=sigma2,
                             bs_beta=bs_beta, s2_delta=s2_delta, k_id=(k - 1))
  }
  
  # update gammas
  for (q in 1:n_q) {
    for (g in 1:n_G_q[q]) {
      gamma[[q]][g] <- 0
    }
  }
  
  # update sigma2
  sigma2 <- sample_sigma2_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
                          z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta,
                          phi1=phi1, psi1=psi1)
  
  # update sigma2_gamma
  for (q in 1:n_q) {
    sigma2_gamma[q] <- sample_sigma2_gamma_scalar(gamma=gamma, phi2=phi2, psi2=psi2, q_id=(q - 1))
  }
  
  # calculate lnL for checks
  lnL <- lnL_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
             z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta)
  
  return(list(alpha=alpha, beta=beta, gamma=gamma, delta=delta, sigma2=sigma2,
              sigma2_gamma=sigma2_gamma, lnL=lnL))
}

update_betas_scalar3 <- function(y, x, z, groups, alpha, beta, gamma, delta, sigma2, sigma2_gamma,
                                 bs_beta, phi1, psi1, phi2, psi2, s2_alpha, s2_beta, s2_delta,
                                 n_q, n_G_q) {
  # update alpha
  alpha <- sample_alpha_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
                        z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta, s2_alpha=s2_alpha)
  
  # update betas
  for (p in 1:length(beta)) {
    beta[p] <- sample_beta_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
                           z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta,
                           s2_beta=s2_beta, p_id=(p - 1))
  }
  
  # update deltas
  for (k in 1:length(delta)) {
    delta[k] <- 0
  }
  
  # update gammas
  for (q in 1:n_q) {
    for (g in 1:(n_G_q[q] - 1)) {
      gamma[[q]][g] <- sample_gamma_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma,
                                    delta=delta, z=z, alpha=alpha, sigma2=sigma2,
                                    bs_beta=bs_beta, sigma2_gamma=sigma2_gamma,
                                    q_id=(q - 1), G_id=(g - 1))  
    }
    gamma[[q]][n_G_q[q]] <- -sum(gamma[[q]][1:(n_G_q[q] - 1)])
  }
  
  # update sigma2
  sigma2 <- sample_sigma2_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
                          z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta,
                          phi1=phi1, psi1=psi1)
  
  # update sigma2_gamma
  for (q in 1:n_q) {
    sigma2_gamma[q] <- sample_sigma2_gamma_scalar(gamma=gamma, phi2=phi2, psi2=psi2, q_id=(q - 1))
  }
  
  # calculate lnL for checks
  lnL <- lnL_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
             z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta)
  
  return(list(alpha=alpha, beta=beta, gamma=gamma, delta=delta, sigma2=sigma2,
              sigma2_gamma=sigma2_gamma, lnL=lnL))
}

update_betas_scalar4 <- function(y, x, z, groups, alpha, beta, gamma, delta, sigma2, sigma2_gamma,
                                 bs_beta, phi1, psi1, phi2, psi2, s2_alpha, s2_beta, s2_delta,
                                 n_q, n_G_q) {
  # update alpha
  alpha <- sample_alpha_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
                        z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta, s2_alpha=s2_alpha)
  
  # update betas
  for (p in 1:length(beta)) {
    beta[p] <- sample_beta_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
                           z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta,
                           s2_beta=s2_beta, p_id=(p - 1))
  }
  
  # update deltas
  for (k in 1:length(delta)) {
    delta[k] <- 0
  }
  
  # update gammas
  for (q in 1:n_q) {
    for (g in 1:n_G_q[q]) {
      gamma[[q]][g] <- 0
    }
  }
  
  # update sigma2
  sigma2 <- sample_sigma2_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
                          z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta,
                          phi1=phi1, psi1=psi1)
  
  # update sigma2_gamma
  for (q in 1:n_q) {
    sigma2_gamma[q] <- sample_sigma2_gamma_scalar(gamma=gamma, phi2=phi2, psi2=psi2, q_id=(q - 1))
  }
  
  # calculate lnL for checks
  lnL <- lnL_scalar(y=y, x=x, groups=groups, beta=beta, gamma=gamma, delta=delta,
             z=z, alpha=alpha, sigma2=sigma2, bs_beta=bs_beta)
  
  return(list(alpha=alpha, beta=beta, gamma=gamma, delta=delta, sigma2=sigma2,
              sigma2_gamma=sigma2_gamma, lnL=lnL))
}