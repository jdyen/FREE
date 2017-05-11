/*
 *  Gibbs update to interface with R
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<R.h>
#include<Rmath.h>

void update_beta_scalar(double *response, double *preds, int *n, int *np, int *nj, double *alpha, double *betas, double *vari, double *psis, double *phis, double *var_beta, double *var_alpha)
{
/* Gibbs update for all parameters. Updates beta and then sigma2 */
  int i, j, k, l, m;
  double cov_sum, cov_sum2, cov_sum4;
  double alpha_mean;
  double alpha_sigma2;
  double beta_mean;
  double beta_sigma2;
  double phi_new, psi_new;

  int num_community = *n, num_bins = *nj, num_covariates = *np;
  double psi = *psis, phi = *phis, sigma2_beta = *var_beta, sigma2_alpha = *var_alpha;
  double alpha_est;
  double sigma2;

  double resp[num_community];
  double covariates[num_community][num_bins][num_covariates];
  double beta[num_covariates][num_bins];

  sigma2 = *vari;
  alpha_est = *alpha;
  for (i = 0; i < num_community; i++) {
    resp[i] = response[i];
    for (j = 0; j < num_bins; j++) {
      for (k = 0; k < num_covariates; k++) {
        covariates[i][j][k] = preds[i * num_bins * num_covariates + k * num_bins + j];
      }
    }
  }
  for (k = 0; k < num_covariates; k++) {
    for (j = 0; j < num_bins; j++) {
      beta[k][j] = betas[k * num_bins + j];
    }
  }

  GetRNGstate();

  //Update alpha_est
  cov_sum = 0;
  for (i = 0; i < num_community; i++) {
    cov_sum4 = 0;
    for (k = 0; k < num_covariates; k++) {
      for (j = 0; j < num_bins; j++) {
        cov_sum4 += beta[k][j] * covariates[i][j][k];
      }
    }
    cov_sum += 2.0 * resp[i] - 2.0 * cov_sum4;
  }
  alpha_mean = cov_sum / (2.0 * sigma2);
  alpha_sigma2 = (sigma2 * sigma2_alpha) / (sigma2 + sigma2_alpha);
  //Sample new alpha_est
  alpha_est = rnorm(alpha_mean, sqrt(alpha_sigma2));
  
  //Update beta
  // j = 1 (index 0 here)
  for (k = 0; k < num_covariates ; k++) {
    cov_sum = 0;
    cov_sum2 = 0;
    for (i = 0; i < num_community; i++) {
      cov_sum4 = 0;
      cov_sum2 += covariates[i][0][k] * covariates[i][0][k];
      for (m = 0; m < num_covariates; m++) {
        for (l = 0; l < num_bins; l++) {
          if ((m != k) | (l != 0)) {
            cov_sum4 += beta[m][l] * covariates[i][l][m];
          }
        }
      }
      cov_sum += 2.0 * resp[i] * covariates[i][0][k] + 2.0 * alpha_est * covariates[i][0][k]
                 - 2.0 * covariates[i][0][k] * cov_sum4;
    }
    beta_mean = (sigma2 * (6.0 * beta[k][1] - 2.0 * beta[k][2])
                 + sigma2_beta * cov_sum) / (6.0 * sigma2 + 2.0 * sigma2_beta * cov_sum2);
    beta_sigma2 = (sigma2 * sigma2_beta) / (3.0 * sigma2 + sigma2_beta * cov_sum2);
    //Sample new beta[k][0]
    beta[k][0] = rnorm(beta_mean, sqrt(beta_sigma2));
  }
  // j = 2 (index 1 here)
  for (k = 0; k < num_covariates ; k++) {
    cov_sum = 0;
    cov_sum2 = 0;
    for (i = 0; i < num_community; i++) {
      cov_sum4 = 0;
      cov_sum2 += covariates[i][1][k] * covariates[i][1][k];
      for (m = 0; m < num_covariates; m++) {
        for (l = 0; l < num_bins; l++) {
          if ((m != k) | (l != 1)) {
            cov_sum4 += beta[m][l] * covariates[i][l][m];
          }
        }
      }
      cov_sum += 2.0 * resp[i] * covariates[i][1][k] + 2.0 * alpha_est * covariates[i][1][k]
                 - 2.0 * covariates[i][1][k] * cov_sum4;
    }
    beta_mean = (sigma2 * (6.0 * beta[k][0] + 8.0 * beta[k][2] - 2.0 * beta[k][3])
                 + sigma2_beta * cov_sum) / (12.0 * sigma2 + 2.0 * sigma2_beta * cov_sum2);
    beta_sigma2 = (sigma2 * sigma2_beta) / (6.0 * sigma2 + sigma2_beta * cov_sum2);
    //Sample new beta[k][1]
    beta[k][1] = rnorm(beta_mean, sqrt(beta_sigma2));
  }
  // 2 < j < max(j) - 1
  for (j = 2; j < (num_bins - 2); j++) {
    for (k = 0; k < num_covariates ; k++) {
      cov_sum = 0;
      cov_sum2 = 0;
      for (i = 0; i < num_community; i++) {
        cov_sum4 = 0;
        cov_sum2 += covariates[i][j][k] * covariates[i][j][k];
        for (m = 0; m < num_covariates; m++) {
          for (l = 0; l < num_bins; l++) {
            if ((m != k) | (l != j)) {
              cov_sum4 += beta[m][l] * covariates[i][l][m];
            }
          }
        }
        cov_sum += 2.0 * resp[i] * covariates[i][j][k] + 2.0 * alpha_est * covariates[i][j][k]
                   - 2.0 * covariates[i][j][k] * cov_sum4;
      }
      beta_mean = (sigma2 * (-2.0 * beta[k][j-2] + 8.0 * beta[k][j-1] + 8.0 * beta[k][j+1]
                             - 2.0 * beta[k][j+2])
                   + sigma2_beta * cov_sum) / (12.0 * sigma2 + 2.0 * sigma2_beta * cov_sum2);
      beta_sigma2 = (sigma2 * sigma2_beta) / (6.0 * sigma2 + sigma2_beta * cov_sum2);
      //Sample new beta[k][j]
      beta[k][j] = rnorm(beta_mean, sqrt(beta_sigma2));
    }
  }
  // j = max(j) - 1 (index j - 2 here)
  for (k = 0; k < num_covariates ; k++) {
    cov_sum = 0;
    cov_sum2 = 0;
    for (i = 0; i < num_community; i++) {
      cov_sum4 = 0;
      cov_sum2 += covariates[i][num_bins - 2][k] * covariates[i][num_bins - 2][k];
      for (m = 0; m < num_covariates; m++) {
        for (l = 0; l < num_bins; l++) {
          if ((m != k) | (l != (num_bins - 2))) {
            cov_sum4 += beta[m][l] * covariates[i][l][m];
          }
        }
      }
      cov_sum += 2.0 * resp[i] * covariates[i][num_bins - 2][k]
                 + 2.0 * alpha_est * covariates[i][num_bins - 2][k]
                 - 2.0 * covariates[i][num_bins - 2][k] * cov_sum4;
    }
    beta_mean = (sigma2 * (-2.0 * beta[k][num_bins - 4] + 8.0 * beta[k][num_bins - 3]
                           + 4.0 * beta[k][num_bins - 1])
                       + sigma2_beta * cov_sum) / (10.0 * sigma2 + 2.0 * sigma2_beta * cov_sum2);
    beta_sigma2 = (sigma2 * sigma2_beta) / (5.0 * sigma2 + sigma2_beta * cov_sum2);
    //Sample new beta[k][num_bins - 2]
    beta[k][num_bins - 2] = rnorm(beta_mean, sqrt(beta_sigma2));
  }
  // j = max(j) (index j - 1 here)
  for (k = 0; k < num_covariates ; k++) {
    cov_sum = 0;
    cov_sum2 = 0;
    for (i = 0; i < num_community; i++) {
      cov_sum4 = 0;
      cov_sum2 += covariates[i][num_bins - 1][k] * covariates[i][num_bins - 1][k];
      for (m = 0; m < num_covariates; m++) {
        for (l = 0; l < num_bins; l++) {
          if ((m != k) | (l != (num_bins - 1))) {
            cov_sum4 += beta[m][l] * covariates[i][l][m];
          }
        }
      }
      cov_sum += 2.0 * resp[i] * covariates[i][num_bins - 1][k]
                 + 2.0 * alpha_est * covariates[i][num_bins - 1][k]
                 - 2.0 * covariates[i][num_bins - 1][k] * cov_sum4;
    }
    beta_mean = (sigma2 * (-2.0 * beta[k][num_bins - 3] + 4.0 * beta[k][num_bins - 2])
                 + sigma2_beta * cov_sum) / (2.0 * sigma2 + 2.0 * sigma2_beta * cov_sum2);
    beta_sigma2 = (sigma2 * sigma2_beta) / (sigma2 + sigma2_beta * cov_sum2);
    //Sample new beta[k][num_bins - 1]
    beta[k][num_bins - 1] = rnorm(beta_mean, sqrt(beta_sigma2));
  }
    
  //Update sigma2
  cov_sum = 0;
  for (i = 0; i < num_community; i++) {
    cov_sum2 = 0;
    for (j = 0; j < num_bins; j++) {
      for (k = 0; k < num_covariates; k++) {
        cov_sum2 += beta[k][j] * covariates[i][j][k];
      }
    }
    cov_sum += (resp[i] - alpha_est - cov_sum2) * (resp[i] - alpha_est - cov_sum2);
  }
  psi_new = psi + num_community / 2.0;
  phi_new = phi + (1.0 / 2.0) * cov_sum;
  //Sample new sigma2
  sigma2 = 1 / rgamma(psi_new, 1 / phi_new);

  //Update and return values
  for (i = 0; i < num_community; i++) {
    response[i] = resp[i];
    for (j = 0; j < num_bins; j++) {
      for (k = 0; k < num_covariates; k++) {
        preds[i * num_covariates * num_bins + j * num_covariates + k] = covariates[i][j][k];
      }
    }
  }
  for (k = 0; k < num_covariates; k++) {
    for (j = 0; j < num_bins; j++) {
      betas[k * num_bins + j] = beta[k][j];
    }
  }
  *vari = sigma2;
  *alpha = alpha_est;

  PutRNGstate();
}