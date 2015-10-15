/*
 *  Gibbs update to interface with R
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<R.h>
#include<Rmath.h>

void update_beta(double *response, double *preds, int *n, int *np, int *nj, double *betas, double *vari, double *psis, double *phis, double *var_beta)
{
/* Gibbs update for all parameters. Updates beta and then sigma2 */
  int i, j, k, m;
  double cov_sum, cov_sum2, cov_sum4;
  double beta_mean;
  double beta_sigma2;
  double phi_new, psi_new;

  int num_community = *n, num_bins = *nj, num_covariates = *np;
  double psi = *psis, phi = *phis, sigma2_beta = *var_beta;
  double sigma2;

  double weight_hist[num_community][num_bins];
  double covariates[num_community][num_covariates];
  double beta[num_covariates][num_bins];

  //Added this line - not sure if it will fix it or not - sigma2 not passing back to R.
  //Consider changing the gamma function - this might fix it.
  sigma2 = *vari;
  for (i = 0; i < num_community; i++) {
    for (j = 0; j < num_bins; j++) {
      weight_hist[i][j] = response[i * num_bins + j];
    }
    for (k = 0; k < num_covariates; k++) {
      covariates[i][k] = preds[i * num_covariates + k];
    }
  }
  for (k = 0; k < num_covariates; k++) {
    for (j = 0; j < num_bins; j++) {
      beta[k][j] = betas[k * num_bins + j];
    }
  }

  GetRNGstate();

  //Update beta
  // j = 0
  for (k = 0; k < num_covariates ; k++) {
    cov_sum = 0;
    cov_sum2 = 0;
    for (i = 0; i < num_community; i++) {
      cov_sum4 = 0;
      cov_sum2 += covariates[i][k] * covariates[i][k];
      for (m = 0; m < num_covariates; m++) {
        if (m != k) {
          cov_sum4 += beta[m][0] * covariates[i][k] * covariates[i][m];
        }
      }
      cov_sum += 2.0 * weight_hist[i][0] * covariates[i][k] - 2.0 * cov_sum4;
    }
    beta_mean = (sigma2 * (6.0 * beta[k][1] - 2.0 * beta[k][2])
                       + sigma2_beta * cov_sum) / (6.0 * sigma2 + 2.0 * sigma2_beta * cov_sum2);
    beta_sigma2 = (sigma2 * sigma2_beta) / (3.0 * sigma2 + sigma2_beta * cov_sum2);
    //Sample new beta[k][0]
    beta[k][0] = rnorm(beta_mean, sqrt(beta_sigma2));
  }
  // j = 1
  for (k = 0; k < num_covariates ; k++) {
    cov_sum = 0;
    cov_sum2 = 0;
    for (i = 0; i < num_community; i++) {
      cov_sum4 = 0;
      cov_sum2 += covariates[i][k] * covariates[i][k];
      for (m = 0; m < num_covariates; m++) {
        if (m != k) {
          cov_sum4 += beta[m][1] * covariates[i][k] * covariates[i][m];
        }
      }
      cov_sum += 2.0 * weight_hist[i][1] * covariates[i][k] - 2.0 * cov_sum4;
    }
    beta_mean = (sigma2 * (6.0 * beta[k][0] + 8.0 * beta[k][2] - 2.0 * beta[k][3])
                       + sigma2_beta * cov_sum) / (12.0 * sigma2 + 2.0 * sigma2_beta * cov_sum2);
    beta_sigma2 = (sigma2 * sigma2_beta) / (6.0 * sigma2 + sigma2_beta * cov_sum2);
    //Sample new beta[k][1]
    beta[k][1] = rnorm(beta_mean, sqrt(beta_sigma2));
  }
  // 1 < j < max(j) - 1
  for (j = 2; j < (num_bins - 2); j++) {
    for (k = 0; k < num_covariates ; k++) {
      cov_sum = 0;
      cov_sum2 = 0;
      for (i = 0; i < num_community; i++) {
        cov_sum4 = 0;
        cov_sum2 += covariates[i][k] * covariates[i][k];
        for (m = 0; m < num_covariates; m++) {
          if (m != k) {
            cov_sum4 += beta[m][j] * covariates[i][k] * covariates[i][m];
          }
        }
        cov_sum += 2.0 * weight_hist[i][j] * covariates[i][k] - 2.0 * cov_sum4;
      }
      beta_mean = (sigma2 * (-2.0 * beta[k][j-2] + 8.0 * beta[k][j-1] + 8.0 * beta[k][j+1] - 2.0 * beta[k][j+2])
                         + sigma2_beta * cov_sum) / (12.0 * sigma2 + 2.0 * sigma2_beta * cov_sum2);
      beta_sigma2 = (sigma2 * sigma2_beta) / (6.0 * sigma2 + sigma2_beta * cov_sum2);
      //Sample new beta[k][j]
      beta[k][j] = rnorm(beta_mean, sqrt(beta_sigma2));
    }
  }
  // j = max(j) - 1
  for (k = 0; k < num_covariates ; k++) {
    cov_sum = 0;
    cov_sum2 = 0;
    for (i = 0; i < num_community; i++) {
      cov_sum4 = 0;
      cov_sum2 += covariates[i][k] * covariates[i][k];
      for (m = 0; m < num_covariates; m++) {
        if (m != k) {
          cov_sum4 += beta[m][num_bins - 2] * covariates[i][k] * covariates[i][m];
        }
      }
      cov_sum += 2.0 * weight_hist[i][num_bins - 2] * covariates[i][k] - 2.0 * cov_sum4;
    }
    beta_mean = (sigma2 * (-2.0 * beta[k][num_bins - 4] + 8.0 * beta[k][num_bins - 3] + 4.0 * beta[k][num_bins - 1])
                       + sigma2_beta * cov_sum) / (10.0 * sigma2 + 2.0 * sigma2_beta * cov_sum2);
    beta_sigma2 = (sigma2 * sigma2_beta) / (5.0 * sigma2 + sigma2_beta * cov_sum2);
    //Sample new beta[k][num_bins - 2]
    beta[k][num_bins - 2] = rnorm(beta_mean, sqrt(beta_sigma2));
  }
  // j = max(j)
  for (k = 0; k < num_covariates ; k++) {
    cov_sum = 0;
    cov_sum2 = 0;
    for (i = 0; i < num_community; i++) {
      cov_sum4 = 0;
      cov_sum2 += covariates[i][k] * covariates[i][k];
      for (m = 0; m < num_covariates; m++) {
        if (m != k) {
          cov_sum4 += beta[m][num_bins - 1] * covariates[i][k] * covariates[i][m];
        }
      }
      cov_sum += 2.0 * weight_hist[i][num_bins - 1] * covariates[i][k] - 2.0 * cov_sum4;
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
    for (j = 0; j < num_bins; j++) {
      cov_sum2 = 0;
      for (k = 0; k < num_covariates; k++) {
        cov_sum2 += beta[k][j] * covariates[i][k];
      }
      cov_sum += (weight_hist[i][j] - cov_sum2) * (weight_hist[i][j] - cov_sum2);
    }
  }
  psi_new = psi + num_covariates * num_bins / 2.0;
  phi_new = phi + (1.0 / 2.0) * cov_sum;
  //Sample new sigma2
  sigma2 = 1 / rgamma(psi_new, 1 / phi_new);

  //Update and return values
  for (i = 0; i < num_community; i++) {
    for (j = 0; j < num_bins; j++) {
      response[i * num_bins + j] = weight_hist[i][j];
    }
    for (k = 0; k < num_covariates; k++) {
      preds[i * num_covariates + k] = covariates[i][k];
    }
  }
  for (k = 0; k < num_covariates; k++) {
    for (j = 0; j < num_bins; j++) {
      betas[k * num_bins + j] = beta[k][j];
    }
  }
  *vari = sigma2;

  PutRNGstate();
}