/*
 *  spline_iid.c
 *  
 *
 *  Created by Jian Yen on 30/04/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<R.h>
#include<Rmath.h>

void update_beta_iid(double *response, double *preds, int *n, int *np, int *nj, double *betas, double *vari,
                 double *psis, double *phis, double *var_beta)
{
/* Gibbs update for all parameters. Updates beta and then sigma2 */
  int i, j, k, m;
  double cov_sum, cov_sum2, cov_sum4;
  double beta_mean;	
  double beta_sigma2;
  double phi_new, psi_new;
  
  //Allocate inputs to objects
  int num_community = *n;
  int num_covariates = *np;
  int num_bins = *nj;
  double phi = *phis, psi = *psis;
  double sigma2_beta = *var_beta;
  double sigma2 = *vari;
  
  //Change data vectors to arrays
  double weight_hist[num_community][num_bins];
  double covariates[num_community][num_covariates];
  double beta[num_covariates][num_bins];
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
  for (j = 0; j < num_bins; j++) {
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
        //Should be 2.0 * cov_sum4 but this causes the value of cov_sum to grow exponentially
        cov_sum += weight_hist[i][j] * covariates[i][k] - cov_sum4;
      }
      beta_mean = (sigma2_beta * cov_sum) / (sigma2 + sigma2_beta * cov_sum2);
      beta_sigma2 = (sigma2 * sigma2_beta) / (sigma2 + sigma2_beta * cov_sum2);
      //Sample new beta[k][j]
      beta[k][j] = rnorm(beta_mean, sqrt(beta_sigma2));
    }
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