#include <Rcpp.h>

using namespace Rcpp;

// function to update parameters for spline function regression model (scalar response)
//
// Jian Yen
// This is a function for Gibbs sampler updates for all the parameters
//  in a B-spline function regression model with a scalar response and
//  function-valued predictor. Model will be called and used in R.
//
//  Created by Jian Yen on 12/01/2015.
//
//

// [[Rcpp::export]]
double lnL_scalar(NumericVector y, List x, NumericMatrix groups,
           NumericMatrix beta, List gamma, NumericVector delta,
           NumericMatrix z, double alpha, double sigma2,
           List bs_beta) {
  double out;
  
  //Calculate mean value for each observation
  double sum_sq = 0;
  for (int i = 0; i < y.size(); i++) {
    double mu = alpha;
    for (int k = 0; k < beta.nrow(); k++) {
      NumericMatrix bs_beta_k = bs_beta[k];
      NumericMatrix x_k = x[k];
      for (int j = 0; j < x_k.ncol(); j++) {
        double mu_tmp = 0;
        for (int jj = 0; jj < beta.ncol(); jj++) {
          mu_tmp += bs_beta_k(j, jj) * beta(k, jj);
        }
        mu += x_k(i, j) * mu_tmp;
      }
    }
    for (int k = 1; k < z.ncol(); k++) {
      mu += delta[k] * z(i, k);
    }
    for (int q = 0; q < gamma.size(); q++) {
      NumericVector gamma_q = gamma[q];
      int group_set = groups(i, q);
      mu += gamma_q(group_set);
    }
    sum_sq += (y(i) - mu) * (y(i) - mu);
  }

  //Calculate loglik
  out = (- (double) y.size() / 2.0) * log(2.0 * PI * sigma2) - (1.0 / (2.0 * sigma2)) * sum_sq;
  return out;
}

// [[Rcpp::export]]
double sample_alpha_scalar(NumericVector y, List x, NumericMatrix groups,
                    NumericMatrix beta, List gamma, NumericVector delta,
                    NumericMatrix z, double alpha, double sigma2,
                    List bs_beta, double s2_alpha) {
  double out;
  NumericVector params(2);
  
  double sum_term = 0.0;
  double sum_term2 = 0.0;
  double sum_term3 = 0.0;
  double sum_term4 = 0.0;
  for (int i = 0; i < y.size(); i++) {
    double term = 0.0;
    for (int k = 0; k < beta.nrow(); k++) {
      NumericMatrix bs_beta_k = bs_beta[k];
      NumericMatrix x_k = x[k];
      for (int j = 0; j < x_k.ncol(); j++) {
        for (int p = 0; p < beta.ncol(); p++) {
          term += beta(k, p) * bs_beta_k(j, p) * x_k(i, j);
        }
      }
    }
    sum_term += term;
    sum_term2 += y[i];
    for (int k = 0; k < z.ncol(); k++) {
      sum_term4 += delta[k] * z(i, k);
    }
    for (int q = 0; q < gamma.size(); q++) {
      NumericVector gamma_q = gamma[q];
      int group_set = groups(i, q);
      sum_term3 += gamma_q[group_set];
    }
  }

  //Calculate params[1]
  params[1] = (sigma2 * s2_alpha) / (s2_alpha * y.size() + sigma2);

  //Calculate params[0]
  params[0] = params[1] * (-1.0 / (2.0 * sigma2)) * (-2.0 * sum_term2 + 2.0 * sum_term +
                                                     2.0 * sum_term4 + 2.0 * sum_term3);

  //Calculate loglik
  out = rnorm(1, params[0], sqrt(params[1]))[0];
  return out;
}

// [[Rcpp::export]]
double sample_beta_scalar(NumericVector y, List x, NumericMatrix groups,
                   NumericMatrix beta, List gamma, NumericVector delta,
                   NumericMatrix z, double alpha, double sigma2,
                   List bs_beta, double s2_beta, int p_id, int k_id) {
  double out;
  NumericVector params(2);
  
  double sum_term = 0.0;
  double sum_term2 = 0.0;
  double sum_term2a = 0.0;
  double sum_term2b = 0.0;
  double sum_term3 = 0.0;
  double sum_term4 = 0.0;
  NumericMatrix x_tmp = x(k_id);
  NumericMatrix bs_beta_tmp = bs_beta(k_id);
  for (int i = 0; i < y.size(); i++) {
    double term = 0.0;
    for (int j = 0; j < x_tmp.ncol(); j++) {
      term += bs_beta_tmp(j, p_id) * x_tmp(i, j);
    }
    sum_term += term * term;
    sum_term2 += y[i] * term;
    sum_term2a += alpha * term;
    for (int k = 0; k < z.ncol(); k++) {
      sum_term2b += term * delta[k] * z(i, k);
    }
    for (int q = 0; q < gamma.size(); q++) {
      NumericVector gamma_q = gamma[q];
      int group_set = groups(i, q);
      sum_term3 += gamma_q[group_set] * term;
    }
    for (int k2 = 0; k2 < beta.nrow(); k2++) {
      NumericMatrix bs_beta_k = bs_beta[k2];
      NumericMatrix x_k = x[k2];
      if (k2 != k_id) {
        for (int s = 0; s < beta.ncol(); s++) {
          for (int j = 0; j < x_k.ncol(); j++) {
            sum_term4 += term * beta(k2, s) * bs_beta_k(j, s) * x_k(i, j);
          }
        }
      } else {
        for (int s = 0; s < beta.ncol(); s++) {
          if (s != p_id) {
            for (int j = 0; j < x_k.ncol(); j++) {
              sum_term4 += term * beta(k2, s) * bs_beta_k(j, s) * x_k(i, j);
            }
          }
        }
      }
    }
  }
  
  //Calculate params[1]
  params[1] = (sigma2 * s2_beta) / (s2_beta * sum_term + sigma2);
  
  //Calculate params[0]
  params[0] = params[1] * (-1.0 / (2.0 * sigma2)) * (-2.0 * sum_term2 + 2.0 * sum_term2a +
                                                     2.0 * sum_term2b + 2.0 * sum_term4 +
                                                     2.0 * sum_term3);

  //Calculate loglik
  out = rnorm(1, params[0], sqrt(params[1]))[0];
  return out;
}

// [[Rcpp::export]]
double sample_delta_scalar(NumericVector y, List x, NumericMatrix groups,
                    NumericMatrix beta, List gamma, NumericVector delta,
                    NumericMatrix z, double alpha, double sigma2,
                    List bs_beta, double s2_delta, int k_id) {
  double out;
  NumericVector params(2);
  
  double sum_term = 0.0;
  double sum_term2 = 0.0;
  double sum_term2a = 0.0;
  double sum_term2b = 0.0;
  double sum_term3 = 0.0;
  double sum_term4 = 0.0;
  for (int i = 0; i < y.size(); i++) {
    double term = 0.0;
    for (int k = 0; k < beta.nrow(); k++) {
      NumericMatrix bs_beta_k = bs_beta[k];
      NumericMatrix x_k = x[k];
      for (int j = 0; j < x_k.ncol(); j++) {
        for (int p = 0; p < beta.ncol(); p++) {
          term += beta(k, p) * bs_beta_k(j, p) * x_k(i, j);
        }
      }
    }
    sum_term += z(i, k_id) * term;
    sum_term2 += y[i] * z(i, k_id);
    sum_term2a += alpha * z(i, k_id);
    sum_term2b += z(i, k_id) * z(i, k_id);
    for (int k = 0; k < z.ncol(); k++) {
      if (k != k_id) {
        sum_term4 += delta[k] * z(i, k) * z(i, k_id);
      }
    }
    for (int q = 0; q < gamma.size(); q++) {
      NumericVector gamma_q = gamma[q];
      int group_set = groups(i, q);
      sum_term3 += gamma_q[group_set] * z(i, k_id);
    }
  }
  
  //Calculate params[1]
  params[1] = (sigma2 * s2_delta) / (s2_delta * sum_term2b + sigma2);
  
  //Calculate params[0]
  params[0] = params[1] * (-1.0 / (2.0 * sigma2)) * (-2.0 * sum_term2 + 2.0 * sum_term +
                                                     2.0 * sum_term2a + 2.0 * sum_term4 +
                                                     2.0 * sum_term3);
  
  //Calculate loglik
  out = rnorm(1, params[0], sqrt(params[1]))[0];
  return out;
}

// [[Rcpp::export]]
double sample_gamma_scalar(NumericVector y, List x, NumericMatrix groups,
                    NumericMatrix beta, List gamma, NumericVector delta,
                    NumericMatrix z, double alpha, double sigma2,
                    List bs_beta, NumericVector sigma2_gamma,
                    int q_id, int G_id) {
  double out;
  NumericVector params(2);
  
  double sum_term = 0.0;
  double sum_term2 = 0.0;
  double sum_term2a = 0.0;
  double sum_term2b = 0.0;
  double sum_term3 = 0.0;
  double sum_term4 = 0.0;
  for (int i = 0; i < y.size(); i++) {
    if (groups(i, q_id) == G_id) {
      double term = 0.0;
      for (int k = 0; k < beta.nrow(); k++) {
        NumericMatrix bs_beta_k = bs_beta[k];
        NumericMatrix x_k = x[k];
        for (int j = 0; j < x_k.ncol(); j++) {
          for (int p = 0; p < beta.ncol(); p++) {
            term += beta(k, p) * bs_beta_k(j, p) * x_k(i, j);
          }
        }
      }
      sum_term += 1;
      sum_term2 += y[i];
      sum_term2a += alpha;
      sum_term2b += term;
      for (int k = 0; k < z.ncol(); k++) {
        sum_term4 += delta[k] * z(i, k);
      }
      for (int q = 0; q < gamma.size(); q++) {
        if (q != q_id) {
          NumericVector gamma_q = gamma[q];
          int group_set = groups(i, q);
          sum_term3 += gamma_q[group_set];
        }
      }
    }
  }
  
  //Calculate params[1]
  params[1] = (sigma2 * sigma2_gamma[q_id]) / (sigma2_gamma[q_id] * sum_term + sigma2);
  
  //Calculate params[0]
  params[0] = params[1] * (-1.0 / (2.0 * sigma2)) * (-2.0 * sum_term2 + 2.0 * sum_term2a +
                                                     2.0 * sum_term2b + 2.0 * sum_term3 +
                                                     2.0 * sum_term4);
  
  //Calculate loglik
  out = rnorm(1, params[0], sqrt(params[1]))[0];
  return out;
}

// [[Rcpp::export]]
double sample_sigma2_scalar(NumericVector y, List x, NumericMatrix groups,
                     NumericMatrix beta, List gamma, NumericVector delta,
                     NumericMatrix z, double alpha, double sigma2,
                     List bs_beta, double phi1, double psi1) {
  double out;
  NumericVector params(2);
    
  //Calculate params[0]
  params[0] = ((double) y.size() / 2.0) + phi1;
    
  //Calculate mean value for each observation
  double sum_sq = 0;
  for (int i = 0; i < y.size(); i++) {
    double mu = alpha;
    for (int k = 0; k < beta.nrow(); k++) {
      NumericMatrix bs_beta_k = bs_beta[k];
      NumericMatrix x_k = x[k];
      for (int j = 0; j < x_k.ncol(); j++) {
        double mu_tmp = 0;
        for (int jj = 0; jj < beta.ncol(); jj++) {
          mu_tmp += bs_beta_k(j, jj) * beta(k, jj);
        }
        mu += x_k(i, j) * mu_tmp;
      }
    }
    for (int k = 1; k < z.ncol(); k++) {
      mu += delta[k] * z(i, k);
    }
    for (int q = 0; q < gamma.size(); q++) {
      NumericVector gamma_q = gamma[q];
      int group_set = groups(i, q);
      mu += gamma_q(group_set);
    }
    sum_sq += (y(i) - mu) * (y(i) - mu);
  }
  
  //Calculate params[1]
  params[1] = psi1 + 0.5 * sum_sq;
    
  //Calculate loglik
  out = 1.0 / (rgamma(1, params[0], 1.0 / params[1])[0]);
  return out;
}

// [[Rcpp::export]]
double sample_sigma2_gamma_scalar(List gamma, double phi2, double psi2, int q_id) {
  double out;
  NumericVector params(2);
    
  //Calculate params[0]
  NumericVector gamma_q = gamma[q_id];
  params[0] = ((double) gamma_q.size() / 2.0) + phi2;
    
  //Calculate params[1]
  double sum_term = 0;
  for (int id = 0; id < gamma_q.size(); id++) {
    sum_term += gamma_q(id) * gamma_q(id);
  }
  params[1] = psi2 + 0.5 * sum_term;
    
  //Calculate loglik
  out = 1.0 / (rgamma(1, params[0], 1.0 / params[1])[0]);
  return out;
}

// [[Rcpp::export]]
NumericMatrix calc_bs_scalar(NumericVector x, NumericVector knots, int degree, NumericVector boundary_knots) {
  int nk = knots.size();
  int len_x = x.size();
  double x_val, term, saved;
  int order = degree + 1;
  NumericMatrix basis(len_x, nk + degree);
  NumericMatrix basis_j(len_x, order);
  int kk_size = 2 * order + nk;
  NumericVector kk(kk_size);
  NumericVector rdel(degree);
  NumericVector ldel(degree);
  
  for (int k = 0; k < order; k++) {
    kk[k] = boundary_knots[0];
    kk[k + order + nk] = boundary_knots[1];
  }
  for (int k = 0; k < nk; k++) {
    kk[k + order] = knots[k];
  }
  
  for (int j = 0; j < len_x; j++) {
    x_val = x[j];
    
    //calculate cursor point
    int curs = -1;
    for (int k = 0; k < kk_size; k++) {
      if (kk[k] >= x_val) curs = k;
      if (kk[k] > x_val) break;
    }
    if (curs > (kk_size - order)) {
      int lastLegit = kk_size - order;
      if (x_val == kk[lastLegit]) {
        curs = lastLegit;
      }
    }
    
    //calculate rdel and ldel
    for (int k = 0; k < degree; k++) {
      rdel[k] = kk[curs + k] - x_val;
      ldel[k] = x_val - kk[curs - k - 1];
    }
    
    basis_j(j, 0) = 1.0;
    for (int k = 1; k < order; k++) {
      saved = 0.0;
      for (int r = 0; r < k; r++) {
        term = basis_j(j, r) / (rdel[r] + ldel[k - 1 - r]);
        basis_j(j, r) = saved + rdel[r] * term;
        saved = ldel[k - 1 - r] * term;
      }
      basis_j(j, k) = saved;
    }
    
    //convert basis_j to the correct size matrix
    int offset = curs - order - 1;
    for (int k = 0; k < order; k++) {
      if ((k + offset) >= 0) {
        basis(j, k + offset) = basis_j(j, k);
      }
    }
  }
  
  return(basis);
}