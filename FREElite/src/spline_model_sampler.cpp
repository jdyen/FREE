// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// slice sampler for spline function regression model
//
// Jian Yen
// Code adapted from https://github.com/mplatzer/BTYDplus/blob/master/src/slice-sampling.cpp
//  and Radford M. Neal 2003. Slice sampling. Annals of Statistics 31: 705-767.
//
// This is a simple slice sampler and the conditional densities needed for a Gibbs sampler update
//  of all parameters in a B-spline function regression model.


double slice_sample(double (*logfn)(double x0, IntegerVector ids, const mat& y,
                                    const mat& x, NumericMatrix groups,
                                    mat& beta, List gamma,
                                    double sigma2, double rho, const mat& b_splines_mat,
                                    double s2_beta, int n, IntegerVector n_j, int n_k, int n_q,
                                    NumericVector sigma2_gamma, List bin_id),
                    double x0, IntegerVector ids, const mat& y, const mat& x,
                    NumericMatrix groups, mat& beta, List gamma,
                    double sigma2, double rho, const mat& b_splines_mat,
                    double s2_beta, int n, IntegerVector n_j, int n_k, int n_q,
                    NumericVector sigma2_gamma, List bin_id,
                    double w_size,
                    double lower,
                    double upper) {
  double u, r0, r1, logy, logz, logys;
  double x1, xs, L, R;
  x1 = x0;
  L = x0;
  R = x0;
  logy = logfn(x1, ids, y, x, groups, beta, gamma,
               sigma2, rho, b_splines_mat, s2_beta,
               n, n_j, n_k, n_q, sigma2_gamma, bin_id);
  
  //draw from [0, y]
  logz = logy - rexp(1)[0];
  
  //expand search range
  u = runif(1)[0] * w_size;
  L = x1 - u;
  R = x1 + (w_size - u);
  while ( L > lower && logfn(L, ids, y, x, groups, beta, gamma,
                             sigma2, rho, b_splines_mat, s2_beta,
                             n, n_j, n_k, n_q, sigma2_gamma, bin_id) > logz )
    L = L - w_size;
  
  while ( R < upper && logfn(R, ids, y, x, groups, beta, gamma,
                             sigma2, rho, b_splines_mat, s2_beta,
                             n, n_j, n_k, n_q, sigma2_gamma, bin_id) > logz )
    R = R + w_size;
  
  //sample until draw is in the correct range
  r0 = std::max(L, lower);
  r1 = std::min(R, upper);
  
  xs = x1;
  int cnt = 0;
  do {
    cnt++;
    xs = runif(1, r0, r1)[0];
    logys = logfn(xs, ids, y, x, groups, beta, gamma,
                  sigma2, rho, b_splines_mat, s2_beta,
                  n, n_j, n_k, n_q, sigma2_gamma, bin_id);
    if ( logys > logz ) {
      break;
    }
    if ( xs < x1 )
      r0 = xs;
    else
      r1 = xs;
    if ((r1 - r0) < 0.0001) {
      xs = r0;
      break;
    }
    if (cnt > 100) {
      Rcout << "slice sampler reached " << cnt << " iterations..." << std::endl;
    }
  } while (cnt < 1e4);
  if (cnt == 1e4) ::Rf_error("slice_sample reached the maximum number of iterations");
  
  x1 = xs;
  logy = logys;
  
  return x1;
}

// [[Rcpp::export]]
double lnL(const mat& y, const mat& x, NumericMatrix groups,
           mat& beta, List gamma, double sigma2,
           double rho, const mat& b_splines_mat,
           double s2_beta, int n, IntegerVector n_j,
           int n_k, int n_q, List bin_id) {
  double loglik;
  vec part_ll(5);
  
  //Calculate part_ll[0]
  part_ll[0] = -((double) n / 2.0) * log (2.0 * PI * sigma2);
  
  //Calculate part_ll[1]
  part_ll[1] = -((double) n / 2.0) * log(1.0 - (rho * rho));
  
  //Calculate mean value for each observation
  mat mu2 = zeros<mat>(n, max(n_j));
  for (int ii = 0; ii < n; ii++) {
    IntegerVector bin_id_i = bin_id[ii];
    for (int j = 0; j < n_j[ii]; j++) {
      mu2(ii, j) += as_scalar((b_splines_mat.row(bin_id_i[j] - 1) * trans(beta)) * trans(x(ii, span::all)));
    }
  }
  for (int q = 0; q < n_q; q++) {
    mat gamma_q = gamma[q];
    for (int ii = 0; ii < n; ii++) {
      IntegerVector bin_id_i = bin_id[ii];
      int group_ID = (int) groups(ii, q);
      for (int j = 0; j < n_j[ii]; j++) {
        mu2(ii, j) += as_scalar(b_splines_mat.row(bin_id_i[j] - 1) *
                                trans(gamma_q.row(group_ID - 1)));
      }
    }
  }
  
  double devi_sq = 0.0;
  double devi_sq2 = 0.0;
  for (int ii = 0; ii < n; ii++) {
    devi_sq2 += (y(ii, 0) - mu2(ii, 0)) * (y(ii, 0) - mu2(ii, 0));
    for (int j = 1; j < n_j[ii]; j++) {
      devi_sq += (y(ii, j) - mu2(ii, j) - rho * (y(ii, j - 1) - mu2(ii, j - 1))) *
      (y(ii, j) - mu2(ii, j) - rho * (y(ii, j - 1) - mu2(ii, j - 1)));
    }
  }
  
  part_ll[2] = -(1.0 / (2.0 * sigma2)) * (1.0 - (rho * rho)) * devi_sq2;
  part_ll[3] = -(((double) sum(n_j)) / 2.0) * log(2.0 * PI * sigma2);
  part_ll[4] = -(1.0 / (2.0 * sigma2)) * devi_sq;
  
  //Calculate loglik
  loglik = sum(part_ll);
  return loglik;
}

// [[Rcpp::export]]
double ln_beta_dens(double x0, IntegerVector ids, const mat& y,
                    const mat& x, NumericMatrix groups, mat& beta,
                    List gamma, double sigma2, double rho,
                    const mat& b_splines_mat, double s2_beta,
                    int n, IntegerVector n_j, int n_k, int n_q,
                    NumericVector sigma2_gamma, List bin_id) {
  double loglik;
  vec part_ll(3);
  int k_id = ids[0];
  int p_id = ids[1];
  
  //Calculate mean value for each observation
  beta(k_id - 1, p_id - 1) = x0;
  mat mu2 = zeros<mat>(n, max(n_j));
  for (int ii = 0; ii < n; ii++) {
    IntegerVector bin_id_i = bin_id[ii];
    for (int j = 0; j < n_j[ii]; j++) {
      mu2(ii, j) += as_scalar((b_splines_mat.row(bin_id_i[j] - 1) * trans(beta)) * trans(x(ii, span::all)));
    }
  }
  for (int q = 0; q < n_q; q++) {
    mat gamma_q = gamma[q];
    for (int ii = 0; ii < n; ii++) {
      IntegerVector bin_id_i = bin_id[ii];
      int group_ID = (int) groups(ii, q);
      for (int j = 0; j < n_j[ii]; j++) {
        mu2(ii, j) += as_scalar(b_splines_mat.row(bin_id_i[j] - 1) *
                                trans(gamma_q.row(group_ID - 1)));
      }
    }
  }
  
  double devi_sq = 0.0;
  double devi_sq2 = 0.0;
  for (int ii = 0; ii < n; ii++) {
    devi_sq2 += (y(ii, 0) - mu2(ii, 0)) * (y(ii, 0) - mu2(ii, 0));
    for (int j = 1; j < n_j[ii]; j++) {
      devi_sq += (y(ii, j) - mu2(ii, j) - rho * (y(ii, j - 1) - mu2(ii, j - 1))) *
      (y(ii, j) - mu2(ii, j) - rho * (y(ii, j - 1) - mu2(ii, j - 1)));
    }
  }
  part_ll[0] = -(1.0 / (2.0 * sigma2)) * (1.0 - (rho * rho)) * devi_sq2;
  part_ll[1] = -(1.0 / (2.0 * sigma2)) * devi_sq;
  
  //Calculate part_ll[2]
  double beta_sum = as_scalar(accu(beta % beta));
  part_ll[2] = -(1.0 / (2.0 * s2_beta)) * beta_sum;
  
  //Calculate loglik
  loglik = sum(part_ll);
  return loglik;
}

// [[Rcpp::export]]
double slice_sample_beta(double x0, IntegerVector ids, const mat& y,
                         const mat& x, NumericMatrix groups,
                         mat& beta, List gamma, double sigma2, double rho,
                         const mat& b_splines_mat, double s2_beta,
                         int n, IntegerVector n_j, int n_k, int n_q,
                         NumericVector sigma2_gamma, List bin_id) {
  return slice_sample(ln_beta_dens, x0, ids, y, x, groups, beta, gamma,
                      sigma2, rho, b_splines_mat, s2_beta, n, n_j, n_k,
                      n_q, sigma2_gamma, bin_id, 0.5, -5.0, 5.0);
}

// [[Rcpp::export]]
double ln_gamma_dens(double x0, IntegerVector ids, const mat& y,
                     const mat& x, NumericMatrix groups,
                     mat& beta, List gamma, double sigma2, double rho,
                     const mat& b_splines_mat, double s2_beta,
                     int n, IntegerVector n_j, int n_k, int n_q,
                     NumericVector sigma2_gamma, List bin_id) {
  double loglik;
  vec part_ll(3);
  int q_id = ids[0];
  int t_id = ids[1];
  int G_id = ids[2];
  mat mu2 = zeros<mat>(n, max(n_j));
  for (int ii = 0; ii < n; ii++) {
    IntegerVector bin_id_i = bin_id[ii];
    for (int j = 0; j < n_j[ii]; j++) {
      mu2(ii, j) += as_scalar((b_splines_mat.row(bin_id_i[j] - 1) * trans(beta)) * trans(x(ii, span::all)));
    }
  }
  for (int q = 0; q < n_q; q++) {
    mat gamma_q = gamma[q];
    if (q == (q_id - 1)) {
      gamma_q(G_id - 1, t_id - 1) = x0;
    }
    for (int ii = 0; ii < n; ii++) {
      IntegerVector bin_id_i = bin_id[ii];
      int group_ID = (int) groups(ii, q);
      for (int j = 0; j < n_j[ii]; j++) {
        mu2(ii, j) += as_scalar(b_splines_mat.row(bin_id_i[j] - 1) *
                                trans(gamma_q.row(group_ID - 1)));
      }
    }
  }
  
  double devi_sq = 0.0;
  double devi_sq2 = 0.0;
  for (int ii = 0; ii < n; ii++) {
    devi_sq2 += (y(ii, 0) - mu2(ii, 0)) * (y(ii, 0) - mu2(ii, 0));
    for (int j = 1; j < n_j[ii]; j++) {
      devi_sq += (y(ii, j) - mu2(ii, j) - rho * (y(ii, j - 1) - mu2(ii, j - 1))) *
      (y(ii, j) - mu2(ii, j) - rho * (y(ii, j - 1) - mu2(ii, j - 1)));
    }
  }
  part_ll[0] = -(1.0 / (2.0 * sigma2)) * (1.0 - (rho * rho)) * devi_sq2;
  part_ll[1] = -(1.0 / (2.0 * sigma2)) * devi_sq;
  
  //Calculate part_ll[2]
  mat gamma_q = gamma[q_id - 1];
  double gamma_sum = as_scalar(accu(gamma_q % gamma_q));
  part_ll[2] = -(1.0 / (2.0 * sigma2_gamma[q_id - 1])) * gamma_sum;
  
  //Calculate loglik
  loglik = sum(part_ll);
  return loglik;
}

// [[Rcpp::export]]
double slice_sample_gamma(double x0, IntegerVector ids, const mat& y,
                          const mat& x, NumericMatrix groups,
                          mat& beta, List gamma, double sigma2, double rho,
                          const mat& b_splines_mat, double s2_beta,
                          int n, IntegerVector n_j, int n_k, int n_q,
                          NumericVector sigma2_gamma, List bin_id) {
  return slice_sample(ln_gamma_dens, x0, ids, y, x, groups, beta, gamma,
                      sigma2, rho, b_splines_mat, s2_beta, n, n_j, n_k,
                      n_q, sigma2_gamma, bin_id, 0.5, -1.0, 1.0);
}

// [[Rcpp::export]]
double ln_rho_dens(double x0, IntegerVector ids, const mat& y,
                   const mat& x, NumericMatrix groups,
                   mat& beta, List gamma, double sigma2, double rho,
                   const mat& b_splines_mat, double s2_beta,
                   int n, IntegerVector n_j, int n_k, int n_q,
                   NumericVector sigma2_gamma, List bin_id) {
  double loglik;
  vec part_ll(3);
  rho = x0;
  
  //Calculate part_ll[0]
  part_ll[0] = -((double) n / 2.0) * log(1.0 - (rho * rho));
  
  //Calculate mean value for each observation
  mat mu2 = zeros<mat>(n, max(n_j));
  for (int ii = 0; ii < n; ii++) {
    IntegerVector bin_id_i = bin_id[ii];
    for (int j = 0; j < n_j[ii]; j++) {
      mu2(ii, j) += as_scalar((b_splines_mat.row(bin_id_i[j] - 1) * trans(beta)) * trans(x(ii, span::all)));
    }
  }
  for (int q = 0; q < n_q; q++) {
    mat gamma_q = gamma[q];
    for (int ii = 0; ii < n; ii++) {
      IntegerVector bin_id_i = bin_id[ii];
      int group_ID = (int) groups(ii, q);
      for (int j = 0; j < n_j[ii]; j++) {
        mu2(ii, j) += as_scalar(b_splines_mat.row(bin_id_i[j] - 1) *
                                trans(gamma_q.row(group_ID - 1)));
      }
    }
  }
  
  double devi_sq = 0.0;
  double devi_sq2 = 0.0;
  for (int ii = 0; ii < n; ii++) {
    devi_sq2 += (y(ii, 0) - mu2(ii, 0)) * (y(ii, 0) - mu2(ii, 0));
    for (int j = 1; j < n_j[ii]; j++) {
      devi_sq += (y(ii, j) - mu2(ii, j) - rho * (y(ii, j - 1) - mu2(ii, j - 1))) *
      (y(ii, j) - mu2(ii, j) - rho * (y(ii, j - 1) - mu2(ii, j - 1)));
    }
  }
  part_ll[1] = -(1.0 / (2.0 * sigma2)) * (1.0 - (rho * rho)) * devi_sq2;
  part_ll[2] = -(1.0 / (2.0 * sigma2)) * devi_sq;
  
  //Calculate loglik
  if (rho < 1 && rho > -1) {
    loglik = sum(part_ll);
  } else {
    loglik = -INFINITY;
  }
  return loglik;
}

// [[Rcpp::export]]
double slice_sample_rho(double x0, IntegerVector ids, const mat& y,
                        const mat& x, NumericMatrix groups,
                        mat& beta, List gamma, double sigma2, double rho,
                        const mat& b_splines_mat, double s2_beta,
                        int n, IntegerVector n_j, int n_k, int n_q,
                        NumericVector sigma2_gamma, List bin_id) {
  return slice_sample(ln_rho_dens, x0, ids, y, x, groups, beta, gamma,
                      sigma2, rho, b_splines_mat, s2_beta, n, n_j, n_k,
                      n_q, sigma2_gamma, bin_id, 0.2, -1.0, 1.0);
}

// [[Rcpp::export]]
NumericMatrix calc_bs(NumericVector x, NumericVector knots,
                      int degree, NumericVector boundary_knots) {
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

// [[Rcpp::export]]
double sample_sigma2(const mat& y, const mat& x, NumericMatrix groups,
                     mat& beta, List gamma, double sigma2, double rho,
                     const mat& b_splines_mat, double s2h_a, double s2h_b,
                     int n, IntegerVector n_j, int n_k, int n_q, List bin_id) {
  double out;
  vec params(2);
  
  //Calculate params[0]
  params[0] = ((double) n / 2.0) + (((double) sum(n_j)) / 2.0) + s2h_a;
  
  //Calculate mean value for each observation
  mat mu2 = zeros<mat>(n, max(n_j));
  for (int ii = 0; ii < n; ii++) {
    IntegerVector bin_id_i = bin_id[ii];
    for (int j = 0; j < n_j[ii]; j++) {
      mu2(ii, j) += as_scalar((b_splines_mat.row(bin_id_i[j] - 1) * trans(beta)) * trans(x(ii, span::all)));
    }
  }
  for (int q = 0; q < n_q; q++) {
    mat gamma_q = gamma[q];
    for (int ii = 0; ii < n; ii++) {
      IntegerVector bin_id_i = bin_id[ii];
      int group_ID = (int) groups(ii, q);
      for (int j = 0; j < n_j[ii]; j++) {
        mu2(ii, j) += as_scalar(b_splines_mat.row(bin_id_i[j] - 1) *
                                trans(gamma_q.row(group_ID - 1)));
      }
    }
  }
  
  double devi_sq = 0.0;
  double devi_sq2 = 0.0;
  for (int ii = 0; ii < n; ii++) {
    devi_sq2 += (y(ii, 0) - mu2(ii, 0)) * (y(ii, 0) - mu2(ii, 0));
    for (int j = 1; j < n_j[ii]; j++) {
      devi_sq += (y(ii, j) - mu2(ii, j) - rho * (y(ii, j - 1) - mu2(ii, j - 1))) *
      (y(ii, j) - mu2(ii, j) - rho * (y(ii, j - 1) - mu2(ii, j - 1)));
    }
  }
  params[1] = s2h_b + 0.5 * (1.0 - (rho * rho)) * devi_sq2 + 0.5 * devi_sq;
  
  //Calculate loglik
  out = 1.0 / (rgamma(1, params[0], 1.0 / params[1])[0]);
  return out;
}

// [[Rcpp::export]]
double sample_sigma2_gamma(List gamma, double s2hg_a, double s2hg_b,
                           int n_q, int q_id) {
  double out;
  
  //Calculate params
  mat gamma_q = gamma[q_id - 1];
  double gamma_sum = as_scalar(accu(gamma_q % gamma_q));
  double num_gammas = as_scalar(gamma_q.n_cols * gamma_q.n_rows);
  
  double params_a = s2hg_a + ((double) num_gammas / 2.0);
  double params_b = s2hg_b + 0.5 * gamma_sum;
  
  //Calculate loglik
  out = 1.0 / (rgamma(1, params_a, params_b)[0]);
  return out;
}