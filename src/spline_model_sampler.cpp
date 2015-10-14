#include <Rcpp.h>

using namespace Rcpp;

// slice sampler for spline function regression model
//
// Jian Yen
// Code adapted from https://github.com/mplatzer/BTYDplus/blob/master/src/slice-sampling.cpp
//  and Radford M. Neal 2003. Slice sampling. Anals of Statistics 31: 705-767.
//
// This is a simple slice sampler and the conditional densities needed for a Gibbs sampler update
//  of all parameters in a B-spline function regression model.


double slice_sample(double (*logfn)(double x0, IntegerVector ids, NumericMatrix y,
                                    NumericMatrix x, NumericVector w, NumericMatrix groups,
                                    List beta, List gamma, List theta1, List theta2,
                                    double sigma2, double rho, List b_splines_beta,
                                    List b_splines_gamma, double s2h_a, double s2h_b,
                                    double s2hg_a, double s2hg_b, double s2_beta,
                                    int n, int n_j, int n_k, int n_q,
                                    NumericVector sigma2_gamma),
                    double x0, IntegerVector ids, NumericMatrix y,
                    NumericMatrix x, NumericVector w, NumericMatrix groups,
                    List beta, List gamma, List theta1, List theta2,
                    double sigma2, double rho, List b_splines_beta,
                    List b_splines_gamma, double s2h_a, double s2h_b,
                    double s2hg_a, double s2hg_b, double s2_beta,
                    int n, int n_j, int n_k, int n_q,
                    NumericVector sigma2_gamma,
                    double w_size,
                    double lower,
                    double upper) {
  double u, r0, r1, logy, logz, logys;
  double x1, xs, L, R;
  x1 = x0;
  L = x0;
  R = x0;
  logy = logfn(x1, ids, y, x, w, groups, beta, gamma, theta1, theta2,
               sigma2, rho, b_splines_beta, b_splines_gamma,
               s2h_a, s2h_b, s2hg_a, s2hg_b, s2_beta, n,
               n_j, n_k, n_q, sigma2_gamma);
  
  //draw from [0, y]
  logz = logy - rexp(1)[0];
    
  //expand search range
  u = runif(1)[0] * w_size;
  L = x1 - u;
  R = x1 + (w_size - u);
  while ( L > lower && logfn(L, ids, y, x, w, groups, beta, gamma, theta1, theta2,
                             sigma2, rho, b_splines_beta, b_splines_gamma,
                             s2h_a, s2h_b, s2hg_a, s2hg_b, s2_beta, n,
                             n_j, n_k, n_q, sigma2_gamma) > logz )
    L = L - w_size;
  while ( R < upper && logfn(R, ids, y, x, w, groups, beta, gamma, theta1, theta2,
                             sigma2, rho, b_splines_beta, b_splines_gamma,
                             s2h_a, s2h_b, s2hg_a, s2hg_b, s2_beta, n,
                             n_j, n_k, n_q, sigma2_gamma) > logz )
    R = R + w_size;
      
  //sample until draw is in the correct range
  r0 = std::max(L, lower);
  r1 = std::min(R, upper);
    
  xs = x1;
  int cnt = 0;
  do {
    cnt++;
    xs = runif(1, r0, r1)[0];
    logys = logfn(xs, ids, y, x, w, groups, beta, gamma, theta1, theta2,
                  sigma2, rho, b_splines_beta, b_splines_gamma,
                  s2h_a, s2h_b, s2hg_a, s2hg_b, s2_beta, n,
                  n_j, n_k, n_q, sigma2_gamma);
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
double lnL(NumericMatrix y, NumericMatrix x, NumericVector w, NumericMatrix groups,
           List beta, List gamma, List theta1, List theta2, double sigma2,
           double rho, List b_splines_beta, List b_splines_gamma,
           double s2h_a, double s2h_b, double s2hg_a, double s2hg_b,
           double s2_beta, int n, int n_j, int n_k, int n_q) {
  double loglik;
  NumericVector part_ll(5);
  
  //Calculate part_ll[0]
  part_ll[0] = -((double) n / 2.0) * log (2.0 * PI * sigma2);
  
  //Calculate part_ll[1]
  part_ll[1] = -((double) n / 2.0) * log(1.0 - (rho * rho));
  
  //Calculate mean value for each observation
  NumericMatrix mu(n, n_j);
  for (int k = 0; k < n_k; k++) {
    NumericMatrix b_splines_k = b_splines_beta[k];
    NumericVector beta_k = beta[k];
    for (int j = 0; j < n_j; j++) {
      for (int ii = 0; ii < n; ii++) {
        double mu_tmp = 0;
        for (int jj = 0; jj < beta_k.size(); jj++) {
          mu_tmp += b_splines_k(j, jj) * beta_k[jj];
        }
        mu(ii, j) += x(ii, k) * mu_tmp;
      }
    }
  }
  for (int q = 0; q < n_q; q++) {
    List b_splines_q = b_splines_gamma[q];
    List gamma_q = gamma[q];
    for (int j = 0; j < n_j; j++) {
      for (int ii = 0; ii < n; ii++) {
        int group_ID = (int) groups(ii, q);
        NumericMatrix bs_qG = b_splines_q(group_ID - 1);
        NumericVector gamma_q_G = gamma_q[group_ID - 1];
        for (int jj = 0; jj < gamma_q_G.size(); jj++) {
          mu(ii, j) += bs_qG(j, jj) * gamma_q_G[jj];
        }
      }
    }
  }
  
  //Calculate part_ll[2]
  part_ll[2] = -(1.0 / (2.0 * sigma2)) * (1.0 - (rho * rho)) *
  sum((y(_, 0) - mu(_, 0)) * (y(_, 0) - mu(_, 0)));
  
  //Calculate part_ll[3]
  part_ll[3] = -(((double) n_j * (double) n) / 2.0) * log(2.0 * PI * sigma2);
  
  //Calculate part_ll[4]
  double devi_sq = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < (n_j - 1); j++) {
      devi_sq += (y(i, j + 1) - mu(i, j + 1) - rho * (y(i, j) - mu(i, j))) *
      (y(i, j + 1) - mu(i, j + 1) - rho * (y(i, j) - mu(i, j)));
    }
  }
  part_ll[4] = -(1.0 / (2.0 * sigma2)) * devi_sq;
  
  //Calculate loglik
  loglik = sum(part_ll);
  return loglik;
}

// [[Rcpp::export]]
double ln_beta_dens(double x0, IntegerVector ids, NumericMatrix y,
                    NumericMatrix x, NumericVector w, NumericMatrix groups,
                    List beta, List gamma, List theta1, List theta2,
                    double sigma2, double rho, List b_splines_beta,
                    List b_splines_gamma, double s2h_a, double s2h_b,
                    double s2hg_a, double s2hg_b, double s2_beta,
                    int n, int n_j, int n_k, int n_q,
                    NumericVector sigma2_gamma) {
  double loglik;
  NumericVector part_ll(3);
  int k_id = ids[0];
  int p_id = ids[1];
  
  //Calculate mean value for each observation
  NumericMatrix mu(n, n_j);
  for (int k = 0; k < n_k; k++) {
    NumericMatrix b_splines_k = b_splines_beta[k];
    NumericVector beta_k = beta[k];
    if (k == (k_id - 1))
      beta_k[p_id - 1] = x0;
    for (int j = 0; j < n_j; j++) {
      for (int ii = 0; ii < n; ii++) {
        double mu_tmp = 0;
        for (int jj = 0; jj < beta_k.size(); jj++) {
          mu_tmp += b_splines_k(j, jj) * beta_k[jj];
        }
        mu(ii, j) += x(ii, k) * mu_tmp;
      }
    }
  }
  for (int q = 0; q < n_q; q++) {
    List b_splines_q = b_splines_gamma[q];
    List gamma_q = gamma[q];
    for (int j = 0; j < n_j; j++) {
      for (int ii = 0; ii < n; ii++) {
        int group_ID = (int) groups(ii, q);
        NumericMatrix bs_qG = b_splines_q[group_ID - 1];
        NumericVector gamma_q_G = gamma_q[group_ID - 1];
        for (int jj = 0; jj < gamma_q_G.size(); jj++) {
          mu(ii, j) += bs_qG(j, jj) * gamma_q_G[jj];
        }
      }
    }
  }
  
  //Calculate part_ll[0]
  part_ll[0] = -(1.0 / (2.0 * sigma2)) * (1.0 - (rho * rho)) *
  sum((y(_, 0) - mu(_, 0)) * (y(_, 0) - mu(_, 0)));
  
  //Calculate part_ll[1]
  double devi_sq = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < (n_j - 1); j++) {
      devi_sq += (y(i, j + 1) - mu(i, j + 1) - rho * (y(i, j) - mu(i, j))) *
      (y(i, j + 1) - mu(i, j + 1) - rho * (y(i, j) - mu(i, j)));
    }
  }
  part_ll[1] = -(1.0 / (2.0 * sigma2)) * devi_sq;
  
  //Calculate part_ll[2]
  double beta_sum = 0;
  for (int k = 0; k < n_k; k++) {
    NumericVector beta_k = beta(k);
    if (k == (k_id - 1))
      beta_k[p_id - 1] = x0;
    beta_sum += sum(beta_k * beta_k);
  }
  part_ll[2] = -(1.0 / (2.0 * s2_beta)) * beta_sum;
  
  //Calculate loglik
  loglik = sum(part_ll);
  return loglik;
}

// [[Rcpp::export]]
double slice_sample_beta(double x0, IntegerVector ids, NumericMatrix y,
                         NumericMatrix x, NumericVector w, NumericMatrix groups,
                         List beta, List gamma, List theta1, List theta2,
                         double sigma2, double rho, List b_splines_beta,
                         List b_splines_gamma, double s2h_a, double s2h_b,
                         double s2hg_a, double s2hg_b, double s2_beta,
                         int n, int n_j, int n_k, int n_q,
                         NumericVector sigma2_gamma) {
  return slice_sample(ln_beta_dens, x0, ids, y, x, w, groups, beta, gamma,
                      theta1, theta2, sigma2, rho, b_splines_beta,
                      b_splines_gamma, s2h_a, s2h_b, s2hg_a, s2hg_b,
                      s2_beta, n, n_j, n_k, n_q, sigma2_gamma,
                      10.0, -1e6, 1e6);
}

// [[Rcpp::export]]
double ln_gamma_dens(double x0, IntegerVector ids, NumericMatrix y,
                     NumericMatrix x, NumericVector w, NumericMatrix groups,
                     List beta, List gamma, List theta1, List theta2,
                     double sigma2, double rho, List b_splines_beta,
                     List b_splines_gamma, double s2h_a, double s2h_b,
                     double s2hg_a, double s2hg_b, double s2_beta,
                     int n, int n_j, int n_k, int n_q,
                     NumericVector sigma2_gamma) {
  double loglik;
  NumericVector part_ll(3);
  int q_id = ids[0];
  int t_id = ids[1];
  int G_id = ids[2];
  
  //Calculate mean value for each observation
  NumericMatrix mu(n, n_j);
  for (int k = 0; k < n_k; k++) {
    NumericMatrix b_splines_k = b_splines_beta(k);
    NumericVector beta_k = beta(k);
    for (int j = 0; j < n_j; j++) {
      for (int ii = 0; ii < n; ii++) {
        double mu_tmp = 0;
        for (int jj = 0; jj < beta_k.size(); jj++) {
          mu_tmp += b_splines_k(j, jj) * beta_k(jj);
        }
        mu(ii, j) += x(ii, k) * mu_tmp;
      }
    }
  }
  for (int q = 0; q < n_q; q++) {
    List b_splines_q = b_splines_gamma(q);
    List gamma_q = gamma(q);
    for (int j = 0; j < n_j; j++) {
      for (int ii = 0; ii < n; ii++) {
        int group_ID = (int) groups(ii, q);
        NumericMatrix bs_qG = b_splines_q(group_ID - 1);
        NumericVector gamma_q_G = gamma_q(group_ID - 1);
        if ((q == (q_id - 1)) && (group_ID == G_id))
          gamma_q_G(t_id - 1) = x0;
        for (int jj = 0; jj < gamma_q_G.size(); jj++) {
          mu(ii, j) += bs_qG(j, jj) * gamma_q_G(jj);
        }
      }
    }
  }
  
  //Calculate part_ll[0]
  part_ll[0] = -(1.0 / (2.0 * sigma2)) * (1.0 - (rho * rho)) *
  sum((y(_, 0) - mu(_, 0)) * (y(_, 0) - mu(_, 0)));
  
  //Calculate part_ll[1]
  double devi_sq = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < (n_j - 1); j++) {
      devi_sq += (y(i, j + 1) - mu(i, j + 1) - rho * (y(i, j) - mu(i, j))) *
      (y(i, j + 1) - mu(i, j + 1) - rho * (y(i, j) - mu(i, j)));
    }
  }
  part_ll[1] = -(1.0 / (2.0 * sigma2)) * devi_sq;
  
  //Calculate part_ll[2]
  double gamma_sum = 0;
  List gamma_q = gamma[q_id - 1];
  for (int i = 0; i < gamma_q.size(); i++) {
    NumericVector gamma_qG = gamma_q[i];
    if (i == (G_id - 1))
      gamma_qG[t_id - 1] = x0;
    gamma_sum += sum(gamma_qG * gamma_qG);
  }
  double s2g_q = sigma2_gamma[q_id - 1];
  part_ll[2] = -(1.0 / (2.0 * s2g_q)) * gamma_sum;
  
  //Calculate loglik
  loglik = sum(part_ll);
  return loglik;
}

// [[Rcpp::export]]
double slice_sample_gamma(double x0, IntegerVector ids, NumericMatrix y,
                          NumericMatrix x, NumericVector w, NumericMatrix groups,
                          List beta, List gamma, List theta1, List theta2,
                          double sigma2, double rho, List b_splines_beta,
                          List b_splines_gamma, double s2h_a, double s2h_b,
                          double s2hg_a, double s2hg_b, double s2_beta,
                          int n, int n_j, int n_k, int n_q,
                          NumericVector sigma2_gamma) {
  return slice_sample(ln_gamma_dens, x0, ids, y, x, w, groups, beta, gamma,
                      theta1, theta2, sigma2, rho, b_splines_beta,
                      b_splines_gamma, s2h_a, s2h_b, s2hg_a, s2hg_b,
                      s2_beta, n, n_j, n_k, n_q, sigma2_gamma,
                      10.0, -1e6, 1e6);
}

// [[Rcpp::export]]
double ln_rho_dens(double x0, IntegerVector ids, NumericMatrix y,
                   NumericMatrix x, NumericVector w, NumericMatrix groups,
                   List beta, List gamma, List theta1, List theta2,
                   double sigma2, double rho, List b_splines_beta,
                   List b_splines_gamma, double s2h_a, double s2h_b,
                   double s2hg_a, double s2hg_b, double s2_beta,
                   int n, int n_j, int n_k, int n_q,
                   NumericVector sigma2_gamma) {
  double loglik;
  NumericVector part_ll(3);
  rho = x0;
  
  //Calculate part_ll[0]
  part_ll[0] = -((double) n / 2.0) * log(1.0 - (rho * rho));
  
  //Calculate mean value for each observation
  NumericMatrix mu(n, n_j);
  for (int k = 0; k < n_k; k++) {
    NumericMatrix b_splines_k = b_splines_beta[k];
    NumericVector beta_k = beta[k];
    for (int j = 0; j < n_j; j++) {
      for (int ii = 0; ii < n; ii++) {
        double mu_tmp = 0;
        for (int jj = 0; jj < beta_k.size(); jj++) {
          mu_tmp += b_splines_k(j, jj) * beta_k(jj);
        }
        mu(ii, j) += x(ii, k) * mu_tmp;
      }
    }
  }
  for (int q = 0; q < n_q; q++) {
    List b_splines_q = b_splines_gamma[q];
    List gamma_q = gamma[q];
    for (int j = 0; j < n_j; j++) {
      for (int ii = 0; ii < n; ii++) {
        int group_ID = (int) groups(ii, q);
        NumericMatrix bs_qG = b_splines_q[group_ID - 1];
        NumericVector gamma_q_G = gamma_q[group_ID - 1];
        for (int jj = 0; jj < gamma_q_G.size(); jj++) {
          mu(ii, j) += bs_qG(j, jj) * gamma_q_G[jj];
        }
      }
    }
  }
  
  //Calculate part_ll[1]
  part_ll[1] = -(1.0 / (2.0 * sigma2)) * (1.0 - (rho * rho)) *
  sum((y(_, 0) - mu(_, 0)) * (y(_, 0) - mu(_, 0)));
  
  //Calculate part_ll[2]
  double devi_sq = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < (n_j - 1); j++) {
      devi_sq += (y(i, j + 1) - mu(i, j + 1) - rho * (y(i, j) - mu(i, j))) *
      (y(i, j + 1) - mu(i, j + 1) - rho * (y(i, j) - mu(i, j)));
    }
  }
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
double slice_sample_rho(double x0, IntegerVector ids, NumericMatrix y,
                        NumericMatrix x, NumericVector w, NumericMatrix groups,
                        List beta, List gamma, List theta1, List theta2,
                        double sigma2, double rho, List b_splines_beta,
                        List b_splines_gamma, double s2h_a, double s2h_b,
                        double s2hg_a, double s2hg_b, double s2_beta,
                        int n, int n_j, int n_k, int n_q,
                        NumericVector sigma2_gamma) {
  return slice_sample(ln_rho_dens, x0, ids, y, x, w, groups, beta, gamma,
                      theta1, theta2, sigma2, rho, b_splines_beta,
                      b_splines_gamma, s2h_a, s2h_b, s2hg_a, s2hg_b,
                      s2_beta, n, n_j, n_k, n_q, sigma2_gamma,
                      0.2, -1.0, 1.0);
}

// [[Rcpp::export]]
NumericMatrix calc_bs(NumericVector x, NumericVector knots, int degree, NumericVector boundary_knots) {
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
double ln_theta1_dens(double x0, IntegerVector ids, NumericMatrix y,
                      NumericMatrix x, NumericVector w, NumericMatrix groups,
                      List beta, List gamma, List theta1, List theta2,
                      double sigma2, double rho, List b_splines_beta,
                      List b_splines_gamma, double s2h_a, double s2h_b,
                      double s2hg_a, double s2hg_b, double s2_beta,
                      int n, int n_j, int n_k, int n_q,
                      NumericVector sigma2_gamma) {
  double loglik;
  NumericVector part_ll(2);
  int k_id = ids[0];
  int p_id = ids[1];
  
  //Calculate mean value for each observation
  NumericMatrix mu(n, n_j);
  for (int k = 0; k < n_k; k++) {
    NumericVector beta_k = beta[k];
    NumericMatrix b_splines_k = b_splines_beta[k];
    if (k == (k_id - 1)) {
      NumericVector theta1_k = theta1[k_id - 1];
      theta1_k[p_id - 1] = x0;
      b_splines_k = calc_bs(w, theta1_k, 3, NumericVector::create(min(w), max(w)));
    }
    for (int j = 0; j < n_j; j++) {
      for (int ii = 0; ii < n; ii++) {
        double mu_tmp = 0;
        for (int jj = 0; jj < beta_k.size(); jj++) {
          mu_tmp += b_splines_k(j, jj) * beta_k[jj];
        }
        mu(ii, j) += x(ii, k) * mu_tmp;
      }
    }
  }
  for (int q = 0; q < n_q; q++) {
    List b_splines_q = b_splines_gamma[q];
    List gamma_q = gamma[q];
    for (int j = 0; j < n_j; j++) {
      for (int ii = 0; ii < n; ii++) {
        int group_ID = (int) groups(ii, q);
        NumericMatrix bs_qG = b_splines_q[group_ID - 1];
        NumericVector gamma_q_G = gamma_q[group_ID - 1];
        for (int jj = 0; jj < gamma_q_G.size(); jj++) {
          mu(ii, j) += bs_qG(j, jj) * gamma_q_G[jj];
        }
      }
    }
  }
  
  //Calculate part_ll[0]
  part_ll[0] = -(1.0 / (2.0 * sigma2)) * (1.0 - (rho * rho)) *
  sum((y(_, 0) - mu(_, 0)) * (y(_, 0) - mu(_, 0)));
  
  //Calculate part_ll[1]
  double devi_sq = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < (n_j - 1); j++) {
      devi_sq += (y(i, j + 1) - mu(i, j + 1) - rho * (y(i, j) - mu(i, j))) *
      (y(i, j + 1) - mu(i, j + 1) - rho * (y(i, j) - mu(i, j)));
    }
  }
  part_ll[1] = -(1.0 / (2.0 * sigma2)) * devi_sq;

  //Calculate loglik
  if (x0 > min(w) && x0 < max(w)) {
    loglik = sum(part_ll);
  } else {
    loglik = -INFINITY;
  }
  return loglik;
}

// [[Rcpp::export]]
double slice_sample_theta1(double x0, IntegerVector ids, NumericMatrix y,
                           NumericMatrix x, NumericVector w, NumericMatrix groups,
                           List beta, List gamma, List theta1, List theta2,
                           double sigma2, double rho, List b_splines_beta,
                           List b_splines_gamma, double s2h_a, double s2h_b,
                           double s2hg_a, double s2hg_b, double s2_beta,
                           int n, int n_j, int n_k, int n_q,
                           NumericVector sigma2_gamma) {
  return slice_sample(ln_theta1_dens, x0, ids, y, x, w, groups, beta, gamma,
                      theta1, theta2, sigma2, rho, b_splines_beta,
                      b_splines_gamma, s2h_a, s2h_b, s2hg_a, s2hg_b,
                      s2_beta, n, n_j, n_k, n_q, sigma2_gamma,
                      1.0, min(w), max(w));
}

// [[Rcpp::export]]
double ln_theta2_dens(double x0, IntegerVector ids, NumericMatrix y,
                      NumericMatrix x, NumericVector w, NumericMatrix groups,
                      List beta, List gamma, List theta1, List theta2,
                      double sigma2, double rho, List b_splines_beta,
                      List b_splines_gamma, double s2h_a, double s2h_b,
                      double s2hg_a, double s2hg_b, double s2_beta,
                      int n, int n_j, int n_k, int n_q,
                      NumericVector sigma2_gamma) {
  double loglik;
  NumericVector part_ll(2);
  int q_id = ids[0];
  int t_id = ids[1];
  int G_id = ids[2];
  
  //Calculate mean value for each observation
  NumericMatrix mu(n, n_j);
  for (int k = 0; k < n_k; k++) {
    NumericVector beta_k = beta[k];
    NumericMatrix b_splines_k = b_splines_beta[k];
    for (int j = 0; j < n_j; j++) {
      for (int ii = 0; ii < n; ii++) {
        double mu_tmp = 0;
        for (int jj = 0; jj < beta_k.size(); jj++) {
          mu_tmp += b_splines_k(j, jj) * beta_k[jj];
        }
        mu(ii, j) += x(ii, k) * mu_tmp;
      }
    }
  }

  for (int q = 0; q < n_q; q++) {
    List b_splines_q = b_splines_gamma[q];
    List gamma_q = gamma[q];
    if (q == (q_id - 1)) {
      List theta2_q = theta2[q];
      NumericVector theta2_qG = theta2_q[G_id - 1];
      theta2_qG[t_id - 1] = x0;
      b_splines_q[G_id - 1] = calc_bs(w, theta2_qG, 3, NumericVector::create(min(w), max(w)));
    }
    for (int j = 0; j < n_j; j++) {
      for (int ii = 0; ii < n; ii++) {
        int group_ID = (int) groups(ii, q);
        NumericMatrix bs_qG = b_splines_q[group_ID - 1];
        NumericVector gamma_q_G = gamma_q[group_ID - 1];
        for (int jj = 0; jj < gamma_q_G.size(); jj++) {
          mu(ii, j) += bs_qG(j, jj) * gamma_q_G[jj];
        }
      }
    }
  }
  
  //Calculate part_ll[0]
  part_ll[0] = -(1.0 / (2.0 * sigma2)) * (1.0 - (rho * rho)) *
  sum((y(_, 0) - mu(_, 0)) * (y(_, 0) - mu(_, 0)));
  
  //Calculate part_ll[1]
  double devi_sq = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < (n_j - 1); j++) {
      devi_sq += (y(i, j + 1) - mu(i, j + 1) - rho * (y(i, j) - mu(i, j))) *
      (y(i, j + 1) - mu(i, j + 1) - rho * (y(i, j) - mu(i, j)));
    }
  }
  part_ll[1] = -(1.0 / (2.0 * sigma2)) * devi_sq;
  
  //Calculate loglik
  if (x0 > min(w) && x0 < max(w)) {
    loglik = sum(part_ll);
  } else {
    loglik = -INFINITY;
  }
  return loglik;
}

// [[Rcpp::export]]
double slice_sample_theta2(double x0, IntegerVector ids, NumericMatrix y,
                           NumericMatrix x, NumericVector w, NumericMatrix groups,
                           List beta, List gamma, List theta1, List theta2,
                           double sigma2, double rho, List b_splines_beta,
                           List b_splines_gamma, double s2h_a, double s2h_b,
                           double s2hg_a, double s2hg_b, double s2_beta,
                           int n, int n_j, int n_k, int n_q,
                           NumericVector sigma2_gamma) {
  return slice_sample(ln_theta2_dens, x0, ids, y, x, w, groups, beta, gamma,
                      theta1, theta2, sigma2, rho, b_splines_beta,
                      b_splines_gamma, s2h_a, s2h_b, s2hg_a, s2hg_b,
                      s2_beta, n, n_j, n_k, n_q, sigma2_gamma,
                      1.0, min(w), max(w));
}

// [[Rcpp::export]]
double sample_sigma2(NumericMatrix y, NumericMatrix x, NumericVector w, NumericMatrix groups,
                           List beta, List gamma, List theta1, List theta2, double sigma2,
                           double rho, List b_splines_beta, List b_splines_gamma,
                           double s2h_a, double s2h_b, double s2hg_a, double s2hg_b,
                           double s2_beta, int n, int n_j, int n_k, int n_q) {
  double out;
  NumericVector params(2);
  
  //Calculate params[0]
  params[0] = ((double) n / 2.0) + (((double) n_j * (double) n) / 2.0) + s2h_a;
  
  //Calculate mean value for each observation
  NumericMatrix mu(n, n_j);
  for (int k = 0; k < n_k; k++) {
    NumericMatrix b_splines_k = b_splines_beta[k];
    NumericVector beta_k = beta[k];
    for (int j = 0; j < n_j; j++) {
      for (int ii = 0; ii < n; ii++) {
        double mu_tmp = 0;
        for (int jj = 0; jj < beta_k.size(); jj++) {
          mu_tmp += b_splines_k(j, jj) * beta_k[jj];
        }
        mu(ii, j) += x(ii, k) * mu_tmp;
      }
    }
  }

  for (int q = 0; q < n_q; q++) {
    List b_splines_q = b_splines_gamma[q];
    List gamma_q = gamma[q];
    for (int j = 0; j < n_j; j++) {
      for (int ii = 0; ii < n; ii++) {
        int group_ID = (int) groups(ii, q);
        NumericMatrix bs_qG = b_splines_q[group_ID - 1];
        NumericVector gamma_q_G = gamma_q(group_ID - 1);
        for (int jj = 0; jj < gamma_q_G.size(); jj++) {
          mu(ii, j) += bs_qG(j, jj) * gamma_q_G(jj);
        }
      }
    }
  }
  
  //Calculate params[1] part 1/2
  params[1] = s2h_b + 0.5 * (1.0 - (rho * rho)) * sum((y(_, 0) - mu(_, 0)) * (y(_, 0) - mu(_, 0)));
  
  //Calculate params[1] part 2/2
  double devi_sq = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < (n_j - 1); j++) {
      devi_sq += (y(i, j + 1) - mu(i, j + 1) - rho * (y(i, j) - mu(i, j))) *
      (y(i, j + 1) - mu(i, j + 1) - rho * (y(i, j) - mu(i, j)));
    }
  }
  params[1] += 0.5 * devi_sq;
  
  //Calculate loglik
  out = 1.0 / (rgamma(1, params[0], 1.0 / params[1])[0]);
  return out;
}

// [[Rcpp::export]]
double sample_sigma2_gamma(List gamma, double s2hg_a, double s2hg_b, int n_q, int q_id) {
  double out;
  
  //Calculate params
  double gamma_sum = 0;
  int num_gammas = 0;
  List gamma_q = gamma(q_id - 1);
  for (int i = 0; i < gamma_q.size(); i++) {
    NumericVector gamma_qG = gamma_q(i);
    num_gammas += gamma_qG.size();
    gamma_sum += sum(gamma_qG * gamma_qG);
  }
  double params_a = s2hg_a + ((double) num_gammas / 2.0);
  double params_b = s2hg_b + 0.5 * gamma_sum;
  
  //Calculate loglik
  out = 1.0 / (rgamma(1, params_a, params_b)[0]);
  return out;
}