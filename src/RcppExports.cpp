// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// lnL_scalar
double lnL_scalar(NumericVector y, List x, NumericMatrix groups, NumericMatrix beta, List gamma, NumericVector delta, NumericMatrix z, double alpha, double sigma2, List bs_beta);
RcppExport SEXP FREE_lnL_scalar(SEXP ySEXP, SEXP xSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP deltaSEXP, SEXP zSEXP, SEXP alphaSEXP, SEXP sigma2SEXP, SEXP bs_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< List >::type bs_beta(bs_betaSEXP);
    __result = Rcpp::wrap(lnL_scalar(y, x, groups, beta, gamma, delta, z, alpha, sigma2, bs_beta));
    return __result;
END_RCPP
}
// sample_alpha_scalar
double sample_alpha_scalar(NumericVector y, List x, NumericMatrix groups, NumericMatrix beta, List gamma, NumericVector delta, NumericMatrix z, double alpha, double sigma2, List bs_beta, double s2_alpha);
RcppExport SEXP FREE_sample_alpha_scalar(SEXP ySEXP, SEXP xSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP deltaSEXP, SEXP zSEXP, SEXP alphaSEXP, SEXP sigma2SEXP, SEXP bs_betaSEXP, SEXP s2_alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< List >::type bs_beta(bs_betaSEXP);
    Rcpp::traits::input_parameter< double >::type s2_alpha(s2_alphaSEXP);
    __result = Rcpp::wrap(sample_alpha_scalar(y, x, groups, beta, gamma, delta, z, alpha, sigma2, bs_beta, s2_alpha));
    return __result;
END_RCPP
}
// sample_beta_scalar
double sample_beta_scalar(NumericVector y, List x, NumericMatrix groups, NumericMatrix beta, List gamma, NumericVector delta, NumericMatrix z, double alpha, double sigma2, List bs_beta, double s2_beta, int p_id, int k_id);
RcppExport SEXP FREE_sample_beta_scalar(SEXP ySEXP, SEXP xSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP deltaSEXP, SEXP zSEXP, SEXP alphaSEXP, SEXP sigma2SEXP, SEXP bs_betaSEXP, SEXP s2_betaSEXP, SEXP p_idSEXP, SEXP k_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< List >::type bs_beta(bs_betaSEXP);
    Rcpp::traits::input_parameter< double >::type s2_beta(s2_betaSEXP);
    Rcpp::traits::input_parameter< int >::type p_id(p_idSEXP);
    Rcpp::traits::input_parameter< int >::type k_id(k_idSEXP);
    __result = Rcpp::wrap(sample_beta_scalar(y, x, groups, beta, gamma, delta, z, alpha, sigma2, bs_beta, s2_beta, p_id, k_id));
    return __result;
END_RCPP
}
// sample_delta_scalar
double sample_delta_scalar(NumericVector y, List x, NumericMatrix groups, NumericMatrix beta, List gamma, NumericVector delta, NumericMatrix z, double alpha, double sigma2, List bs_beta, double s2_delta, int k_id);
RcppExport SEXP FREE_sample_delta_scalar(SEXP ySEXP, SEXP xSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP deltaSEXP, SEXP zSEXP, SEXP alphaSEXP, SEXP sigma2SEXP, SEXP bs_betaSEXP, SEXP s2_deltaSEXP, SEXP k_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< List >::type bs_beta(bs_betaSEXP);
    Rcpp::traits::input_parameter< double >::type s2_delta(s2_deltaSEXP);
    Rcpp::traits::input_parameter< int >::type k_id(k_idSEXP);
    __result = Rcpp::wrap(sample_delta_scalar(y, x, groups, beta, gamma, delta, z, alpha, sigma2, bs_beta, s2_delta, k_id));
    return __result;
END_RCPP
}
// sample_gamma_scalar
double sample_gamma_scalar(NumericVector y, List x, NumericMatrix groups, NumericMatrix beta, List gamma, NumericVector delta, NumericMatrix z, double alpha, double sigma2, List bs_beta, NumericVector sigma2_gamma, int q_id, int G_id);
RcppExport SEXP FREE_sample_gamma_scalar(SEXP ySEXP, SEXP xSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP deltaSEXP, SEXP zSEXP, SEXP alphaSEXP, SEXP sigma2SEXP, SEXP bs_betaSEXP, SEXP sigma2_gammaSEXP, SEXP q_idSEXP, SEXP G_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< List >::type bs_beta(bs_betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma2_gamma(sigma2_gammaSEXP);
    Rcpp::traits::input_parameter< int >::type q_id(q_idSEXP);
    Rcpp::traits::input_parameter< int >::type G_id(G_idSEXP);
    __result = Rcpp::wrap(sample_gamma_scalar(y, x, groups, beta, gamma, delta, z, alpha, sigma2, bs_beta, sigma2_gamma, q_id, G_id));
    return __result;
END_RCPP
}
// sample_sigma2_scalar
double sample_sigma2_scalar(NumericVector y, List x, NumericMatrix groups, NumericMatrix beta, List gamma, NumericVector delta, NumericMatrix z, double alpha, double sigma2, List bs_beta, double phi1, double psi1);
RcppExport SEXP FREE_sample_sigma2_scalar(SEXP ySEXP, SEXP xSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP deltaSEXP, SEXP zSEXP, SEXP alphaSEXP, SEXP sigma2SEXP, SEXP bs_betaSEXP, SEXP phi1SEXP, SEXP psi1SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< List >::type bs_beta(bs_betaSEXP);
    Rcpp::traits::input_parameter< double >::type phi1(phi1SEXP);
    Rcpp::traits::input_parameter< double >::type psi1(psi1SEXP);
    __result = Rcpp::wrap(sample_sigma2_scalar(y, x, groups, beta, gamma, delta, z, alpha, sigma2, bs_beta, phi1, psi1));
    return __result;
END_RCPP
}
// sample_sigma2_gamma_scalar
double sample_sigma2_gamma_scalar(List gamma, double phi2, double psi2, int q_id);
RcppExport SEXP FREE_sample_sigma2_gamma_scalar(SEXP gammaSEXP, SEXP phi2SEXP, SEXP psi2SEXP, SEXP q_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type phi2(phi2SEXP);
    Rcpp::traits::input_parameter< double >::type psi2(psi2SEXP);
    Rcpp::traits::input_parameter< int >::type q_id(q_idSEXP);
    __result = Rcpp::wrap(sample_sigma2_gamma_scalar(gamma, phi2, psi2, q_id));
    return __result;
END_RCPP
}
// calc_bs_scalar
NumericMatrix calc_bs_scalar(NumericVector x, NumericVector knots, int degree, NumericVector boundary_knots);
RcppExport SEXP FREE_calc_bs_scalar(SEXP xSEXP, SEXP knotsSEXP, SEXP degreeSEXP, SEXP boundary_knotsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< int >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type boundary_knots(boundary_knotsSEXP);
    __result = Rcpp::wrap(calc_bs_scalar(x, knots, degree, boundary_knots));
    return __result;
END_RCPP
}
// lnL
double lnL(const mat& y, const mat& x, NumericMatrix groups, mat& beta, List gamma, double sigma2, double rho, const mat& b_splines_mat, double s2_beta, int n, IntegerVector n_j, int n_k, int n_q, List bin_id);
RcppExport SEXP FREE_lnL(SEXP ySEXP, SEXP xSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP sigma2SEXP, SEXP rhoSEXP, SEXP b_splines_matSEXP, SEXP s2_betaSEXP, SEXP nSEXP, SEXP n_jSEXP, SEXP n_kSEXP, SEXP n_qSEXP, SEXP bin_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const mat& >::type b_splines_mat(b_splines_matSEXP);
    Rcpp::traits::input_parameter< double >::type s2_beta(s2_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< int >::type n_k(n_kSEXP);
    Rcpp::traits::input_parameter< int >::type n_q(n_qSEXP);
    Rcpp::traits::input_parameter< List >::type bin_id(bin_idSEXP);
    __result = Rcpp::wrap(lnL(y, x, groups, beta, gamma, sigma2, rho, b_splines_mat, s2_beta, n, n_j, n_k, n_q, bin_id));
    return __result;
END_RCPP
}
// ln_beta_dens
double ln_beta_dens(double x0, IntegerVector ids, const mat& y, const mat& x, NumericMatrix groups, mat& beta, List gamma, double sigma2, double rho, const mat& b_splines_mat, double s2_beta, int n, IntegerVector n_j, int n_k, int n_q, NumericVector sigma2_gamma, List bin_id);
RcppExport SEXP FREE_ln_beta_dens(SEXP x0SEXP, SEXP idsSEXP, SEXP ySEXP, SEXP xSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP sigma2SEXP, SEXP rhoSEXP, SEXP b_splines_matSEXP, SEXP s2_betaSEXP, SEXP nSEXP, SEXP n_jSEXP, SEXP n_kSEXP, SEXP n_qSEXP, SEXP sigma2_gammaSEXP, SEXP bin_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< const mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const mat& >::type b_splines_mat(b_splines_matSEXP);
    Rcpp::traits::input_parameter< double >::type s2_beta(s2_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< int >::type n_k(n_kSEXP);
    Rcpp::traits::input_parameter< int >::type n_q(n_qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma2_gamma(sigma2_gammaSEXP);
    Rcpp::traits::input_parameter< List >::type bin_id(bin_idSEXP);
    __result = Rcpp::wrap(ln_beta_dens(x0, ids, y, x, groups, beta, gamma, sigma2, rho, b_splines_mat, s2_beta, n, n_j, n_k, n_q, sigma2_gamma, bin_id));
    return __result;
END_RCPP
}
// slice_sample_beta
double slice_sample_beta(double x0, IntegerVector ids, const mat& y, const mat& x, NumericMatrix groups, mat& beta, List gamma, double sigma2, double rho, const mat& b_splines_mat, double s2_beta, int n, IntegerVector n_j, int n_k, int n_q, NumericVector sigma2_gamma, List bin_id);
RcppExport SEXP FREE_slice_sample_beta(SEXP x0SEXP, SEXP idsSEXP, SEXP ySEXP, SEXP xSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP sigma2SEXP, SEXP rhoSEXP, SEXP b_splines_matSEXP, SEXP s2_betaSEXP, SEXP nSEXP, SEXP n_jSEXP, SEXP n_kSEXP, SEXP n_qSEXP, SEXP sigma2_gammaSEXP, SEXP bin_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< const mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const mat& >::type b_splines_mat(b_splines_matSEXP);
    Rcpp::traits::input_parameter< double >::type s2_beta(s2_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< int >::type n_k(n_kSEXP);
    Rcpp::traits::input_parameter< int >::type n_q(n_qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma2_gamma(sigma2_gammaSEXP);
    Rcpp::traits::input_parameter< List >::type bin_id(bin_idSEXP);
    __result = Rcpp::wrap(slice_sample_beta(x0, ids, y, x, groups, beta, gamma, sigma2, rho, b_splines_mat, s2_beta, n, n_j, n_k, n_q, sigma2_gamma, bin_id));
    return __result;
END_RCPP
}
// ln_gamma_dens
double ln_gamma_dens(double x0, IntegerVector ids, const mat& y, const mat& x, NumericMatrix groups, mat& beta, List gamma, double sigma2, double rho, const mat& b_splines_mat, double s2_beta, int n, IntegerVector n_j, int n_k, int n_q, NumericVector sigma2_gamma, List bin_id);
RcppExport SEXP FREE_ln_gamma_dens(SEXP x0SEXP, SEXP idsSEXP, SEXP ySEXP, SEXP xSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP sigma2SEXP, SEXP rhoSEXP, SEXP b_splines_matSEXP, SEXP s2_betaSEXP, SEXP nSEXP, SEXP n_jSEXP, SEXP n_kSEXP, SEXP n_qSEXP, SEXP sigma2_gammaSEXP, SEXP bin_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< const mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const mat& >::type b_splines_mat(b_splines_matSEXP);
    Rcpp::traits::input_parameter< double >::type s2_beta(s2_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< int >::type n_k(n_kSEXP);
    Rcpp::traits::input_parameter< int >::type n_q(n_qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma2_gamma(sigma2_gammaSEXP);
    Rcpp::traits::input_parameter< List >::type bin_id(bin_idSEXP);
    __result = Rcpp::wrap(ln_gamma_dens(x0, ids, y, x, groups, beta, gamma, sigma2, rho, b_splines_mat, s2_beta, n, n_j, n_k, n_q, sigma2_gamma, bin_id));
    return __result;
END_RCPP
}
// slice_sample_gamma
double slice_sample_gamma(double x0, IntegerVector ids, const mat& y, const mat& x, NumericMatrix groups, mat& beta, List gamma, double sigma2, double rho, const mat& b_splines_mat, double s2_beta, int n, IntegerVector n_j, int n_k, int n_q, NumericVector sigma2_gamma, List bin_id);
RcppExport SEXP FREE_slice_sample_gamma(SEXP x0SEXP, SEXP idsSEXP, SEXP ySEXP, SEXP xSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP sigma2SEXP, SEXP rhoSEXP, SEXP b_splines_matSEXP, SEXP s2_betaSEXP, SEXP nSEXP, SEXP n_jSEXP, SEXP n_kSEXP, SEXP n_qSEXP, SEXP sigma2_gammaSEXP, SEXP bin_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< const mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const mat& >::type b_splines_mat(b_splines_matSEXP);
    Rcpp::traits::input_parameter< double >::type s2_beta(s2_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< int >::type n_k(n_kSEXP);
    Rcpp::traits::input_parameter< int >::type n_q(n_qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma2_gamma(sigma2_gammaSEXP);
    Rcpp::traits::input_parameter< List >::type bin_id(bin_idSEXP);
    __result = Rcpp::wrap(slice_sample_gamma(x0, ids, y, x, groups, beta, gamma, sigma2, rho, b_splines_mat, s2_beta, n, n_j, n_k, n_q, sigma2_gamma, bin_id));
    return __result;
END_RCPP
}
// ln_rho_dens
double ln_rho_dens(double x0, IntegerVector ids, const mat& y, const mat& x, NumericMatrix groups, mat& beta, List gamma, double sigma2, double rho, const mat& b_splines_mat, double s2_beta, int n, IntegerVector n_j, int n_k, int n_q, NumericVector sigma2_gamma, List bin_id);
RcppExport SEXP FREE_ln_rho_dens(SEXP x0SEXP, SEXP idsSEXP, SEXP ySEXP, SEXP xSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP sigma2SEXP, SEXP rhoSEXP, SEXP b_splines_matSEXP, SEXP s2_betaSEXP, SEXP nSEXP, SEXP n_jSEXP, SEXP n_kSEXP, SEXP n_qSEXP, SEXP sigma2_gammaSEXP, SEXP bin_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< const mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const mat& >::type b_splines_mat(b_splines_matSEXP);
    Rcpp::traits::input_parameter< double >::type s2_beta(s2_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< int >::type n_k(n_kSEXP);
    Rcpp::traits::input_parameter< int >::type n_q(n_qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma2_gamma(sigma2_gammaSEXP);
    Rcpp::traits::input_parameter< List >::type bin_id(bin_idSEXP);
    __result = Rcpp::wrap(ln_rho_dens(x0, ids, y, x, groups, beta, gamma, sigma2, rho, b_splines_mat, s2_beta, n, n_j, n_k, n_q, sigma2_gamma, bin_id));
    return __result;
END_RCPP
}
// slice_sample_rho
double slice_sample_rho(double x0, IntegerVector ids, const mat& y, const mat& x, NumericMatrix groups, mat& beta, List gamma, double sigma2, double rho, const mat& b_splines_mat, double s2_beta, int n, IntegerVector n_j, int n_k, int n_q, NumericVector sigma2_gamma, List bin_id);
RcppExport SEXP FREE_slice_sample_rho(SEXP x0SEXP, SEXP idsSEXP, SEXP ySEXP, SEXP xSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP sigma2SEXP, SEXP rhoSEXP, SEXP b_splines_matSEXP, SEXP s2_betaSEXP, SEXP nSEXP, SEXP n_jSEXP, SEXP n_kSEXP, SEXP n_qSEXP, SEXP sigma2_gammaSEXP, SEXP bin_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< const mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const mat& >::type b_splines_mat(b_splines_matSEXP);
    Rcpp::traits::input_parameter< double >::type s2_beta(s2_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< int >::type n_k(n_kSEXP);
    Rcpp::traits::input_parameter< int >::type n_q(n_qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma2_gamma(sigma2_gammaSEXP);
    Rcpp::traits::input_parameter< List >::type bin_id(bin_idSEXP);
    __result = Rcpp::wrap(slice_sample_rho(x0, ids, y, x, groups, beta, gamma, sigma2, rho, b_splines_mat, s2_beta, n, n_j, n_k, n_q, sigma2_gamma, bin_id));
    return __result;
END_RCPP
}
// calc_bs
NumericMatrix calc_bs(NumericVector x, NumericVector knots, int degree, NumericVector boundary_knots);
RcppExport SEXP FREE_calc_bs(SEXP xSEXP, SEXP knotsSEXP, SEXP degreeSEXP, SEXP boundary_knotsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< int >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type boundary_knots(boundary_knotsSEXP);
    __result = Rcpp::wrap(calc_bs(x, knots, degree, boundary_knots));
    return __result;
END_RCPP
}
// sample_sigma2
double sample_sigma2(const mat& y, const mat& x, NumericMatrix groups, mat& beta, List gamma, double sigma2, double rho, const mat& b_splines_mat, double s2h_a, double s2h_b, int n, IntegerVector n_j, int n_k, int n_q, List bin_id);
RcppExport SEXP FREE_sample_sigma2(SEXP ySEXP, SEXP xSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP sigma2SEXP, SEXP rhoSEXP, SEXP b_splines_matSEXP, SEXP s2h_aSEXP, SEXP s2h_bSEXP, SEXP nSEXP, SEXP n_jSEXP, SEXP n_kSEXP, SEXP n_qSEXP, SEXP bin_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const mat& >::type b_splines_mat(b_splines_matSEXP);
    Rcpp::traits::input_parameter< double >::type s2h_a(s2h_aSEXP);
    Rcpp::traits::input_parameter< double >::type s2h_b(s2h_bSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< int >::type n_k(n_kSEXP);
    Rcpp::traits::input_parameter< int >::type n_q(n_qSEXP);
    Rcpp::traits::input_parameter< List >::type bin_id(bin_idSEXP);
    __result = Rcpp::wrap(sample_sigma2(y, x, groups, beta, gamma, sigma2, rho, b_splines_mat, s2h_a, s2h_b, n, n_j, n_k, n_q, bin_id));
    return __result;
END_RCPP
}
// sample_sigma2_gamma
double sample_sigma2_gamma(List gamma, double s2hg_a, double s2hg_b, int n_q, int q_id);
RcppExport SEXP FREE_sample_sigma2_gamma(SEXP gammaSEXP, SEXP s2hg_aSEXP, SEXP s2hg_bSEXP, SEXP n_qSEXP, SEXP q_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type s2hg_a(s2hg_aSEXP);
    Rcpp::traits::input_parameter< double >::type s2hg_b(s2hg_bSEXP);
    Rcpp::traits::input_parameter< int >::type n_q(n_qSEXP);
    Rcpp::traits::input_parameter< int >::type q_id(q_idSEXP);
    __result = Rcpp::wrap(sample_sigma2_gamma(gamma, s2hg_a, s2hg_b, n_q, q_id));
    return __result;
END_RCPP
}
