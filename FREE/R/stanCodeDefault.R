stanCodeDefault <-
function(){
  # Define a default stan code for stanFit analysis
  stan.code <- "data{
    int<lower=0> N;    // number of observations
    int<lower=0> D;    // number of bins
    int<lower=0> k;    // number of preds

    int Kt;            // number of basis functions

    vector[D] Y[N];    // response matrix
    vector[k] X[N];    // predictor matrix
    matrix[D,Kt] BS;   // B-splines pre-evaluated

    cov_matrix[Kt] CovMat;    // prior precision for spline effects
}

transformed data{
    vector[Kt] mu_beta;    // prior mean for splines

    for(i in 1:Kt){
        mu_beta[i] <- 0;
    }
}

parameters{
    matrix[k,Kt] beta;    // matrix of fixed effect spline coefficients
    real<lower=0,upper=10> y_sig;  // regression variance
    vector<lower=0,upper=10>[k] beta_sig;    // tuning parameter
}

transformed parameters{
    vector<lower=0>[k] beta_tau2;

    for(kcur in 1:k){
        beta_tau2[kcur] <- pow(beta_sig[kcur], -2);
    }
}

model{
    for(kcur in 1:k){
        (beta[kcur])' ~ multi_normal_prec(mu_beta, beta_tau2[kcur] * CovMat);
    }

    // outcome likelihood
    for(n in 1:N){
        Y[n] ~ normal((BS*beta')*X[n], y_sig);
    }
}"
  return(stan.code)
}
