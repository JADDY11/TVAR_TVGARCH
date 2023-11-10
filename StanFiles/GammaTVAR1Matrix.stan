// The .Stan file for the Gamma TV-AR(1) process with multivariate Horseshoe priors
data {
  int<lower=0> N;
  int<lower=0> KTVAR;
  vector[KTVAR] vec_null;
  vector[N] y;
  vector[N] ytm1;
  matrix[N, 1] modmat;
  matrix[N, KTVAR] B_TVAR;
}

parameters {
  vector[1] beta;
  vector[KTVAR] beta_TVAR;
  vector[1] tau0;
  vector[KTVAR] tau;
  vector<lower=0>[KTVAR] phi_TVAR;
  vector<lower=0>[KTVAR] phi_tau_TVAR;
  real<lower=0> lambda_TVAR;
  real<lower=0> lambda_tau_TVAR;
  cholesky_factor_corr[KTVAR] L_beta_TVAR;
  cholesky_factor_corr[KTVAR] L_tau_TVAR;
  real<lower=0> eta_beta_TVAR;
  real<lower=0> eta_tau_TVAR;
}

transformed parameters{
  vector[N] yhat1;
  vector[N] yhat2;
  vector[N] yhat;
  vector[N] vari;
  
  vector[N] gam_alpha;
  vector[N] gam_beta;
  
  cov_matrix[KTVAR] K_est_beta_TVAR;
  cov_matrix[KTVAR] K_est_tau_TVAR;
  
  yhat1 = to_vector(modmat*beta);
  yhat2 =  log(to_vector(ytm1)) .* to_vector(B_TVAR * beta_TVAR);
  yhat = exp(yhat1 + yhat2);
  
  vari = exp(to_vector(modmat*tau0) + to_vector(B_TVAR * tau));
  
  gam_alpha = (to_vector(yhat).*to_vector(yhat))./vari;
  gam_beta = to_vector(yhat)./vari;
  
  K_est_beta_TVAR = phi_TVAR * to_row_vector(phi_TVAR) .* (L_beta_TVAR * L_beta_TVAR');
  K_est_tau_TVAR = phi_tau_TVAR * to_row_vector(phi_tau_TVAR) .* (L_tau_TVAR * L_tau_TVAR');
}

model {
  beta_TVAR ~ multi_normal(vec_null, K_est_beta_TVAR*lambda_TVAR^2);
  tau ~ multi_normal(vec_null, K_est_tau_TVAR*lambda_tau_TVAR^2);
  phi_TVAR ~ cauchy(0, 1);
  phi_tau_TVAR ~ cauchy(0, 1);  
  lambda_TVAR ~ cauchy(0, 1);
  lambda_tau_TVAR ~ cauchy(0, 1);
  L_beta_TVAR ~ lkj_corr_cholesky(eta_beta_TVAR);
  L_tau_TVAR ~ lkj_corr_cholesky(eta_tau_TVAR);
  eta_beta_TVAR ~ cauchy(0, 1);
  eta_tau_TVAR ~ cauchy(0, 1);
  
  y ~ gamma(gam_alpha, gam_beta);
}

generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;
  for (n in 1:N) {
    y_rep[n] = gamma_rng(gam_alpha[n], gam_beta[n]);
    log_lik[n] = gamma_lpdf(y[n]| gam_alpha[n], gam_beta[n]);
  }
}
