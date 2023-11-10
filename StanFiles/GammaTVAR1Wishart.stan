// The .Stan file for the Gamma TV-AR(1) process with a Inverse-Wishard prior
data {
  int<lower=0> N;
  int<lower=0> KTVAR;
  vector[KTVAR] vec_null;
  vector[N] y;
  vector[N] ytm1;
  matrix[KTVAR, KTVAR] Sig_Idnt;
  matrix[N, 1] modmat;
  matrix[N, KTVAR] B_TVAR;
}

parameters {
  vector[1] beta;
  vector[KTVAR] beta_TVAR;
  vector[1] tau0;
  vector[KTVAR] tau;
  cov_matrix[KTVAR] K_est_beta_TVAR;
  cov_matrix[KTVAR] K_est_tau_TVAR;
}

transformed parameters{
  vector[N] yhat1;
  vector[N] yhat2;
  vector[N] yhat;
  vector[N] vari;
  
  vector[N] gam_alpha;
  vector[N] gam_beta;
  
  yhat1 = to_vector(modmat*beta);
  yhat2 =  log(to_vector(ytm1)) .* to_vector(B_TVAR * beta_TVAR);
  yhat = exp(yhat1 + yhat2);
  
  vari = exp(to_vector(modmat*tau0) + to_vector(B_TVAR * tau));
  
  gam_alpha = (to_vector(yhat).*to_vector(yhat))./vari;
  gam_beta = to_vector(yhat)./vari;
}

model {
  beta_TVAR ~ multi_normal(vec_null, K_est_beta_TVAR);
  K_est_beta_TVAR ~ inv_wishart(20, Sig_Idnt);
  tau ~ multi_normal(vec_null, K_est_tau_TVAR);
  K_est_tau_TVAR ~ inv_wishart(20, Sig_Idnt);
  
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
