// The .Stan file for the Gamma TV-AR(1) process with Horseshoe priors
data {
  int<lower=0> N;
  int<lower=0> KTVAR;
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
  beta_TVAR ~ normal(0, (phi_TVAR[KTVAR])*(lambda_TVAR));
  phi_TVAR ~ cauchy(0, 1);
  lambda_TVAR ~ cauchy(0, 1);
  
  tau ~ normal(0, (phi_tau_TVAR[KTVAR])*(lambda_tau_TVAR));
  phi_tau_TVAR ~ cauchy(0, 1);
  lambda_tau_TVAR ~ cauchy(0, 1);
  
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
