// The .Stan file for the Gamma TV-GARCH(1) process with multivariate Horseshoe priors
data {
  int<lower=0> N;
  int<lower=0> KARCH;
  int<lower=0> KGARCH;
  vector[KARCH] vec_null;
  vector[N] y;
  vector[N] pred;
  real resi2[N];
  matrix[N, KARCH] B_ARCH;
  matrix[N, KARCH] B_GARCH;
}

parameters {
  real tau_0;
  vector[KARCH] tau_ARCH;
  vector[KGARCH] tau_GARCH;
  vector<lower=0>[KARCH] phi_ARCH;
  vector<lower=0>[KGARCH] phi_GARCH;
  real<lower=0> lambda_ARCH;
  real<lower=0> lambda_GARCH;
  cholesky_factor_corr[KARCH] L_ARCH;
  cholesky_factor_corr[KGARCH] L_GARCH;
  real<lower=0> eta_ARCH;
  real<lower=0> eta_GARCH;
}

transformed parameters{
  real<lower = 0> varhat[N];
  vector[N] pred_ARCH;
  vector[N] pred_GARCH;
  cov_matrix[KARCH] K_est_ARCH;
  cov_matrix[KGARCH] K_est_GARCH;
  
  pred_ARCH = to_vector(B_ARCH * tau_ARCH);
  pred_GARCH = to_vector(B_GARCH * tau_GARCH);

  varhat[1] = exp(tau_0);
  for(n in 2:N)
    varhat[n] = exp(tau_0
    + pred_ARCH[n]*(pow(resi2[n - 1]/sqrt(varhat[n - 1]), 2))
    + pred_GARCH[n]*log(varhat[n - 1]));
    
  K_est_ARCH = phi_ARCH * to_row_vector(phi_ARCH) .* (L_ARCH * L_ARCH');
  K_est_GARCH = phi_GARCH * to_row_vector(phi_GARCH) .* (L_GARCH * L_GARCH');
}

model {
  tau_ARCH ~ multi_normal(vec_null, K_est_ARCH*lambda_ARCH^2);
  tau_GARCH ~ multi_normal(vec_null, K_est_GARCH*lambda_GARCH^2);
  phi_ARCH ~ cauchy(0, 1);
  phi_GARCH ~ cauchy(0, 1);
  lambda_ARCH ~ cauchy(0, 1);
  lambda_GARCH ~ cauchy(0, 1);
  L_ARCH ~ lkj_corr_cholesky(eta_ARCH);
  L_GARCH ~ lkj_corr_cholesky(eta_GARCH);
  eta_ARCH ~ cauchy(0, 1);
  eta_GARCH ~ cauchy(0, 1);
  
  for(n in 2:N)
    y[n] ~ gamma(pow(pred[n], 2)/varhat[n], pred[n]/varhat[n]);
}

generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;
  for (n in 1:N) {
    y_rep[n] = gamma_rng(pow(pred[n], 2)/varhat[n], pred[n]/varhat[n]);
    log_lik[n] = gamma_lpdf(y[n]| pow(pred[n], 2)/varhat[n], pred[n]/varhat[n]);
  }
}
