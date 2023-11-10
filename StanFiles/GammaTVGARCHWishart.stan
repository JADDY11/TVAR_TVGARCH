// The .Stan file for the Gamma TV-GARCH(1) process with a Inverse-Wishard prior
data {
  int<lower=0> N;
  int<lower=0> KARCH;
  int<lower=0> KGARCH;
  vector[KARCH] vec_null;
  vector[N] y;
  vector[N] pred;
  real resi2[N];
  matrix[KARCH, KARCH] Sig_Idnt;
  matrix[N, KARCH] B_ARCH;
  matrix[N, KARCH] B_GARCH;
}

parameters {
  real tau_0;
  vector[KARCH] tau_ARCH;
  vector[KGARCH] tau_GARCH;
  cov_matrix[KARCH] K_est_ARCH;
  cov_matrix[KARCH] K_est_GARCH;
}

transformed parameters{
  real<lower = 0> varhat[N];
  vector[N] pred_ARCH;
  vector[N] pred_GARCH;
  
  pred_ARCH = to_vector(B_ARCH * tau_ARCH);
  pred_GARCH = to_vector(B_GARCH * tau_GARCH);

  varhat[1] = exp(tau_0);
  for(n in 2:N)
    varhat[n] = exp(tau_0
    + pred_ARCH[n]*(pow(resi2[n - 1]/sqrt(varhat[n - 1]), 2))
    + pred_GARCH[n]*log(varhat[n - 1]));
}

model {
  tau_ARCH ~ multi_normal(vec_null, K_est_ARCH);
  tau_GARCH ~ multi_normal(vec_null, K_est_GARCH);
  K_est_ARCH ~ inv_wishart(20, Sig_Idnt);
  K_est_GARCH ~ inv_wishart(20, Sig_Idnt);
  
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
