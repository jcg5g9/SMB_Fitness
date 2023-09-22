// saved as vbgf.stan
data {
  int<lower=0> Nobs;                    // number of observations 
  int<lower=0> Nages;                    // number of observations 
  real length[Nobs];                    // length
  real age[Nobs];                       // length
  
  int<lower=0> Nind;                    // number of individuals
  int<lower=0> Ncoef;                   // Number of predictors
  matrix[Nind, Ncoef] X;                // Predictor matrix
  real q[Nind];                         // Proportion SMB
  
  int<lower=1, upper=Nind> id[Nobs];    // Sample ID
}
parameters {
  // VBGF Params ----
  vector[3] mu;                     // VBGF params
  matrix[2, 3] eta_group;           // Group level deviation
  matrix[Nind, 3] eta_ind;          // Individual deviation from VBGF parameters
  
  // Regressors
  vector[Ncoef] beta_linf;
  vector[Ncoef] beta_k;
  vector[Ncoef] beta_t0;
  
  // Variance ----
  corr_matrix[3] Omega_ind;          // Prior correlation for individual-level variation
  vector<lower=0>[3] tau_ind;        // Prior scale for individual-level variation
  
  corr_matrix[3] Omega_group;        // Prior correlation for group-level variation
  vector<lower=0>[3] tau_group;      // Prior scale for group-level variation
  
  real<lower=0> sigma;               // prediction error scale
}
transformed parameters {
  // Predicted length
  vector[Nobs] length_hat;
  
  // Global parameters
  real mu_linf = exp(mu[1]);
  real mu_k = exp(mu[2]);
  real mu_t0 = mu[3];
  
  // Group level parameters
  vector[2] linf_group = exp(mu[1] + eta_group[,1]);
  vector[2] k_group = exp(mu[2] + eta_group[,2]);
  vector[2] t0_group = mu[3] + eta_group[,3];
  
  // Individual level parameters
  vector[Nind] linf_ind = exp(X * beta_linf); 
  vector[Nind] k_ind = exp(X * beta_k);
  vector[Nind] t0_ind = X * beta_t0;
  
  for(i in 1:Nind){
    linf_ind[i] = linf_ind[i] * exp(mu[1] + eta_group[1,1] * q[i] + eta_group[2,1] * (1-q[i]) + eta_ind[i,1]);
    k_ind[i] = k_ind[i] * exp(mu[2] + eta_group[1,1] * q[i] + eta_group[2,2] * (1-q[i]) + eta_ind[i,2]);
    t0_ind[i] = t0_ind[i] + mu[3] + eta_group[1,1] * q[i] + eta_group[2,3] * (1-q[i]) + eta_ind[i,3];
  }
  
  // Predicted length
  for(i in 1:Nobs){
    length_hat[i] = linf_ind[id[i]] * (1-exp(-k_ind[id[i]]* (age[i] - (t0_ind[id[i]]))));
  }
}
model {
  vector[3] mu_mvn;
  for(i in 1:3){
    mu_mvn[i] = 0;
  }
  
  // Group and individual variation priors
  tau_ind ~ cauchy(0, 2.5);
  Omega_ind ~ lkj_corr(2);
  
  tau_group ~ cauchy(0, 2.5);
  Omega_group ~ lkj_corr(2);
  
  for(i in 1:Nind){
    eta_ind[i,] ~ multi_normal(mu_mvn, quad_form_diag(Omega_ind, tau_ind));
  }
  for(i in 1:2){
    eta_group[i,] ~ multi_normal(mu_mvn, quad_form_diag(Omega_group, tau_group));
  }
  
  // Regessor priors
  beta_linf ~ normal(0, 1);
  beta_k ~ normal(0, 1);
  beta_t0 ~ normal(0, 1);
  
  // Likelihood
  length ~ normal(length_hat, sigma);
} 
generated quantities{
  // Predicted length
  matrix[Nages, 7] length_pred;
  
  for(i in 1:Nages){
    length_pred[i,1] = mu_linf * (1 - exp(-mu_k * (i - mu_t0))); # Global
    length_pred[i,2] = linf_group[1] * (1 - exp(-k_group[1] * (i - t0_group[1]))); # SMB - females
    length_pred[i,3] = linf_group[2] * (1 - exp(-k_group[2] * (i - t0_group[2]))); # Neosho - females
    length_pred[i,4] = linf_group[1] * exp(beta_linf[1]) * (1 - exp(-k_group[1] * exp(beta_k[1]) * (i - t0_group[1] - beta_t0[1]))); # SMB - males
    length_pred[i,5] = linf_group[2] * exp(beta_linf[1]) * (1 - exp(-k_group[2] * exp(beta_k[1]) * (i - t0_group[2] - beta_t0[1]))); # Neosho - males
  }
}
