//  VBGF with sex/river effects saved as vbgf4.stan
//  Includes ancestry level random effects
data {
  // Length-at-age data
  int<lower=0> Nobs;                    // number of observations 
  int<lower=0> Nages;                   // number of observations 
  real length[Nobs];                    // length
  real age[Nobs];                       // length
  
  // Individual level data
  int<lower=0> Nind;                    // number of individuals
  int<lower=0> Ncoef;                   // Number of predictors
  matrix[Nind, Ncoef] X;                // Design matrix
  matrix[4,Ncoef] Xhat;                     // Prediction matrix
  vector[Nind] q;                       // Proportion SMB
  
  int<lower=1, upper=Nind> id[Nobs];    // Sample ID
  vector[3] Zero;                       // Vector of zero for mean
}
parameters {
  // VBGF Params ----
  real<lower=0> mu_linf;                // asymptotic length
  real<lower=0> mu_k;                   // growth coef
  real mu_t0;                           // age at length 0
  matrix[2, 3] eta_lineage;             // Ancestry level deviation
  matrix[Nind, 3] eta_ind;              // Individual deviation from VBGF parameters
  
  // Regressors
  vector[Ncoef] beta_linf;
  vector[Ncoef] beta_k;
  vector[Ncoef] beta_t0;
  
  // Variance ----
  cholesky_factor_corr[3] Lcorr_ind;          // Prior correlation for individual-level variation
  vector<lower=0>[3] sigma_ind;        // Prior scale for individual-level variation
  
  cholesky_factor_corr[3] Lcorr_group;           // Prior correlation for group-level variation
  vector<lower=0>[3] sigma_group;         // Prior scale for group-level variation
  
  real<lower=0> sigma;                  // observation error
}
transformed parameters {
  // Predicted length
  vector[Nobs] length_hat;
  
  // Group level parameters
  vector[2] linf_lineage = mu_linf * exp(eta_lineage[,1]);
  vector[2] k_lineage = mu_k * exp(eta_lineage[,2]);
  vector[2] t0_lineage = mu_t0 + eta_lineage[,3];
  
  real eta_smb_linf = eta_lineage[1,1];
  real eta_smb_k = eta_lineage[1,2];
  real eta_smb_t0 = eta_lineage[1,3];
  
  real eta_n_linf = eta_lineage[2,1];
  real eta_n_k = eta_lineage[2,2];
  real eta_n_t0 = eta_lineage[2,3];
  
  // Individual level parameters
  vector[Nind] linf_ind = mu_linf * exp(X * beta_linf + eta_smb_linf * q + eta_n_linf * (1-q)); 
  vector[Nind] k_ind = mu_k * exp(X * beta_k + eta_smb_k * q + eta_n_k * (1-q));
  vector[Nind] t0_ind = mu_t0 + X * beta_t0 + eta_smb_t0 * q + eta_n_t0 * (1-q);
  
  // Predicted length
  for(i in 1:Nobs){
    length_hat[i] = linf_ind[id[i]] * exp(eta_ind[id[i],1]) * (1-exp(-k_ind[id[i]] * exp(eta_ind[id[i],2]) * (age[i] - (t0_ind[id[i]] + eta_ind[id[i],3]))));
  }
}
model {
  // Priors
  // - Global parameters
  // - Starks, T. A., & Rodger, A. W. (2020). Otolith and scale‐based growth standards for lotic Smallmouth Bass. North American Journal of Fisheries Management, 40(4), 986-994.
  mu_linf ~ lognormal(log(578), 0.02432599);
  mu_k ~ lognormal(log(0.125), 0.04580929);
  mu_t0 ~ normal(-1.79, 0.0625);
  
  // - Ancestry level variation priors
  sigma_group ~ cauchy(0, 0.5);
  Lcorr_group ~ lkj_corr_cholesky(3); // Centered around 0 https://mjskay.github.io/ggdist/reference/lkjcorr_marginal.html

  for(i in 1:2){
    eta_lineage[i,] ~ multi_normal_cholesky(Zero, diag_pre_multiply(sigma_group, Lcorr_group));
  }
  
  // - Individual variation priors
  sigma_ind ~ cauchy(0, 0.5);
  Lcorr_ind ~ lkj_corr_cholesky(3); // Centered around 0 https://mjskay.github.io/ggdist/reference/lkjcorr_marginal.html

  for(i in 1:Nind){
    eta_ind[i,] ~ multi_normal_cholesky(Zero, diag_pre_multiply(sigma_ind, Lcorr_ind));
  }

  // - Regessor priors
  beta_linf ~ normal(0, 1);
  beta_k ~ normal(0, 1);
  beta_t0 ~ normal(0, 1);
  
  
  // Likelihood
  length ~ normal(length_hat, sigma);
} 
generated quantities{
  // Predicted length
  matrix[9, Nages] length_pred;
  
  // Correlation matrices
  matrix[3, 3] Omega_group;
  matrix[3, 3] Omega_ind;
  
  
  // Parameters for prediction
  vector[4] linf_pred = exp(Xhat * beta_linf) * mu_linf; 
  vector[4] k_pred = exp(Xhat * beta_k) * mu_k;
  vector[4] t0_pred = Xhat * beta_t0 + mu_t0;
  
  
  for(i in 1:Nages){
    // Global
    length_pred[1,i] = mu_linf * (1 - exp(-mu_k * (i - mu_t0))); // Global
    
    // SMB
    length_pred[2,i] = linf_pred[1] * exp(eta_smb_linf) * (1 - exp(-k_pred[1] * exp(eta_smb_k) * (i - t0_pred[1] - eta_smb_t0))); // Females river 1
    length_pred[3,i] = linf_pred[2] * exp(eta_smb_linf) * (1 - exp(-k_pred[2] * exp(eta_smb_k) * (i - t0_pred[2] - eta_smb_t0))); // Male river 1
    length_pred[4,i] = linf_pred[3] * exp(eta_smb_linf) * (1 - exp(-k_pred[3] * exp(eta_smb_k) * (i - t0_pred[3] - eta_smb_t0))); // Females river 2
    length_pred[5,i] = linf_pred[4] * exp(eta_smb_linf) * (1 - exp(-k_pred[4] * exp(eta_smb_k) * (i - t0_pred[4] - eta_smb_t0))); // Males river 2
    
    // Neosho
    length_pred[6,i] = linf_pred[1] * exp(eta_n_linf) * (1 - exp(-k_pred[1] * exp(eta_n_k) * (i - t0_pred[1] - eta_n_t0))); // Females river 1
    length_pred[7,i] = linf_pred[2] * exp(eta_n_linf) * (1 - exp(-k_pred[2] * exp(eta_n_k) * (i - t0_pred[2] - eta_n_t0))); // Male river 1
    length_pred[8,i] = linf_pred[3] * exp(eta_n_linf) * (1 - exp(-k_pred[3] * exp(eta_n_k) * (i - t0_pred[3] - eta_n_t0))); // Females river 2
    length_pred[9,i] = linf_pred[4] * exp(eta_n_linf) * (1 - exp(-k_pred[4] * exp(eta_n_k) * (i - t0_pred[4] - eta_n_t0))); // Males river 2
  }
  
  // Correlation matrices
  Omega_group = Lcorr_group * Lcorr_group ; // Correlation matrix of ancestry distribution
  Omega_ind = Lcorr_ind * Lcorr_ind ; // Correlation matrix of ind re distribution
}
