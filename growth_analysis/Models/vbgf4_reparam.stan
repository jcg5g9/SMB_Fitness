//  VBGF with sex/river effects saved as vbgf4.stan
//  Includes ancestry level random effects
//  Reparameterizes the VBGF to reduce parameter correlation
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
  matrix[4,2] Xhat;                     // Prediction matrix
  vector[Nind] q;                       // Proportion SMB
  
  real amax;                            // max age
  real amin;                            // min age
  
  int<lower=1, upper=Nind> id[Nobs];    // Sample ID
  vector[3] Zero;                       // Vector of zero for mean
}
parameters {
  // VBGF Params ----
  real mu_lmin;                         // estimated length at min age
  real mu_lmax;                         // estimated length at max age
  real mu_k;                            // growth coef
  vector[2] eta_lineage;             // Ancestry level deviation
  vector[Nind] eta_ind;              // Individual deviation from VBGF parameters
  
  // Regressors
  //vector[Ncoef] beta_lmin;
  //vector[Ncoef] beta_lmax;
  vector[Ncoef] beta_k;
  
  // Variance ----
  //cholesky_factor_corr[3] Lcorr_ind;          // Prior correlation for individual-level variation
  vector<lower=0>[3] sigma_ind;        // Prior scale for individual-level variation
  
  //cholesky_factor_corr[3] Lcorr_group;           // Prior correlation for group-level variation
  vector<lower=0>[3] sigma_group;         // Prior scale for group-level variation
  
  real<lower=0> sigma;                  // observation error
}
transformed parameters {
  // Predicted length
  vector[Nobs] length_hat;
  
  // Global parameters (natural scale)
  real mu_linf = mu_lmin + (mu_lmax - mu_lmin)/(1-exp(-mu_k * (amax - amin)));
  real mu_t0 = log(mu_linf/(mu_linf-mu_lmin))/(-mu_k) + amin;
  
  // Group level parameters (transformed scale)
  //vector[2] lmin_lineage = mu_lmin * exp(eta_lineage[,1]);
  //vector[2] lmax_lineage = mu_lmax * exp(eta_lineage[,2]);
  vector[2] k_lineage = mu_k * exp(eta_lineage);
  
  // Group level parameters (natural scale)
  vector[2] linf_lineage; 
  vector[2] t0_lineage; 
  
  // Individual levels parameters (transformed scale)
  //real eta_smb_lmin = eta_lineage[1,1];
  //real eta_smb_lmax = eta_lineage[1,2];
  real eta_smb_k = eta_lineage[1];
  
  //real eta_n_lmin = eta_lineage[2,1];
  //real eta_n_lmax = eta_lineage[2,2];
  real eta_n_k = eta_lineage[2];
  
  // Individual level (transformed scale)
  // vector[Nind] lmin_ind = mu_lmin * exp(X * beta_lmin + eta_smb_lmin * q + eta_n_lmin * (1-q) + eta_ind[id, 1]); 
  // vector[Nind] lmax_ind = mu_lmax * exp(X * beta_lmax + eta_smb_lmax * q + eta_n_lmax * (1-q) + eta_ind[id, 2]);
  vector[Nind] k_ind;
  
  // Individual level (natural scale)
  vector[Nind] linf_ind;
  vector[Nind] t0_ind;
  
  
  // Put parameters on natural scale
  // - Group
  for(i in 1:2){
    linf_lineage[i] = mu_lmin + (mu_lmax - mu_lmin)/(1-exp(-k_lineage[i] * (amax - amin)));
    t0_lineage[i] = log(mu_lmax/(mu_lmax-mu_lmin))/(-k_lineage[i]) + amin;
  }
  
  // - Individual
  for(i in 1:Nind){
    k_ind[i] = mu_k * exp(X[i,] * beta_k + eta_smb_k * q[i] + eta_n_k * (1-q[i]) +  + eta_ind[id[i]]);
    linf_ind[i] = mu_lmin + (mu_lmax - mu_lmin)/(1-exp(-k_ind[i] * (amax - amin)));
    t0_ind[i] = log(mu_lmax/(mu_lmax-mu_lmin))/(-k_ind[i]) + amin;
  }
  
  
  // Predict length
  for(i in 1:Nobs){
    length_hat[i] = linf_ind[id[i]] * (1-exp(-k_ind[id[i]] * (age[i] - (t0_ind[id[i]]))));
  }
}
model {
  // Priors
  // - Global parameters
  // - Starks, T. A., & Rodger, A. W. (2020). Otolith and scale‐based growth standards for lotic Smallmouth Bass. North American Journal of Fisheries Management, 40(4), 986-994.
  mu_lmax ~ lognormal(log(578), 0.02432599);
  mu_lmin ~ uniform(0, mu_lmax);
  mu_k ~ lognormal(log(0.125), 0.04580929);
  
  // - Ancestry level variation priors
  //sigma_group ~ cauchy(0, 0.5);
  //Lcorr_group ~ lkj_corr_cholesky(3); // Centered around 0 https://mjskay.github.io/ggdist/reference/lkjcorr_marginal.html
  
  for(i in 1:2){
    eta_lineage[i] ~ normal(0, sigma_group);// diag_pre_multiply(sigma_group, Lcorr_group));
  }
  
  // - Individual variation priors
  //sigma_ind ~ cauchy(0, 0.5);
  //Lcorr_ind ~ lkj_corr_cholesky(3); // Centered around 0 https://mjskay.github.io/ggdist/reference/lkjcorr_marginal.html
  
  for(i in 1:Nind){
    eta_ind[i] ~ normal(0, sigma_ind); //multi_normal_cholesky(Zero, diag_pre_multiply(sigma_ind, Lcorr_ind));
  }
  
  // - Regessor priors
  //beta_lmin ~ normal(0, 1);
  //beta_lmax ~ normal(0, 1);
  beta_k ~ normal(0, 1);
  
  
  // Likelihood
  length ~ normal(length_hat, sigma);
} 
generated quantities{
  // Predicted length
  matrix[9, Nages] length_pred;
  
  // Correlation matrices
  //matrix[3, 3] Omega_group;
  //matrix[3, 3] Omega_ind;
  
  
  // Parameters for prediction
  // - SMB
  //vector[4] lmin_pred_smb = mu_lmin; //exp(Xhat * beta_lmin) * mu_lmin * exp(eta_smb_lmin); 
  //vector[4] lmax_pred_smb = mu_lmax; //exp(Xhat * beta_lmax) * mu_lmax * exp(eta_smb_lmax); 
  vector[4] k_pred_smb = exp(Xhat * beta_k) * mu_k * exp(eta_smb_k);
  
  // -- Transform to natural scale
  vector[4] linf_pred_smb;
  vector[4] t0_pred_smb;
  
  // - Neosho
  //vector[4] lmin_pred_n = mu_lmin; // exp(Xhat * beta_lmin) * mu_lmin * exp(eta_n_lmin); 
  //vector[4] lmax_pred_n = mu_lmax; //exp(Xhat * beta_lmax) * mu_lmax * exp(eta_n_lmax); 
  vector[4] k_pred_n = exp(Xhat * beta_k) * mu_k * exp(eta_n_k);
  
  // -- Transform to natural scale
  vector[4] linf_pred_n;
  vector[4] t0_pred_n;
  
  
  // Transform
  for(i in 1:4){
    linf_pred_smb[i] = mu_lmin + (mu_lmax - mu_lmin)/(1-exp(-k_pred_smb[i] * (amax - amin)));
    t0_pred_smb[i] = log(mu_lmax/(mu_lmax-mu_lmin))/(-k_pred_smb[i]) + amin;
    
    linf_pred_n[i] = mu_lmin + (mu_lmax - mu_lmin)/(1-exp(-k_pred_n[i] * (amax - amin)));
    t0_pred_n[i] = log(mu_lmax/(mu_lmax-mu_lmin))/(-k_pred_n[i]) + amin;
  }
  
  
  for(i in 1:Nages){
    // Global
    length_pred[1,i] = mu_linf * (1 - exp(-mu_k * (i - mu_t0))); // Global
    
    // SMB
    length_pred[2,i] = linf_pred_smb[1] * (1 - exp(-k_pred_smb[1] * (i - t0_pred_smb[1]))); // Females river 1
    length_pred[3,i] = linf_pred_smb[2] * (1 - exp(-k_pred_smb[2] * (i - t0_pred_smb[2]))); // Male river 1
    length_pred[4,i] = linf_pred_smb[3] * (1 - exp(-k_pred_smb[3] * (i - t0_pred_smb[3]))); // Females river 2
    length_pred[5,i] = linf_pred_smb[4] * (1 - exp(-k_pred_smb[4] * (i - t0_pred_smb[4]))); // Males river 2
    
    // Neosho
    length_pred[6,i] = linf_pred_n[1] * (1 - exp(-k_pred_n[1] * (i - t0_pred_n[1]))); // Females river 1
    length_pred[7,i] = linf_pred_n[2] * (1 - exp(-k_pred_n[2] * (i - t0_pred_n[2]))); // Male river 1
    length_pred[8,i] = linf_pred_n[3] * (1 - exp(-k_pred_n[3] * (i - t0_pred_n[3]))); // Females river 2
    length_pred[9,i] = linf_pred_n[4] * (1 - exp(-k_pred_n[4] * (i - t0_pred_n[4]))); // Males river 2
  }
  
  // Correlation matrices
  //Omega_group = Lcorr_group * Lcorr_group ; // Correlation matrix of ancestry distribution
  //Omega_ind = Lcorr_ind * Lcorr_ind ; // Correlation matrix of ind re distribution
}
