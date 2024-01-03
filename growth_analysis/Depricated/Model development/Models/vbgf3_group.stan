//  VBGF with sex/river effects saved as vbgf3_group.stan
//  Includes ancestry level random effects as groups
data {
  // Length-at-age data
  int<lower=0> Nobs;                    // number of observations 
  int<lower=0> G;                       // number of Groups 
  int<lower=0> Nages;                   // number of observations 
  real length[Nobs];                    // length
  real age[Nobs];                       // length
  
  // Individual level data
  int<lower=0> Nind;                    // number of individuals
  int<lower=0> Ncoef;                   // Number of predictors
  matrix[Nind, Ncoef] X;                // Design matrix
  matrix[4,2] Xhat;                     // Prediction matrix
  
  int<lower=1, upper=G> group[Nobs];    // Sample group
  int<lower=1, upper=Nind> id[Nobs];    // Sample ID
  vector[3] Zero;                       // Vector of zero for mean
}
parameters {
  // VBGF Params ----
  real mu_linf;                         // asymptotic length
  real mu_k;                            // growth coef
  real mu_t0;                           // age at length 0
  matrix[G, 3] eta_lineage;             // Ancestry level deviation
  
  // Regressors
  vector[Ncoef] beta_linf;
  vector[Ncoef] beta_k;
  vector[Ncoef] beta_t0;
  
  // Variance ----
  cholesky_factor_corr[3] Lcorr_group;           // Prior correlation for group-level variation
  vector<lower=0>[3] sigma_group;         // Prior scale for group-level variation
  
  real<lower=0> sigma;                  // observation error
}
transformed parameters {
  // Predicted length
  vector[Nobs] length_hat;
  
  // Group level parameters
  vector[G] linf_lineage = mu_linf * exp(eta_lineage[,1]);
  vector[G] k_lineage = mu_k * exp(eta_lineage[,2]);
  vector[G] t0_lineage = mu_t0 + eta_lineage[,3];
  
  // Individual level parameters
  vector[Nind] linf_ind = mu_linf * exp(X * beta_linf); 
  vector[Nind] k_ind = mu_k * exp(X * beta_k);
  vector[Nind] t0_ind = mu_t0 + X * beta_t0;
  
  // Predicted length
  for(i in 1:Nobs){
    length_hat[i] = linf_ind[id[i]] * exp(eta_lineage[group[i], 1]) * (1-exp(-k_ind[id[i]] * exp(eta_lineage[group[i], 2]) * (age[i] - (t0_ind[id[i]] + eta_lineage[group[i], 3]))));
  }
}
model {
  // Priors
  // - Global parameters
  // - Jackson, Z. J., Quist, M. C., & Larscheid, J. G. (2008). Growth standards for nine North American fish species. Fisheries Management and Ecology, 15(2), 107-118.
  mu_linf ~ lognormal(log(498.6), 0.4);
  mu_k ~ lognormal(log(0.229), 0.5);
  mu_t0 ~ normal(0.141, 1);
  
  // - Ancestry level variation priors
  sigma_group ~ cauchy(0, 0.5);
  Lcorr_group ~ lkj_corr_cholesky(1);
  
  for(i in 1:G){
    eta_lineage[i,] ~ multi_normal_cholesky(Zero, diag_pre_multiply(sigma_group, Lcorr_group));
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
  matrix[1+(G*4), Nages] length_pred;
  
  // Parameters for prediction (Note need to adjust for more predictors)
  vector[4] linf_pred = exp(Xhat * beta_linf) * mu_linf; 
  vector[4] k_pred = exp(Xhat * beta_k) * mu_k;
  vector[4] t0_pred = Xhat * beta_t0 + mu_t0;
  
  // Global
  for(i in 1:Nages){
    length_pred[1,i] = mu_linf * (1 - exp(-mu_k * (i - mu_t0))); // Global
  }
  
  // Loop through groups
  int ind = 2;
  for(g in 1:G){
    // Loop through ages
    for(i in 1:Nages){
      length_pred[ind,i] = linf_pred[1] * exp(eta_lineage[g,1]) * (1 - exp(-k_pred[1] * exp(eta_lineage[g,2]) * (i - t0_pred[1] - eta_lineage[g,3]))); // Females river 1
      length_pred[ind+1,i] = linf_pred[2] * exp(eta_lineage[g,1]) * (1 - exp(-k_pred[2] * exp(eta_lineage[g,2]) * (i - t0_pred[2] - eta_lineage[g,3]))); // Male river 1
      length_pred[ind+2,i] = linf_pred[3] * exp(eta_lineage[g,1]) * (1 - exp(-k_pred[3] * exp(eta_lineage[g,2]) * (i - t0_pred[3] - eta_lineage[g,3]))); // Females river 2
      length_pred[ind+3,i] = linf_pred[4] * exp(eta_lineage[g,1]) * (1 - exp(-k_pred[4] * exp(eta_lineage[g,2]) * (i - t0_pred[4] - eta_lineage[g,3]))); // Males river 2
    }
    ind += 4; // Update index
  }
}
