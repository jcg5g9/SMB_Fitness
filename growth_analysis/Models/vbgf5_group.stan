//  VBGF with sex/river effects saved as vbgf5_group.stan
//  Includes ancestry level random effects as groups and individual random effects
//  No sex/river effect
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
  matrix[Nind, 3] eta_ind;              // Individual deviation from VBGF parameters
  
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
  vector[G] linf_lineage = mu_linf * exp(eta_lineage[,1]);
  vector[G] k_lineage = mu_k * exp(eta_lineage[,2]);
  vector[G] t0_lineage = mu_t0 + eta_lineage[,3];
  
  // Predicted length
  for(i in 1:Nobs){
    length_hat[i] = mu_linf * exp(eta_ind[id[i],1] + eta_lineage[group[i], 1]) * (1-exp(-mu_k * exp(eta_ind[id[i],2]  + eta_lineage[group[i], 2]) * (age[i] - (mu_t0 + eta_ind[id[i],3] + eta_lineage[group[i], 3]))));
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
  
  // - Individual variation priors
  sigma_ind ~ cauchy(0, 0.5);
  Lcorr_ind ~ lkj_corr_cholesky(2);
  
  for(i in 1:Nind){
    eta_ind[i,] ~ multi_normal_cholesky(Zero, diag_pre_multiply(sigma_ind, Lcorr_ind));
  }
  
  // Likelihood
  length ~ normal(length_hat, sigma);
} 
generated quantities{
  // Predicted length
  matrix[1+G, Nages] length_pred;
  
  // Global
  for(i in 1:Nages){
    length_pred[1,i] = mu_linf * (1 - exp(-mu_k * (i - mu_t0))); // Global
  }
  
  // Loop through groups
  for(g in 1:G){
    // Loop through ages
    for(i in 1:Nages){
      length_pred[g+1,i] = mu_linf * exp(eta_lineage[g,1]) * (1 - exp(-mu_k * exp(eta_lineage[g,2]) * (i - mu_t0 - eta_lineage[g,3]))); 
    }
  }
}
