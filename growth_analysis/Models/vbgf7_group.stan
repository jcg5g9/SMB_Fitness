//  VBGF with ancestry level random effects
// No individual and sex/river effects
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
  
  // Predicted length
  for(i in 1:Nobs){
    length_hat[i] = linf_lineage[group[i]] * (1-exp(-k_lineage[group[i]] * (age[i] - t0_lineage[group[i]])));
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
  Lcorr_group ~ lkj_corr_cholesky(1);
  
  for(i in 1:G){
    eta_lineage[i,] ~ multi_normal_cholesky(Zero, diag_pre_multiply(sigma_group, Lcorr_group));
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
