//  VBGF with ancestry fixed and individual level random effects
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
  // Group level parameters
  vector[G] linf_lineage ;
  vector[G] k_lineage;
  vector[G] t0_lineage;
  matrix[Nind, 3] eta_ind;              // Individual deviation from VBGF parameters
  
  // Variance ----
  cholesky_factor_corr[3] Lcorr_ind;          // Prior correlation for individual-level variation
  vector<lower=0>[3] sigma_ind;        // Prior scale for individual-level variation
  
  real<lower=0> sigma;                  // observation error
}
transformed parameters {
  // Predicted length
  vector[Nobs] length_hat;
  
  // Predicted length
  for(i in 1:Nobs){
    length_hat[i] = linf_lineage[group[i]] * exp(eta_ind[id[i],1]) * (1-exp(-k_lineage[group[i]] * exp(eta_ind[id[i],2]) * (age[i] - (t0_lineage[group[i]] + eta_ind[id[i],3]))));
  }
}
model {
  // Priors
  // - Global parameters
  // - Jackson, Z. J., Quist, M. C., & Larscheid, J. G. (2008). Growth standards for nine North American fish species. Fisheries Management and Ecology, 15(2), 107-118.
  linf_lineage ~ lognormal(log(498.6), 0.4);
  k_lineage ~ lognormal(log(0.229), 0.5);
  t0_lineage ~ normal(0.141, 1);
  
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
  matrix[G, Nages] length_pred;
  
  // Loop through groups
  for(g in 1:G){
    // Loop through ages
    for(i in 1:Nages){
      length_pred[g,i] = linf_lineage[g] * (1 - exp(-k_lineage[g] * (i - t0_lineage[g])));
    }
  }
}
