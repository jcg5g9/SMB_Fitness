//  VBGF with ancestry fixed effects
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
  // Group level parameters
  vector[G] linf_lineage ;
  vector[G] k_lineage;
  vector[G] t0_lineage;
  
  // Variance ----
  real<lower=0> sigma;                  // observation error
}
transformed parameters {
  // Predicted length
  vector[Nobs] length_hat;
  
  // Predicted length
  for(i in 1:Nobs){
    length_hat[i] = linf_lineage[group[i]] * (1-exp(-k_lineage[group[i]] * (age[i] - (t0_lineage[group[i]]))));
  }
}
model {
  // Priors
  // - Global parameters
  // - Starks, T. A., & Rodger, A. W. (2020). Otolith and scale‐based growth standards for lotic Smallmouth Bass. North American Journal of Fisheries Management, 40(4), 986-994.
  linf_lineage ~ lognormal(log(578), 0.02432599);
  k_lineage ~ lognormal(log(0.125), 0.04580929);
  t0_lineage ~ normal(-1.79, 0.0625);
  
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