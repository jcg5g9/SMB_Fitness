// Global VBGF model saved as vbgf1.stan
data {
  // Length-at-age data
  int<lower=0> Nobs;                    // number of observations 
  int<lower=0> Nages;                   // number of observations 
  real length[Nobs];                    // length
  vector[Nobs] age;                     // length
  
  // Individual level data
  int<lower=0> Nind;                    // number of individuals
  int<lower=0> Ncoef;                   // Number of predictors
  matrix[Nind, Ncoef] X;                // Design matrix
  real q[Nind];                         // Proportion SMB
  
  int<lower=1, upper=Nind> id[Nobs];    // Sample ID
}
parameters {
  // VBGF Params
  real mu_linf;                         // asymptotic length
  real mu_k;                            // growth coef
  real mu_t0;                           // age at length 0
  real<lower=0> sigma;                  // observation error
}
transformed parameters {
  // Predicted length
  vector[Nobs] length_hat = mu_linf * (1-exp(-mu_k * (age - mu_t0)));
}
model {
  // Priors
  // - Global parameters
  // - Starks, T. A., & Rodger, A. W. (2020). Otolith and scale‐based growth standards for lotic Smallmouth Bass. North American Journal of Fisheries Management, 40(4), 986-994.
  mu_linf ~ lognormal(log(578), 0.02432599);
  mu_k ~ lognormal(log(0.125), 0.04580929);
  mu_t0 ~ normal(-1.79, 0.0625);
  
  // Likelihood
  length ~ normal(length_hat, sigma);
} 
generated quantities{
  // Predicted length
  matrix[Nages, 1] length_pred;
  
  for(i in 1:Nages){
    length_pred[i,1] = mu_linf * (1 - exp(-mu_k * (i - mu_t0))); // Global
  }
}
