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
  
  real amax;                            // max age
  real amin;                            // min age
  
  int<lower=1, upper=Nind> id[Nobs];    // Sample ID
}
parameters {
  // VBGF Params
  real mu_lmin;                         // estimated length at min age
  real mu_lmax;                         // estimated length at max age
  real mu_k;                            // growth coef
  real<lower=0> sigma;                  // observation error
}
transformed parameters {
  // Global parameters (natural scale)
  real mu_linf = mu_lmin + (mu_lmax - mu_lmin)/(1-exp(-mu_k * (amax - amin)));
  real mu_t0 = log(mu_linf/(mu_linf-mu_lmin))/(-mu_k) + amin;
  
  // Linf <- L1 + (L2-L1)/(1-exp(-K*(A2-A1)))
  // A0 <- log(Linf/(Linf-L1))/-K + A1
  
  // Predicted length
  vector[Nobs] length_hat = mu_linf * (1-exp(-mu_k * (age - mu_t0)));
}
model {
  // Priors
  // - Global parameters
  // - Starks, T. A., & Rodger, A. W. (2020). Otolith and scale‐based growth standards for lotic Smallmouth Bass. North American Journal of Fisheries Management, 40(4), 986-994.
  mu_lmax ~ lognormal(log(400), 0.05);
  mu_lmin ~ uniform(0, mu_lmax);
  mu_k ~ lognormal(log(0.125), 0.04580929);
  
  // Likelihood
  length ~ lognormal(log(length_hat), sigma);
} 
generated quantities{
  // Predicted length
  matrix[Nages, 1] length_pred;
  
  for(i in 1:Nages){
    length_pred[i,1] = mu_linf * (1 - exp(-mu_k * (i - mu_t0))); // Global
  }
}
