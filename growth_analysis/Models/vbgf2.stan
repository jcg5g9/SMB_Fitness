// Single parameter VBGF with sex/river effects saved as vbgf2.stan
data {
  // Length-at-age data
  int<lower=0> Nobs;                    // number of observations 
  int<lower=0> Nages;                    // number of observations 
  real length[Nobs];                    // length
  real age[Nobs];                       // length
  
  // Individual level data
  int<lower=0> Nind;                    // number of individuals
  int<lower=0> Ncoef;                   // Number of predictors
  matrix[Nind, Ncoef] X;                // Design matrix
  matrix[4,2] Xhat;                     // Prediction matrix
  real q[Nind];                         // Proportion SMB
  
  int<lower=1, upper=Nind> id[Nobs];    // Sample ID
}
parameters {
  // VBGF Params ----
  real mu_linf;                         // asymptotic length
  real mu_k;                            // growth coef
  real mu_t0;                           // age at length 0

  // Regressors
  vector[Ncoef] beta_linf;
  vector[Ncoef] beta_k;
  vector[Ncoef] beta_t0;
 
  real<lower=0> sigma;                  // observation error
}
transformed parameters {
  // Predicted length
  vector[Nobs] length_hat;
  
  // Individual level parameters
  vector[Nind] linf_ind = exp(X * beta_linf) * mu_linf; 
  vector[Nind] k_ind = exp(X * beta_k) * mu_k;
  vector[Nind] t0_ind = X * beta_t0 + mu_t0;
  
  // Predicted length
  for(i in 1:Nobs){
    length_hat[i] = linf_ind[id[i]] * (1-exp(-k_ind[id[i]]* (age[i] - (t0_ind[id[i]]))));
  }
}
model {
  // Priors
  // - Global parameters
  // - Starks, T. A., & Rodger, A. W. (2020). Otolith and scale‐based growth standards for lotic Smallmouth Bass. North American Journal of Fisheries Management, 40(4), 986-994.
  mu_linf ~ lognormal(log(578), 0.02432599);
  mu_k ~ lognormal(log(0.125), 0.04580929);
  mu_t0 ~ normal(-1.79, 0.0625);

  // Regessor priors
  beta_linf ~ normal(0, 1);
  beta_k ~ normal(0, 1);
  beta_t0 ~ normal(0, 1);

  // Likelihood
  length ~ normal(length_hat, sigma);
} 
generated quantities{
  // Predicted length
  matrix[5, Nages] length_pred;
  
  // Parameters for prediction
  vector[4] linf_pred = exp(Xhat * beta_linf) * mu_linf; 
  vector[4] k_pred = exp(Xhat * beta_k) * mu_k;
  vector[4] t0_pred = Xhat * beta_t0 + mu_t0;


  for(i in 1:Nages){
    length_pred[1,i] = mu_linf * (1 - exp(-mu_k * (i - mu_t0))); // Global
    length_pred[2,i] = linf_pred[1] * (1 - exp(-k_pred[1] * (i - t0_pred[1]))); // Females river 1
    length_pred[3,i] = linf_pred[2] * (1 - exp(-k_pred[2] * (i - t0_pred[2]))); // Male river 1
    length_pred[4,i] = linf_pred[3] * (1 - exp(-k_pred[3] * (i - t0_pred[3]))); // Females river 2
    length_pred[5,i] = linf_pred[4] * (1 - exp(-k_pred[4] * (i - t0_pred[4]))); // Males river 2
  }
}
