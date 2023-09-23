//  VBGF with ancestry fixed and individual level random effects
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
  
  int<lower=1, upper=Nind> id[Nobs];    // Sample ID
  vector[3] Zero;                       // Vector of zero for mean
}
parameters {
  // VBGF Params ----
  
  // Group level parameters
  vector[2] linf_lineage;
  vector[2] k_lineage;
  vector[2] t0_lineage;
  
  // Variance ----
  real<lower=0> sigma;                  // observation error
}
transformed parameters {
  // Predicted length
  vector[Nobs] length_hat;

  real mu_smb_linf = linf_lineage[1];
  real mu_smb_k = k_lineage[1];
  real mu_smb_t0 = t0_lineage[1];
  
  real mu_n_linf = linf_lineage[2];
  real mu_n_k = k_lineage[2];
  real mu_n_t0 = t0_lineage[2];
  
  // Individual level parameters
  vector[Nind] linf_ind =  mu_smb_linf * q + mu_n_linf * (1-q); 
  vector[Nind] k_ind = mu_smb_k * q + mu_n_k * (1-q);
  vector[Nind] t0_ind = mu_smb_t0 * q + mu_n_t0 * (1-q);
  
  // Predicted length
  for(i in 1:Nobs){
    length_hat[i] = linf_ind[id[i]] * (1-exp(-k_ind[id[i]] * (age[i] - (t0_ind[id[i]]))));
  }
}
model {
  // Priors
  // - Global parameters
  // - Jackson, Z. J., Quist, M. C., & Larscheid, J. G. (2008). Growth standards for nine North American fish species. Fisheries Management and Ecology, 15(2), 107-118.
  linf_lineage ~ lognormal(log(498.6), 0.4);
  k_lineage ~ lognormal(log(0.229), 0.5);
  t0_lineage ~ normal(0.141, 1);
  
  // Likelihood
  length ~ normal(length_hat, sigma);
} 
generated quantities{
  // Predicted length
  matrix[2, Nages] length_pred;
  
  for(i in 1:Nages){
    
    // SMB
    length_pred[1,i] = (mu_smb_linf) * (1 - exp(-(mu_smb_k) * (i - mu_smb_t0))); 
    
    // Neosho
    length_pred[2,i] = (mu_n_linf) * (1 - exp(-(mu_n_k) * (i - mu_n_t0))); 
  }
}
