//  VBGF with sex/river effects saved as vbgf4.stan
//  Includes ancestry level random effects
functions{
  real normal_lb_rng(real mu, real sigma, real lb) {
    real p = normal_cdf(lb, mu, sigma);  // cdf for bounds
    real u = uniform_rng(p, 1);
    return (sigma * inv_Phi(u)) + mu;  // inverse cdf for value
  }
}
data {
  // * Length-at-age data ----
  int<lower=0> Nobs;                    // number of observations 
  int<lower=0> Nages;                   // number of observations 
  real length[Nobs];                    // length
  real age[Nobs];                       // length
  
  // * Priors ----
  real eta_scale_prior;                    // scale cauchy prior
  real cholesky_prior;                  // prior for lkj 
  real beta_scale;                      // scale beta sigma
  
  // * Individual level data ----
  int<lower=0> Nind;                    // number of individuals
  int<lower=0> Ncoef;                   // Number of predictors
  matrix[Nind, Ncoef] X;                // Design matrix
  matrix[4,Ncoef] Xhat;                     // Prediction matrix
  vector[Nind] q;                       // Proportion SMB
  
  int<lower=1, upper=Nind> id[Nobs];    // Sample ID
  vector[3] Zero;                       // Vector of zero for mean
}
parameters {
  // * VBGF Params ----
  real<lower=0> mu_linf;                // asymptotic length
  real<lower=0> mu_k;                   // growth coef
  real mu_t0;                           // age at length 0
  matrix[2, 3] eta_lineage;             // Ancestry level deviation
  matrix[Nind, 3] eta_ind;              // Individual deviation from VBGF parameters
  
  // * Regressors ----
  vector[Ncoef] beta_linf;
  vector[Ncoef] beta_k;
  vector[Ncoef] beta_t0;
  
  // * Variance ----
  vector<lower=0>[3] sigma_ind;        // Prior scale for individual-level variation
  vector<lower=0>[3] sigma_group;      // Prior scale for group-level variation
  cholesky_factor_corr[3] Lcorr;       // Prior correlation for group- and individual-level variation
  real<lower=0> sigma;                // observation error
}
transformed parameters {
  // * Correlation matrices ----
  matrix[3, 3] Omega = Lcorr * Lcorr; // Correlation matrix of ancestry distribution
  
  // * Predicted length ----
  vector[Nobs] length_hat;
  
  // * Group level parameters ----
  vector[2] linf_lineage = mu_linf * exp(eta_lineage[,1]);
  vector[2] k_lineage = mu_k * exp(eta_lineage[,2]);
  vector[2] t0_lineage = mu_t0 + eta_lineage[,3];
  
  real post_linf_diff = linf_lineage[1] - linf_lineage[2];
  real post_k_diff = k_lineage[1] - k_lineage[2];
  real post_t0_diff = t0_lineage[1] - t0_lineage[2];
  
  real eta_smb_linf = eta_lineage[1,1];
  real eta_smb_k = eta_lineage[1,2];
  real eta_smb_t0 = eta_lineage[1,3];
  
  real eta_n_linf = eta_lineage[2,1];
  real eta_n_k = eta_lineage[2,2];
  real eta_n_t0 = eta_lineage[2,3];
  
  // * Individual level parameters ----
  vector[Nind] linf_ind = mu_linf * exp(X * beta_linf + eta_smb_linf * q + eta_n_linf * (1-q)); 
  vector[Nind] k_ind = mu_k * exp(X * beta_k + eta_smb_k * q + eta_n_k * (1-q));
  vector[Nind] t0_ind = mu_t0 + X * beta_t0 + eta_smb_t0 * q + eta_n_t0 * (1-q);
  
  // * Predicted length ----
  for(i in 1:Nobs){
    length_hat[i] = linf_ind[id[i]] * exp(eta_ind[id[i],1]) * (1-exp(-k_ind[id[i]] * exp(eta_ind[id[i],2]) * (age[i] - (t0_ind[id[i]] + eta_ind[id[i],3]))));
  }
}
model {
  // * Priors ----
  // - Global parameters
  // - Starks, T. A., & Rodger, A. W. (2020). Otolith and scale‐based growth standards for lotic Smallmouth Bass. North American Journal of Fisheries Management, 40(4), 986-994.
  mu_linf ~ lognormal(log(578), 0.02432599 * 10);
  mu_k ~ lognormal(log(0.125), 0.04580929 * 10);
  mu_t0 ~ normal(-1.79, 0.0625 * 10);
  
  // - Ancestry level variation priors
  sigma_group ~ normal(0, eta_scale_prior);
  Lcorr ~ lkj_corr_cholesky(cholesky_prior); // Centered around 0 https://mjskay.github.io/ggdist/reference/lkjcorr_marginal.html
  
  for(i in 1:2){
    eta_lineage[i,] ~ multi_normal_cholesky(Zero, diag_pre_multiply(sigma_group, Lcorr));
  }
  
  // - Individual variation priors
  sigma_ind ~ normal(0, eta_scale_prior);
  
  for(i in 1:Nind){
    eta_ind[i,] ~ multi_normal_cholesky(Zero, diag_pre_multiply(sigma_ind, Lcorr));
  }
  
  // - Regessor priors
  beta_linf ~ normal(0, beta_scale);
  beta_k ~ normal(0, beta_scale);
  beta_t0 ~ normal(0, beta_scale);
  
  // - Observation error
  sigma ~ normal(0, 5);
  
  
  // * Likelihood ----
  length ~ normal(length_hat, sigma);
} 
generated quantities{
  // * Prior predicted ----
  // - Only doing group level
  
  // - VBGF Params 
  real prior_mu_linf = lognormal_rng(log(578), 0.02432599 * 10);      // asymptotic length
  real prior_mu_k = lognormal_rng(log(0.125), 0.04580929 * 10);       // growth coef
  real prior_mu_t0 = normal_rng(-1.79, 0.0625 * 10);                  // age at length 0
  matrix[2, 3] prior_eta_lineage;                                // Ancestry level deviation
  // matrix[Nind, 3] prior_eta_ind;                                 // Individual deviation from VBGF parameters
  
  // - Regressors
  vector[Ncoef] prior_beta_linf;
  vector[Ncoef] prior_beta_k;
  vector[Ncoef] prior_beta_t0;
  
  // - Variance 
  vector[3] prior_sigma_group;                                        // Prior scale for group-level variation
  matrix[3, 3] prior_Lcorr = lkj_corr_cholesky_rng(3, cholesky_prior);// Prior correlation for group-level variation
  real prior_sigma = normal_lb_rng(0, 5, 0);                          // observation error (left bounded at 0)
  
  // - Group level transformed parameters
  vector[2] prior_linf_lineage;
  vector[2] prior_k_lineage;
  vector[2] prior_t0_lineage;
  
  real prior_eta_smb_linf;
  real prior_eta_smb_k;
  real prior_eta_smb_t0;
  
  real prior_eta_n_linf;
  real prior_eta_n_k;
  real prior_eta_n_t0;
  
  
  // * Posterior/prior predicted length ----
  matrix[9, Nages] pred_length;
  matrix[9, Nages] prior_length;
  
  // - Parameters for prediction
  vector[4] pred_linf = exp(Xhat * beta_linf) * mu_linf; 
  vector[4] pred_k = exp(Xhat * beta_k) * mu_k;
  vector[4] pred_t0 = Xhat * beta_t0 + mu_t0;
  
  vector[4] prior_linf;
  vector[4] prior_k;
  vector[4] prior_t0; 
  
  
  // * Simulate prior ----
  for(i in 1:3){
    prior_sigma_group[i] = normal_lb_rng(0, eta_scale_prior, 0); // Left bounded at 0
  }
  
  for(i in 1:Ncoef){
    prior_beta_linf[i] = normal_rng(0, beta_scale);
    prior_beta_k[i] = normal_rng(0, beta_scale);
    prior_beta_t0[i] = normal_rng(0, beta_scale);
  }
  
  for(i in 1:2){
    prior_eta_lineage[i,] = to_row_vector(multi_normal_cholesky_rng(Zero, diag_pre_multiply(prior_sigma_group, prior_Lcorr)));
  }
  
  // - Group level transformed parameters
  prior_linf_lineage = prior_mu_linf * exp(prior_eta_lineage[,1]);
  prior_k_lineage = prior_mu_k * exp(prior_eta_lineage[,2]);
  prior_t0_lineage = prior_mu_t0 + prior_eta_lineage[,3];
  
  real prior_linf_diff = prior_linf_lineage[1] - prior_linf_lineage[2];
  real prior_k_diff = prior_k_lineage[1] - prior_k_lineage[2];
  real prior_t0_diff = prior_t0_lineage[1] - prior_t0_lineage[2];
  
  prior_eta_smb_linf = prior_eta_lineage[1,1];
  prior_eta_smb_k = prior_eta_lineage[1,2];
  prior_eta_smb_t0 = prior_eta_lineage[1,3];
  
  prior_eta_n_linf = prior_eta_lineage[2,1];
  prior_eta_n_k = prior_eta_lineage[2,2];
  prior_eta_n_t0 = prior_eta_lineage[2,3];
  
  prior_linf = exp(Xhat * prior_beta_linf) * prior_mu_linf; 
  prior_k = exp(Xhat * prior_beta_k) * prior_mu_k;
  prior_t0 = Xhat * prior_beta_t0 + prior_mu_t0;
  
  
  // * Predict length (posterior and prior) ----
  for(i in 1:Nages){
    
    // ** Global ----
    pred_length[1,i] = mu_linf * (1 - exp(-mu_k * (i - mu_t0))); // Global posterior
    prior_length[1,i] = mu_linf * (1 - exp(-mu_k * (i - mu_t0))); // Global prior
    
    
    // ** SMB ----
    // - Females river 1
    pred_length[2,i] = pred_linf[1] * exp(eta_smb_linf) * (1 - exp(-pred_k[1] * exp(eta_smb_k) * (i - pred_t0[1] - eta_smb_t0))); // Post
    prior_length[2,i] = prior_linf[1] * exp(prior_eta_smb_linf) * (1 - exp(-prior_k[1] * exp(prior_eta_smb_k) * (i - prior_t0[1] - prior_eta_smb_t0))); // Prior
    
    // - Male river 1
    pred_length[3,i] = pred_linf[2] * exp(eta_smb_linf) * (1 - exp(-pred_k[2] * exp(eta_smb_k) * (i - pred_t0[2] - eta_smb_t0))); // Post
    prior_length[3,i] = prior_linf[2] * exp(prior_eta_smb_linf) * (1 - exp(-prior_k[2] * exp(prior_eta_smb_k) * (i - prior_t0[2] - prior_eta_smb_t0))); // Prior
    
    // - Females river 2
    pred_length[4,i] = pred_linf[3] * exp(eta_smb_linf) * (1 - exp(-pred_k[3] * exp(eta_smb_k) * (i - pred_t0[3] - eta_smb_t0))); // Post
    prior_length[4,i] = prior_linf[3] * exp(prior_eta_smb_linf) * (1 - exp(-prior_k[3] * exp(prior_eta_smb_k) * (i - prior_t0[3] - prior_eta_smb_t0))); // Prior
    
    // - Males river 2
    pred_length[5,i] = pred_linf[4] * exp(eta_smb_linf) * (1 - exp(-pred_k[4] * exp(eta_smb_k) * (i - pred_t0[4] -eta_smb_t0))); // Post
    prior_length[5,i] = prior_linf[4] * exp(prior_eta_smb_linf) * (1 - exp(-prior_k[4] * exp(prior_eta_smb_k) * (i - prior_t0[4] - prior_eta_smb_t0))); // Prior
    
    
    // ** Neosho ----
    // - Females river 1
    pred_length[6,i] = pred_linf[1] * exp(eta_n_linf) * (1 - exp(-pred_k[1] * exp(eta_n_k) * (i - pred_t0[1] - eta_n_t0))); // Post
    prior_length[6,i] = prior_linf[1] * exp(prior_eta_n_linf) * (1 - exp(-prior_k[1] * exp(prior_eta_n_k) * (i - prior_t0[1] - prior_eta_n_t0))); // Prior
    
    // - Male river 1
    pred_length[7,i] = pred_linf[2] * exp(eta_n_linf) * (1 - exp(-pred_k[2] * exp(eta_n_k) * (i - pred_t0[2] - eta_n_t0))); // Post
    prior_length[7,i] = prior_linf[2] * exp(prior_eta_n_linf) * (1 - exp(-prior_k[2] * exp(prior_eta_n_k) * (i - prior_t0[2] - prior_eta_n_t0))); // Prior
    
    // - Females river 2
    pred_length[8,i] = pred_linf[3] * exp(eta_n_linf) * (1 - exp(-pred_k[3] * exp(eta_n_k) * (i - pred_t0[3] - eta_n_t0))); // Post
    prior_length[8,i] = prior_linf[3] * exp(prior_eta_n_linf) * (1 - exp(-prior_k[3] * exp(prior_eta_n_k) * (i - prior_t0[3] - prior_eta_n_t0))); // Prior
    
    // - Males river 2
    pred_length[9,i] = pred_linf[4] * exp(eta_n_linf) * (1 - exp(-pred_k[4] * exp(eta_n_k) * (i - pred_t0[4] -eta_n_t0))); // Post
    prior_length[9,i] = prior_linf[4] * exp(prior_eta_n_linf) * (1 - exp(-prior_k[4] * exp(prior_eta_n_k) * (i - prior_t0[4] - prior_eta_n_t0))); // Prior
  }
}
