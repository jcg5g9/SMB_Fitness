

//
data {
  int<lower=0> N; // Number of individuals
  vector[N] smb; // proportion of smb 
  vector[N] con; // body condition index
}

parameters {
  real beta1;
  real alpha;
  real <lower=0> sigma;
  //real beta2;
}

model {
  // Priors
  alpha ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  //beta2 ~ normal(0, 100);
  
  // Likelihood
  // Distribution for random part
  for (n in 1:N){
    con[n] ~ normal(alpha + beta1*smb[n], sigma);
  }
  
}




