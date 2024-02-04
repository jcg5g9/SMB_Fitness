// condition model 3 - with interactions
//
data {
  int<lower=0> N; // Number of individuals
  vector[N] smb; // proportion of smb 
  vector[N] con; // body condition index
  vector[N] stream; // river
  vector[N] sex;// sex
}

parameters {
  real beta1;
  real beta2;
  real beta3;
  real beta4;
  real alpha;
  real <lower=0> sigma;
}

model {
  // Priors
  alpha ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  beta2 ~ normal(0, 100);
  beta3 ~ normal(0, 100);
  beta4 ~ normal(0, 100);
  
  // Likelihood
  // Distribution for random part
  for (n in 1:N){
    con[n] ~ normal(alpha + beta1*smb[n] + beta2*stream[n] +beta3*sex[n] + beta4*(smb[n]*stream[n]), sigma);
  }
  
}



