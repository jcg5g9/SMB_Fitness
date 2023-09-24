# This code simulates length-at-age data for individuals from a population where there are multiple sub-groups and there are sex and location (two rivers) specific differences in growth. Multiple length-at-age samples are "simulated/sampled" from each individual, representing capture-recapture data.

# The data are then fit to a hierarchical von bertalanfy growth curve with the following structure:

# Length.obs.i = Linf.i * (1-exp(-K.i * (age.obs.i - t0.i))) + eps.i
# Linf_i = exp(mu.linf + linf.group + linf.i + beta.linf.river * river.i + beta.linf.sex * sex.i)
# K_i  = exp(mu.k + k.group + k.i + beta.k.river * river.i + beta.k.sex * sex.i)
# t0_i = mu.t0 + t0.group + t0.i + beta.t0.river * river.i + beta.t0.sex * sex.i

# Where the "betas" are river and sex effects depending on the river of origin and sex of individual "i" 
# The group and and individual i level parameters are assumed to be MVN with mean 0 and conjugate prior variance
# Eps.i is the normally distributed error term


####################################################################################
# Simulate length-at-age data ######################################################
####################################################################################
# Set Seed
library(MASS)
library(dplyr)
set.seed(1234)

## Specify data size
Nind = 200
G <- 3
ages = seq(from=1,to=20, by = 1)

group <- sample(1:G, Nind, replace = TRUE) # Percent SMB for individual X
nsamples <- sample(1:10, Nind, replace = TRUE) # Randomly select number of samples per individual X
sex <- sample(0:1, Nind, replace = TRUE) # Randomly select sex of individual X
river <- sample(0:1, Nind, replace = TRUE) # Randomly select river of individual X
N <- sum(nsamples) # Total number of samples
sample.id <- rep(1:Nind, times = nsamples) # ID of individual X from sample Y
sample.sex <- rep(sex, times = nsamples) # sex of individual X from sample Y
sample.river <- rep(river, times = nsamples) # river of individual X from sample Y
sample.group <- rep(group, times = nsamples) # Group of individual X from sample Y


## Mu VBGM hyperparameters 
mu.Linf = 500
mu.k = 0.3 
mut.t0 = 0.5
mu.parms <- c(mu.Linf, mu.k, mut.t0)
sigma = 10 # Observation error

# River and sex effects (e.g. difference for sex 1 and river 1 from sex 0 and river 0, respectively)
loglinf.beta.sex <- 0.05
loglinf.beta.river <- 0.01

logk.beta.sex <- 0.01
logk.beta.river <- 0.1

t0.beta.sex <- 0.1
t0.beta.river <- 0.1


## Group level random effects
sigma.group = c(.03, 0.05, 0.2) #NOTE: Sigma is important
rho = 0.3 # Correlation between group level parameters
cor.group.mat = matrix(rho, 3, 3)
diag(cor.group.mat) <- 2
cov.group.mat <- diag(sigma.group) %*% cor.group.mat %*% diag(sigma.group) # Get covariance


## Individual level random effects
sigma.ind = c(0.1, 0.1, 0.1)
rho = 0.3 # Correlation between group level parameters
cor.ind.mat = matrix(rho, 3, 3)
diag(cor.ind.mat) <- 1
cov.ind.mat <- diag(sigma.ind) %*% cor.ind.mat %*% diag(sigma.ind) # Get covariance


## Simulate parameters for groups/individuals ----
# - Empty matrix and vectors to fill with parameters and data, respectively
group.param.mat <- group.re.mat <- matrix(NA,G,3,byrow = T)
ind.param.mat <- ind.re.mat <- matrix(NA,Nind,3,byrow = T)

# - Random effects
colnames(group.re.mat) <- c("log.Linf.group.re", "log.k.group.re", "t0.group.re")
colnames(ind.re.mat) <- c("log.Linf.ind.re", "log.k.ind.re", "t0.in.re")

# - On VBGF scalge
colnames(group.param.mat) <- c("Linf.group", "k.group", "t0.group")
colnames(ind.param.mat) <- c("Linf.ind", "k.ind", "t0.ind")


# - Simulate group level parameters
for(i in 1:G){
  group.re.mat[i,] <- mvrnorm(1, rep(0,3), cov.group.mat)
  group.param.mat[i,1:2] <- mu.parms[1:2] * exp(group.re.mat[i,1:2]) # Log to natural scale
  group.param.mat[i,3] <- mu.parms[3] + group.re.mat[i,3]
}

# - Simulate individual level parameters
for(i in 1:Nind){
  ind.re.mat[i,] <- mvrnorm(1, rep(0,3), cov.ind.mat)
  ind.param.mat[i,1] <- mu.parms[1] * exp(group.re.mat[group[i],1] + ind.re.mat[i,1] + loglinf.beta.sex * sex[i] + loglinf.beta.river * river[i]) # Log to natural scale Linf
  ind.param.mat[i,2] <- mu.parms[2] * exp(group.re.mat[group[i],2] + ind.re.mat[i,2] + logk.beta.sex * sex[i] + logk.beta.river * river[i]) # Log to natural scale K
  ind.param.mat[i,3] <- mu.parms[3] + group.re.mat[group[i],3] + ind.re.mat[i,3] + t0.beta.sex * sex[i] + t0.beta.river * river[i]
}

## Simulate length-at-age data
age = c()
length = c()
ind <- 1
for(i in 1:Nind) {
  ages_temp = sample(ages, nsamples[i], replace = FALSE) # Sample random age from age range w/o replacement
  age[ind:(ind+nsamples[i]-1)] <- ages_temp
  length[ind:(ind+nsamples[i]-1)] = (ind.param.mat[i,1] * (1 - exp(-ind.param.mat[i,2]*(ages_temp-ind.param.mat[i,3])))) + rnorm(nsamples[i],0,sigma)
  
  ind <- ind + nsamples[i]
}


# Assign data to data frame ----

# - Assign to list
dat = list(
  G = G,
  Nobs = length(length),
  Nages = round(max(age)),
  length = length,
  age = age,
  Zero = rep(0, 3),
  
  Nind = Nind,
  Ncoef = 2,
  group = group,
  X = model.matrix(~ sex + river)[,2:3], # No intercept
  Xhat = matrix(c(0,0,1,0,0,1,1,1), nrow = 4, ncol = 2, byrow = TRUE), # Matrix for prediction
  
  id = sample.id
)

print(group.param.mat)
