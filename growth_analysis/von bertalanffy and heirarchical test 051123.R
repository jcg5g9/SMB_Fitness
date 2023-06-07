# von  Bertalanffy equation test - Jags
# code from G Adams github (modifed to run with jagsUI)
# https://github.com/grantdadams/Growth-Models
# the code on his website didn't work, the one on his github did


library(boot)
library(jagsUI)

setwd('G:\\My Drive\\other_peoples_projects\\joe bass fitness')


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Data Prep
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# *****Simulate length-at-age data - From G. Adams Github *****

# Set Seed
set.seed(1234)

# Number of Groups
Groups=6

# True VBGM hyperparameters c(mean, sd)
true.Linf = c(500,20)
true.k = c(.3,.05)
true.t0 = c(1.5,.4)

sigma <- 15

# Age range
ages = seq(from=1,to=15, by = .05)

# Empty matrix and vectors to fill with parameters and data, respectively
param.mat = matrix(NA,Groups,3,byrow = T)
colnames(param.mat) <- c("Linf", "k", "t0")
ctr = 0
age = c()
length = c()
group = c()

# Simulate group level parameters
for(i in 1:Groups){
  param.mat[i,1] = rnorm(1,true.Linf[1],true.Linf[2]) # Assign group level Linf
  param.mat[i,2] = rnorm(1,true.k[1],true.k[2]) # Assign group level k
  param.mat[i,3] = rnorm(1,true.t0[1],true.t0[2]) # Assign group level t0
  n.samples = sample(200:1000, 1) # Number of samples per group s
  
  # Simulate data
  for(j in 1:n.samples) {
    ctr = ctr + 1 # Indexing variable
    age[ctr] = sample(ages, 1) # Sample randon age from age range
    length[ctr] = (param.mat[i,1] * (1 - exp(-param.mat[i,2]*(age[ctr]-param.mat[i,3])))) + rnorm(1,0,sigma)
    group[ctr] = i
  }
}

# Assign data to data frame
dat = data.frame(age = age, length = length, group = group, N = length(age), G = length(unique(group)))
dat <- dat[which(dat$length > 0),]

# Plot the data
for (i in 1:length(unique(dat$group))){
  sub = dat[which(dat$group==i),]
  if ( i == 1){
    plot(sub$age, sub$length, col = i, xlab = "Age", ylab = "Length")
  } else {
    points(sub$age, sub$length , col = i)
  }
}

# Assign data to list
dat = list(age = age, length = length, group = group, N = length(age), G = length(unique(group)))


#  ******** SMB Data ********
load('growth_analysis/data/bc_data/full_bc_data.rda')
unique(as.factor(full_bc_data$ancestry_group))
# [1] Admixed         Neosho_Bass     Smallmouth_Bass
# Levels: Admixed Neosho_Bass Smallmouth_Bass

age2 <- full_bc_data$annulus
length2 <- full_bc_data$bc_tl
group2 <- as.numeric(as.factor(full_bc_data$ancestry_group)) # it didn't converge when separating the groups, probably there is a way to fix..
group3 <- rep(1, length(age2)) # this makes all in one group so you can run all at once without changing the model
smb <- as.numeric(full_bc_data$smb)
smb_s <- as.numeric(scale(full_bc_data$smb))
river <- as.numeric(as.factor(full_bc_data$river)) #[1] big sugar [2] elk river
sex <- as.numeric(as.factor(full_bc_data$sex)) #[1] female [2] male
id <- as.numeric(as.factor(full_bc_data$sample_id))
river[river==3] <- 2

# Data set for models!
bc_data <-  list(age = age2, length = length2, length2=length2, smb=smb_s, group = group3, 
                 river=river, sex=sex, id=id,
                 N = length(age2), G = length(unique(group3)))

# group sample sizes - smb = 44, nb = 134, adm = 159

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Jags Models
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# MCMC settings
ni = 50000
nb = 20000
na = 20000
nt = 1
nc = 3

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  01 JAGS VBGM from G Adams Github "test" and "test0" 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This one runs with the simulated data and converges
# Also runs and converges with Joe's data without groups
# Doesn't quite converge with Joe's data with groups but it's close - and probably fixable
# for analysis with covariates wouldn't want to include smb and ancestry group in the same model anyway

# SPECIFY JAGS MODEL CODE 

sink('vb_test.jags')
cat(
  "model{
for(i in 1:N){
length[i] ~ dnorm(y.hat[i], tau.y)
y.hat[i] = Linf[group[i]] * (1-exp(-k[group[i]] * (age[i] - t0[group[i]] )))
}
 
# SD
tau.y = pow(sigma, -2)
sigma ~ dunif(0,100)
 
# Level-2 parameters
for(j in 1:G){ 
Linf[j] ~ dnorm(mu.Linf, tau.Linf)
k[j] ~ dnorm(mu.k, tau.k)
t0[j] ~ dnorm(mu.t0, tau.t0)
}
 
# Priors for level-2 parameters
log.mu.Linf ~ dnorm(0,0.0001)
log.mu.k ~ dnorm(0,0.0001)
mu.t0 ~ dnorm(0,0.0001)
 
# Get hyperparameters on untransformed scale
mu.Linf = exp(log.mu.Linf)
mu.k = exp(log.mu.k)
 
# Precision
tau.Linf = pow(sig.Linf,-2)
tau.k = pow(sig.k,-2)
tau.t0 = pow(sig.t0,-2)
 
# SD of parameters
sig.Linf ~ dunif(0,10)
sig.k ~ dunif(0,10)
sig.t0 ~ dunif(0,10)
}"
)
sink()


#PARAMETERS TO MONITOR
params = c("Linf", "k", "t0", "mu.Linf", "mu.k", "mu.t0", "mu.Linf", 
           "mu.k", "mu.t0", "sig.Linf","sig.k","sig.t0","sigma" )
inits = function(){list()}


# dat = simulated data from G Adams, bc_data = Joe's data w/back calculated tl
test <- jagsUI::jags(data=bc_data, inits=inits, parameters.to.save=params, model.file="vb_test.jags", 
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel=T)
print(test, digits=3)
jagsUI::traceplot(test)

test0 <- jagsUI::jags(data=bc_data, inits=inits, parameters.to.save=params, model.file="vb_test.jags", 
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel=T)
print(test0, digits=3)
jagsUI::traceplot(test0)

# Plot the curve
# Linf*(1-exp(-K*(annulus[i]-t0)))

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 02 Heirarchical Model Attempt "vb_test2"
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sink('vb_test2.jags')
cat(
  "model{
for(i in 1:N){
length[i] ~ dnorm(y.hat[i], tau.y)
y.hat[i] = Linf[group[i]] * (1-exp(-k[group[i]] * (age[i] - t0[group[i]] )))

# residuals from VB equation
res[i] <- (length[i]-y.hat[i])

smb[i] ~ dnorm(mu.r[i], tau.r)
mu.r[i] <- beta1*res[i] +alpha
}


# SD
tau.y = pow(sigma, -2)
sigma ~ dunif(0,100)

tau.r = pow(sigma.r, -2)
sigma.r ~ dunif(0,100)

alpha ~ dnorm(0,0.001)
 
# Level-2 parameters
for(j in 1:G){ 
Linf[j] ~ dnorm(mu.Linf, tau.Linf)
k[j] ~ dnorm(mu.k, tau.k)
t0[j] ~ dnorm(mu.t0, tau.t0)
}

beta1 ~ dnorm(0,0.0001)

# Priors for level-2 parameters
log.mu.Linf ~ dnorm(0,0.0001)
log.mu.k ~ dnorm(0,0.0001)
mu.t0 ~ dnorm(0,0.0001)
 
# Get hyperparameters on untransformed scale
mu.Linf = exp(log.mu.Linf)
mu.k = exp(log.mu.k)
 
# Precision
tau.Linf = pow(sig.Linf,-2)
tau.k = pow(sig.k,-2)
tau.t0 = pow(sig.t0,-2)
 
# SD of parameters
sig.Linf ~ dunif(0,10)
sig.k ~ dunif(0,10)
sig.t0 ~ dunif(0,10)

}"
)
sink()

params = c("beta1","Linf", "k", "t0", "mu.Linf", "mu.k", "mu.t0",
           "mu.Linf", "mu.k", "mu.t0", "sig.Linf","sig.k","sig.t0","sigma","res" )
inits = function(){list()}



# dat = simulated data from G Adams, bc_data = Joe's data w/back calculated tl
test2 <- jagsUI::jags(data=bc_data, inits=inits, parameters.to.save=params, model.file="vb_test2.jags", 
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel=T)
print(test2, digits=3)
jagsUI::traceplot(test2)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 03 putting smb straight into VB equation "vb_test3"
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sink('vb_test3.jags')
cat(
  "model{
for(i in 1:N){
length[i] ~ dnorm(y.hat[i], tau.y)
y.hat[i] = Linf[group[i]] * (1-exp(-k[group[i]] * (age[i] - t0[group[i]] ))) + beta1*smb[i]
}
 
# SD
tau.y = pow(sigma, -2)
sigma ~ dunif(0,100)
 
# Level-2 parameters
for(j in 1:G){ 
Linf[j] ~ dnorm(mu.Linf, tau.Linf)
k[j] ~ dnorm(mu.k, tau.k)
t0[j] ~ dnorm(mu.t0, tau.t0)
}

beta1 ~ dnorm(0,0.0001)
 
# Priors for level-2 parameters
log.mu.Linf ~ dnorm(0,0.0001)
log.mu.k ~ dnorm(0,0.0001)
mu.t0 ~ dnorm(0,0.0001)
 
# Get hyperparameters on untransformed scale
mu.Linf = exp(log.mu.Linf)
mu.k = exp(log.mu.k)
 
# Precision
tau.Linf = pow(sig.Linf,-2)
tau.k = pow(sig.k,-2)
tau.t0 = pow(sig.t0,-2)
 
# SD of parameters
sig.Linf ~ dunif(0,10)
sig.k ~ dunif(0,10)
sig.t0 ~ dunif(0,10)
}"
)
sink()


params = c("beta1","Linf", "k", "t0", "mu.Linf", "mu.k", "mu.t0", "mu.Linf",
           "mu.k", "mu.t0", "sig.Linf","sig.k","sig.t0","sigma" )
inits = function(){list()}

test3 <- jagsUI::jags(data=bc_data, inits=inits, parameters.to.save=params, model.file="vb_test3.jags", 
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel=T)
print(test3, digits=3)
jagsUI::traceplot(test)


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 04 putting all covariates straight into VB equation "vb_test4"
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# it runs, but what does it mean?????

sink('vb_test4.jags')
cat(
  "model{
for(i in 1:N){
length[i] ~ dnorm(y.hat[i], tau.y)
y.hat[i] = (Linf[group[i]] * (1-exp(-k[group[i]] * (age[i] - t0[group[i]] )))) + betaSMB*smb[i] 
 + betaRIV[river[i]] + betaSEX[sex[i]] +
  betaSMBxRIV[river[i]]*smb[i] + betaSMBxSEX[river[i]]*sex[i]
}
 
# SD
tau.y = pow(sigma, -2)
sigma ~ dunif(0,100)
 
# Level-2 parameters
for(j in 1:G){ 
Linf[j] ~ dnorm(mu.Linf, tau.Linf)
k[j] ~ dnorm(mu.k, tau.k)
t0[j] ~ dnorm(mu.t0, tau.t0)
}

betaSMB ~ dnorm(0,0.0001)
betaRIV[1] <- 0 
betaRIV[2] ~ dnorm(0,0.0001) 
betaSEX[1] <- 0 
betaSEX[2] ~ dnorm(0,0.0001)
betaSMBxRIV[1] <- 0 
betaSMBxRIV[2] ~ dnorm(0,0.001) 
betaSMBxSEX[1] <- 0 
betaSMBxSEX[2] ~ dnorm(0,0.001) 
 
# Priors for level-2 parameters
log.mu.Linf ~ dnorm(0,0.0001)
log.mu.k ~ dnorm(0,0.0001)
mu.t0 ~ dnorm(0,0.0001)
 
# Get hyperparameters on untransformed scale
mu.Linf = exp(log.mu.Linf)
mu.k = exp(log.mu.k)
 
# Precision
tau.Linf = pow(sig.Linf,-2)
tau.k = pow(sig.k,-2)
tau.t0 = pow(sig.t0,-2)
 
# SD of parameters
sig.Linf ~ dunif(0,10)
sig.k ~ dunif(0,10)
sig.t0 ~ dunif(0,10)
}"
)
sink()


params = c("betaSMB","betaRIV", "betaSMBxRIV", "betaSEX", "betaSMBxSEX","Linf", "k", "t0", "mu.Linf", 
           "mu.k", "mu.t0", "mu.Linf", "mu.k", "mu.t0", "sig.Linf","sig.k","sig.t0","sigma" )
inits = inits <- function (){list(betaRIV=c(NA, 0), betaSEX=c(NA,0), betaSMBxRIV=c(NA,0), betaSMBxSEX=c(NA,0))}

test4 <- jagsUI::jags(data=bc_data, inits=inits, parameters.to.save=params, model.file="vb_test4.jags", 
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel=T)
print(test4, digits=3)
jagsUI::traceplot(test)