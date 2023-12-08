
library(jagsUI)

#################################################################################### 
# JAGS VBGM ######################################################################## 
####################################################################################


##### SPECIFY JAGS MODEL CODE #####
jags.mod =
  "model{

## Data likelihood ----
for(i in 1:N){
  length[i] ~ dnorm(y.hat[i], tau.y) # May want this to be lognormal
  y.hat[i] = Linf.i[sample.id[i]] * (1-exp(-k.i[sample.id[i]] * (age[i] - t0.i[sample.id[i]])))
}

# - Observation error
tau.y = pow(sigma, -2)
sigma ~ dunif(0,100)


## Group level random effects ----
for (j in 1:G){ # - Loop through groups
  group.re.pars[j,1:3] ~ dmnorm(mu.re, group.inv.cov.mat)  # The `0` may need to be a vector of 3 zeros?

  log.Linf.group.re[j] = group.re.pars[j,1]
  log.k.group.re[j] = group.re.pars[j,2] 
  t0.group.re[j] = group.re.pars[j,3]
  
  # - Get group level params on natural scale for reporting
  Linf.group[j] = exp(log.mu.Linf + log.Linf.group.re[j])
  k.group[j] = exp(log.mu.k + log.k.group.re[j])
  t0.group[j] = mu.t0 + t0.group.re[j]
}

# - Group level covariance matrix
group.inv.cov.mat ~ dwish( z.group.mat, 4 ) # Add `z.group.mat = diag(x=1,nrow=3)` to the data
group.cov.mat = inverse( group.inv.cov.mat ) # you can get out the variance and correlation coeficients using `cov2cor` in R or add code here


## Individual level random effects ----
for (k in 1:Nind){ # - Loop through individuals
  ind.re.pars[k,1:3] ~ dmnorm(mu.re, ind.inv.cov.mat) # The `0` may need to be a vector of 3 zeros?

  log.Linf.ind.re[k] = ind.re.pars[k,1]
  log.k.ind.re[k] = ind.re.pars[k,2] 
  t0.ind.re[k] = ind.re.pars[k,3]
  
  # - Get individual level params on natural scale 
  # NOTE: group is now of length Nind and will have to be adjusted in the data
  Linf.i[k] = exp(log.mu.Linf + log.Linf.group.re[group[k]] + log.Linf.ind.re[k] + log.linf.beta.sex * sex[k] + log.linf.beta.river * river[k]) # Exponent here to keep Linf positive
  k.i[k] = exp(log.mu.k + log.k.group.re[group[k]] + log.k.ind.re[k] + log.k.beta.sex * sex[k] + log.k.beta.river * river[k]) # Exponent here to keep K positive
  t0.i[k] = mu.t0 + t0.group.re[group[k]] + t0.ind.re[k] + t0.beta.sex * sex[k] + t0.beta.river * river[k]
}

# - Individual level covariance matrix
ind.inv.cov.mat ~ dwish( z.ind.mat, 4) # Add `z.ind.mat = diag(x=1,nrow=3)` to the data
ind.cov.mat = inverse( ind.inv.cov.mat ) # you can get out the variance and correlation coeficients using =`cov2cor` in R or add code here

# Sex and river effet priors
log.linf.beta.river ~ dnorm(0, 0.01)
log.k.beta.river ~ dnorm(0, 0.01)
t0.beta.river ~ dnorm(0, 0.01)

log.linf.beta.sex ~ dnorm(0, 0.01)
log.k.beta.sex ~ dnorm(0, 0.01)
t0.beta.sex ~ dnorm(0, 0.01)

## Population level parameters ----
# - These could probably be more informative given previous studies
log.mu.Linf ~ dmnorm(0,0.0001)
log.mu.k ~ dunif(0,1)
mu.t0 ~ dmnorm(0,0.0001)

# - Put on natural scale for reporting
mu.Linf = exp(log.mu.Linf)
mu.k = exp(log.mu.k)
 
}"

# write model to a text file
writeLines(jags.mod, "jags_model.txt")

##### PARAMETERS TO MONITOR #####
params = c("mu.Linf", "mu.k", "mu.t0", "Linf.group", "k.group", "t0.group","Linf.i", "k.i", "t0.i", 
           "ind.cov.mat","group.cov.mat","sigma", 
           "log.linf.beta.river", "log.k.beta.river", "t0.beta.river", "log.linf.beta.sex", "log.k.beta.sex", "t0.beta.sex" )


##### Assign data to list ##### 
# Real data
load('growth_analysis/data/bc_data/full_bc_data.rda')
unique(as.factor(full_bc_data$ancestry_group))
# [1] Admixed         Neosho_Bass     Smallmouth_Bass
# Levels: Admixed Neosho_Bass Smallmouth_Bass

age <- full_bc_data$annulus
length <- full_bc_data$bc_tl
group <- as.numeric(as.factor(full_bc_data$ancestry_group)) 
river <- as.numeric(as.factor(full_bc_data$river)) #[1] big sugar [2] elk river
sex <- as.numeric(as.factor(full_bc_data$sex)) #[1] female [2] male
sample.id <- as.numeric(as.factor(full_bc_data$sample_id))
river[river==3] <- 2
Nind <- length(unique(full_bc_data$sample_id))


dat = list(age = age, length = length, sex = sex, river = river, sample.id = sample.id, 
           N = Nind, G = 3, Nind = Nind, group = group, z.group.mat = diag(x=1,nrow=3), 
           z.ind.mat = diag(x=1,nrow=3), mu.re = rep(0, 3))


##### MCMC DIMENSIONS #####
ni = 500000 # this was how many it took to sort of converge with simulated data (and runJags)
#ni=5000
nb = 100000
na = 100000
nt = 1000
nc = 3
n.iter = ni + nb



##### RUN THE MODEL IN JAGS #####
runJagsOut <- jags( model.file="jags_model.txt" ,
                    parameters.to.save =params ,
                    data=dat ,
                    n.chains=nc ,
                    n.adapt=na ,
                    n.burnin=nb ,
                    n.iter=n.iter ,
                    n.thin=nt ,
                    parallel=T)
print(runJagsOut)
