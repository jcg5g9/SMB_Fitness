# Ancestry-Condition correlation
# SMB Fitness paper
# SJC
# Oct 2023

library(rstan)
library(lme4)
library(jagsUI)

# check data
# load condition data
head(condition_data)
cd <- condition_data

# Basic model
# condition ~ smb

con <- cd$condition * 100
smb <- as.numeric(scale(cd$smb))
n <- length(c)
riv <- as.numeric(as.factor(cd$river))
riv[which(riv==3)] <- 2 # recode bc 1 and 3 was causing issues with jags
sex <- as.numeric(as.factor(cd$sex))
# make riv and sex 0 or 1 for stan
riv1 <- riv-1
sex1 <- sex-1

dat <- list(con=con, smb=smb, N=n, riv=riv1, sex=sex1)


# frequentist models
test1 <- lm(condition ~ smb, data=cd)
summary(test1)

test2 <- lm(condition ~ smb + river + sex, data=cd)
summary(test2)

# significant negative effect of smb on condition, positive effect of elk vs big sugar, no sex difference

# Bayesian - Jags (testing because stan turned out weird and I think I screwed up)
# MCMC settings 
ni <- 5000 

nt <- 5 
nb <- 2000 
nc <- 3 

sink("m1.jags") 
cat("
    # likelihood
    model{
     for (i in 1:N){
       con[i] ~ dnorm(mu[i],tau)
       mu[i] <- beta1*smb[i] + beta2[riv[i]] + beta3[sex[i]]+ alpha 
     }
    
    # Priors
     alpha ~ dnorm(0,0.001) 
     beta1 ~ dnorm(0,0.001) 
     beta2[1] <- 0
     beta2[2]  ~ dnorm(0,0.001)
     beta3[1] <- 0 
     beta3[2] ~ dnorm(0,0.001) 
    
     tau <- pow(sig, -2) # precision for mu normal distribution
     sig ~ dunif(0, 10) # variance (to get precision for mu normal distribution, this is a hyperprior) 
    }
    ",fill = TRUE)

sink() 

# put data in a list to give to jags ("bundle data")
jags.data <- list(con=con, smb=smb, N=n, riv=riv, sex=sex)

# Initial values
inits <- function (){list(beta2= c(NA, 0), beta3=c(NA, 0), alpha = 0)}

# Parameters 

parameters <- c("beta1","beta2", "beta3","alpha","sig")


# Call jags from R
m1 <- jags(data=jags.data, inits=inits, parameters.to.save=parameters, model.file="m1.jags", 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=F)


# Summarize posteriors
print(m1, digits = 3) 

# Jags results: significant negative effect of smb, positive effect of elk vs big sugar, no sex difference


# Bayesian - Stan

fit1 <- stan(
  file = "ac_analysis/condition_model_test_1.stan",  
  data = dat,    
  chains = 4,             
  warmup = 2000,          
  iter = 4000,            
  cores = 4,              
)

print(fit1,digits=5)

fit2 <- stan(
  file = "ac_analysis/condition_model_2.stan",  
  data = dat,    
  chains = 4,             
  warmup = 2000,          
  iter = 4000,            
  cores = 4,              
)

print(fit2,digits=5)

# Stan resultsL significant negative effect of smb, positive effect of elk vs big sugar, no sex difference

# looks consistent across all 