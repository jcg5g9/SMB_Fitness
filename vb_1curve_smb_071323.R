# von  Bertalanffy equation - Jags
# 1 curve per total, with individual random effect and smb covariate


# code from G Adams github (modifed a bit)
# https://github.com/grantdadams/Growth-Models



library(boot)
library(jagsUI)
library(ggplot2)

setwd('G:\\My Drive\\other_peoples_projects\\joe bass fitness')


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Data Prep
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#  ******** SMB Data ********

load('full_bc_data.rda')
unique(as.factor(full_bc_data$ancestry_group))
# [1] Admixed         Neosho_Bass     Smallmouth_Bass
# Levels: Admixed Neosho_Bass Smallmouth_Bass

age2 <- full_bc_data$annulus
length2 <- full_bc_data$bc_tl
#group2 <- as.numeric(as.factor(full_bc_data$ancestry_group)) # [1] admixed [2] neosho [3] smallmouth
group3 <- rep(1, length(age2)) # this makes all in one group so you can run all at once without changing the model
smb <- as.numeric(full_bc_data$smb)
smb_s <- as.numeric(scale(full_bc_data$smb))
river <- as.numeric(as.factor(full_bc_data$river)) #[1] big sugar [2] elk river
sex <- as.numeric(as.factor(full_bc_data$sex)) #[1] female [2] male
id <- as.numeric(as.factor(full_bc_data$sample_id))
river[river==3] <- 2

# Data set for models!
bc_data <-  list(age = age2, length = length2, length2=length2,group = group3, 
                 river=river, sex=sex, id=id, smb=smb_s,
                 N = length(age2), G = length(unique(group3)))

# group sample sizes - smb = 44, nb = 134, adm = 159

# MODEL OUTPUT IS SAVED IN GITHUB, CAN SKIP TO LINE 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Jags Model
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# MCMC settings
ni = 50000
nb = 20000
na = 20000
nt = 5
nc = 3

setwd('C:\\Users\\sarah.clements1\\Documents\\GitHub\\SMB_Fitness')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  01 Stardard model, 3 genetic groups (3 curves) random individual effect
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# JAGS VBGM modified from G Adams Github 
# This converges with groups
# changed prior for mu.k and removed group-specific t0
# added random group effect


sink('vb_1curve_smb_071323.jags')
cat(
  "model{
for(i in 1:N){
length[i] ~ dnorm(y.hat[i], tau.y)
y.hat[i] = Linf[group[i]] * (1-exp(-k[group[i]] * (age[i] - t0))) + betaSMB*smb[i] + idr[id[i]]
}
 
# SD
tau.y = pow(sigma, -2)
sigma ~ dunif(0,100)
 
# Level-2 parameters
for(j in 1:G){ 
Linf[j] ~ dnorm(mu.Linf, tau.Linf)
k[j] ~ dnorm(mu.k, tau.k)
}
t0 ~ dnorm(mu.t0, tau.t0)

for (q in 1:N){
idr[q] ~ dnorm(0,tau.idr)
}


# Priors for level-2 parameters
log.mu.Linf ~ dnorm(0,0.0001)
log.mu.k ~ dunif(0,1)
mu.t0 ~ dnorm(0,0.0001)

betaSMB ~ dnorm(0,0.0001)
 
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

tau.idr <- pow(sig.idr, -2) # tau and sig for random effect
sig.idr ~ dunif(0, 10)

}"
)
sink()

# changed log.mu.k to ~dunif(0,1) following Schofeld et al. 2013 https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.12069

#PARAMETERS TO MONITOR
params = c("betaSMB","Linf", "k", "t0", "mu.Linf", "mu.k", "mu.t0", "mu.Linf", 
           "mu.k", "mu.t0", "sig.Linf","sig.k","sig.t0","sig.idr","sigma")
inits = function(){list()}


# dat = simulated data from G Adams, bc_data = Joe's data w/back calculated tl
vbc <- jagsUI::jags(data=bc_data, inits=inits, parameters.to.save=params, model.file="vb_1curve_smb_071323.jags", 
                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel=T)
print(vbc, digits=3)
# a lot of credible interval overlap in both parameters among the 3 groups
jagsUI::traceplot(vbc)

save(vbc, file='vb_1curve_smb_071323.rda')

load('C:\\Users\\sarah.clements1\\Documents\\GitHub\\SMB_Fitness\\vb_1curve_smb_071323.rda')


# plot curves
# equation
# Y ~ Linf[group[i]] * (1-exp(-k[group[i]] * (age[i] - t0))) 

# pull out posterior samples (using capital letters for model generated data)
T0 <- c(vbc$samples[[1]][,which(row.names(vbc$summary)=='t0')], 
        vbc$samples[[2]][,which(row.names(vbc$summary)=='t0')], 
        vbc$samples[[3]][,which(row.names(vbc$summary)=='t0')])


LINF1 <-  c(vbc$samples[[1]][,which(row.names(vbc$summary)=='Linf')], 
            vbc$samples[[2]][,which(row.names(vbc$summary)=='Linf')], 
            vbc$samples[[3]][,which(row.names(vbc$summary)=='Linf')])



K1 <-  c(vbc$samples[[1]][,which(row.names(vbc$summary)=='k')], 
         vbc$samples[[2]][,which(row.names(vbc$summary)=='k')], 
         vbc$samples[[3]][,which(row.names(vbc$summary)=='k')])






# organize data for prediction plots (in base r, could be done in ggplot just as easy)
LINF <- data.frame(LINF1=LINF1)
K <- data.frame(K1=K1)

age_seq <- seq(min(age2), max(age2), by=1)
age_seq <- seq(min(age2), 16, by=1)

c1  <- data.frame(mean=rep(NA, length(age_seq)), lci=rep(NA, length(age_seq)), hci=rep(NA, length(age_seq)))

curves <- c1



# make curves
# I may or may not be making the curve correctly (quantile linf and t0 too or just k)??

# way #1 calculate size once per age with quantiles of parameter estimates
for (i in 1:1){
  for (j in 1:length(age_seq)){
    curves$mean[j] <- mean(LINF[,i]) * (1-exp(-(mean(K[,i])) * (age_seq[j] - mean(T0))))
    curves$lci[j] <- quantile(LINF[,i],probs=0.025) * (1-exp(-(quantile(K[,i], probs=0.025)) * (age_seq[j] - quantile(T0, probs=0.025))))
    curves$hci[j] <- quantile(LINF[,i],probs=0.975) * (1-exp(-(quantile(K[,i], probs=0.975)) * (age_seq[j] - quantile(T0, probs=0.975))))
  }
} 


# make transparent colors
col2rgb('gray80')
#trpurple <- rgb(147,112,219,90,maxColorValue=255)
#trblue <- rgb(0,191,255,90,maxColorValue=255)
trgray80 <- rgb(204,204,204,75,maxColorValue=255)
trgray20 <- rgb(51,51,51,75, maxColorValue=255)
trgray <- rgb(0,0,0,75, maxColorValue=255)
# look at curves
plot(x=age_seq, y=curves$mean, type='n', ylim=c(100,450),
     xlab="Age", ylab="Total Length") 
# raw data
points(x=bc_data$age, y=bc_data$length, pch=16, col=trgray)
# ribbons
polygon(x=c(age_seq,rev(age_seq)), y=c(curves$lci,rev(curves$hci)), col=trgray20, border=NA)

# lines
lines(x=age_seq, y=curves$mean, type='l', col='gray20', lwd=2) 

# linf lines
abline(h=mean(LINF[,1]), lwd=2, lty='dotted', col='gray20')

