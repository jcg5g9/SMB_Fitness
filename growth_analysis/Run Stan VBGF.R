# TODO

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
