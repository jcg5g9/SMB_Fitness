# Fitting hierarchical growth models using simulated length-at-age

##### Setup #####
library(rstan)
library(dplyr)
rstan_options(threads_per_chain = 1)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)
cols <- c("#86BBD8", "#F6AE2D", "#F26419","#2F4858") # Colors for the VBGF lines (Females/Males)

##### Assign simulated data to list ##### 
source("growth_analysis/Simulate VBGF Q Data.R")

##### Model 1 #####
# - Single parameter model
fit1 <- stan(
  file = "growth_analysis/Models/vbgf1.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 2000,          # number of warmup iterations per chain
  iter = 4000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = list(max_treedepth = 12, adapt_delta = 0.9)
)

saveRDS(fit1, file = "growth_analysis/Models/Fits/vbgf_sim_fit1.rds")

# - Print and plot MCMC
print(fit1, pars=c("mu_linf", "mu_k", "mu_t0"), probs=c(.1,.5,.9))
traceplot(fit1, pars = c("mu_linf", "mu_k", "mu_t0"), inc_warmup = FALSE, nrow = 2)
pairs(fit1, pars = c("mu_linf", "mu_k", "mu_t0"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit1, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

# - Plot fitted model
plot(y = dat$length , x = dat$age, ylab = "Total length (mm)", xlab = "Age (yr)", cex = 1, cex.lab = 1.25, pch = 20, main = "Model 1")
draws <- as.data.frame(fit1)

# - Plot median curve
length_at_age <-draws[,grepl("length_pred",colnames(draws))]
median.lengths <- apply(length_at_age, 2, median)
lines(1:max(age), median.lengths, col = 1, lty = 1, lwd = 2)


##### Model 2 #####
# - Single parameter model with sex/river effects
fit2 <- stan(
  file = "growth_analysis/Models/vbgf2.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 2000,          # number of warmup iterations per chain
  iter = 4000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = list(max_treedepth = 12, adapt_delta = 0.9)
)

saveRDS(fit2, file = "growth_analysis/Models/Fits/vbgf_sim_fit2.rds")

# - Print and plot MCMC
print(fit2, pars=c("mu_linf", "mu_k", "mu_t0", "beta_linf", "beta_k", "beta_t0"), probs=c(.1,.5,.9)) # None of the betas are sig
traceplot(fit2, pars = c("mu_linf", "mu_k", "mu_t0"), inc_warmup = FALSE, nrow = 2)
pairs(fit2, pars = c("mu_linf", "mu_k", "mu_t0"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit2, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

# - Plot fitted model
plot(y = dat$length , x = dat$age, ylab = "Total length (mm)", xlab = "Age (yr)", cex = 2, cex.lab = 1.25, 
     col = cols[factor(sample.sex)], pch = c(17, 19)[factor(sample.river)], main = "Model 2")
draws <- as.data.frame(fit2)

# - Plot median curve
laa_global <-draws[,grepl("length_pred\\[1",colnames(draws))] # Global
laa_fem_r1 <-draws[,grepl("length_pred\\[2,",colnames(draws))] # Females river 1
laa_mal_r1 <-draws[,grepl("length_pred\\[3,",colnames(draws))] # Males river 1
laa_fem_r2 <-draws[,grepl("length_pred\\[4,",colnames(draws))] # Females river 2
laa_mal_r2 <-draws[,grepl("length_pred\\[5,",colnames(draws))] # Males river 2

lines(1:max(age), apply(laa_global, 2, median), col = 1, lty = 1, lwd = 3) # Essentially same as females river 1
lines(1:max(age), apply(laa_fem_r1, 2, median), col = cols[1], lty = 1, lwd = 2)
lines(1:max(age), apply(laa_mal_r1, 2, median), col = cols[2], lty = 1, lwd = 2)
lines(1:max(age), apply(laa_fem_r2, 2, median), col = cols[1], lty = 2, lwd = 2)
lines(1:max(age), apply(laa_mal_r2, 2, median), col = cols[2], lty = 2, lwd = 2)
legend("bottomright", c("Females river 1", "Females river 2","Males river 1", "Males river 2"), col = rep(cols, each = 2), lty = c(1,2,1,2), bty = "n", lwd = 2)


##### Model 3 #####
# -  VBGF model with sex/river effects and ancestry level random effects
fit3 <- stan(
  file = "growth_analysis/Models/vbgf3.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 2000,          # number of warmup iterations per chain
  iter = 4000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = list(max_treedepth = 12, adapt_delta = 0.9)
)

saveRDS(fit3, file = "growth_analysis/Models/Fits/vbgf_sim_fit3.rds")

# - Print and plot MCMC
print(fit3, pars=c("mu_linf", "mu_k", "mu_t0", "beta_linf", "beta_k", "beta_t0", "linf_lineage", "k_lineage", "t0_lineage"), probs=c(.1,.5,.9)) # None of the betas are sig
traceplot(fit3, pars = c("mu_linf", "mu_k", "mu_t0"), inc_warmup = FALSE, nrow = 2)
pairs(fit3, pars = c("mu_linf", "mu_k", "mu_t0"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit3, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

# - Plot fitted model
plot(y = dat$length , x = dat$age, ylab = "Total length (mm)", xlab = "Age (yr)", cex = 2, cex.lab = 1.25, 
     col = cols[(sample.sex+1)*2-1 + (sample.river)], pch = c(17, 19)[factor(sample.sex)], main = "Model 3")
draws <- as.data.frame(fit3)

# - Plot median curve
lines(1:max(age), apply(draws[,grepl("length_pred\\[1",colnames(draws))], 2, median), col = 1, lty = 1, lwd = 4) # Global

# - SMB Females
lines(1:max(age), apply(draws[,grepl("length_pred\\[2,",colnames(draws))], 2, median), col = cols[1], lty = 1, lwd = 2) # Females river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[4,",colnames(draws))], 2, median), col = cols[31], lty = 2, lwd = 2) # Females river 2

# - Neosho Females
lines(1:max(age), apply(draws[,grepl("length_pred\\[6,",colnames(draws))], 2, median), col = cols[2], lty = 1, lwd = 2) # Females river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[8,",colnames(draws))], 2, median), col = cols[2], lty = 2, lwd = 2) # Females river 2

# - SMB males
lines(1:max(age), apply(draws[,grepl("length_pred\\[3,",colnames(draws))], 2, median), col = cols[3], lty = 1, lwd = 2) # Males river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[5,",colnames(draws))], 2, median), col = cols[3], lty = 2, lwd = 2) # Males river 2

# - Neosho Males
lines(1:max(age), apply(draws[,grepl("length_pred\\[7,",colnames(draws))], 2, median), col = cols[4], lty = 1, lwd = 2) # Males river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[9,",colnames(draws))], 2, median), col = cols[4], lty = 2, lwd = 2) # Males river 2

legend("bottomright", c("SMB Females", "Neosho Females","SMB Males", "Neosho Males", "River 1", "River 2"), col = c(cols,1,1), lty = c(1,1,1,1,1,2), bty = "n", lwd = 2)


##### Model 4 #####
# -  VBGF model with sex/river effects and ancestry and individual level random effects
fit4 <- stan(
  file = "growth_analysis/Models/vbgf4.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 4000,          # number of warmup iterations per chain
  iter = 7000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = list(max_treedepth = 12, adapt_delta = 0.9)
)

saveRDS(fit4, file = "growth_analysis/Models/Fits/vbgf_sim_fit4.rds")

# - Print and plot MCMC
print(fit4, pars=c("mu_linf", "mu_k", "mu_t0", "beta_linf", "beta_k", "beta_t0", "linf_lineage", "k_lineage", "t0_lineage"), probs=c(.1,.5,.9)) # None of the betas are sig
traceplot(fit4, pars = c("mu_linf", "mu_k", "mu_t0"), inc_warmup = FALSE, nrow = 2)
pairs(fit4, pars = c("mu_linf", "mu_k", "mu_t0"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit4, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

# - Plot fitted model
plot(y = dat$length , x = dat$age, ylab = "Total length (mm)", xlab = "Age (yr)", cex = 0.5, cex.lab = 1.25, 
     col = cols[(sample.sex+1)*2-1 + (sample.river)], pch = c(17, 19)[factor(sample.sex)], main = "Model 4")
draws <- as.data.frame(fit4)

# - Plot median curve
lines(1:max(age), apply(draws[,grepl("length_pred\\[1",colnames(draws))], 2, median), col = 1, lty = 1, lwd = 4) # Global

# - SMB Females
lines(1:max(age), apply(draws[,grepl("length_pred\\[2,",colnames(draws))], 2, median), col = cols[1], lty = 1, lwd = 2) # Females river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[4,",colnames(draws))], 2, median), col = cols[31], lty = 2, lwd = 2) # Females river 2

# - Neosho Females
lines(1:max(age), apply(draws[,grepl("length_pred\\[6,",colnames(draws))], 2, median), col = cols[2], lty = 1, lwd = 2) # Females river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[8,",colnames(draws))], 2, median), col = cols[2], lty = 2, lwd = 2) # Females river 2

# - SMB males
lines(1:max(age), apply(draws[,grepl("length_pred\\[3,",colnames(draws))], 2, median), col = cols[3], lty = 1, lwd = 2) # Males river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[5,",colnames(draws))], 2, median), col = cols[3], lty = 2, lwd = 2) # Males river 2

# - Neosho Males
lines(1:max(age), apply(draws[,grepl("length_pred\\[7,",colnames(draws))], 2, median), col = cols[4], lty = 1, lwd = 2) # Males river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[9,",colnames(draws))], 2, median), col = cols[4], lty = 2, lwd = 2) # Males river 2

legend("bottomright", c("SMB Females", "Neosho Females","SMB Males", "Neosho Males", "River 1", "River 2"), col = c(cols,1,1), lty = c(1,1,1,1,1,2), bty = "n", lwd = 2)


# * Test Lineage Parameters #####
# - No sig difference
par(mfrow = c(1,3))
hist(draws$`linf_lineage[1]`-draws$`linf_lineage[2]`, xlab = "Linf diff", main = NA)
hist(draws$`k_lineage[1]`-draws$`k_lineage[2]`, xlab = "K diff", main = "Model 4")
hist(draws$`t0_lineage[1]`-draws$`t0_lineage[2]`, xlab = "t0 diff", main = NA)
par(mfrow = c(1,1))



##### Model 5 #####
# -  VBGF model with  ancestry and individual level random effects
# - No river/sex effects
fit5 <- stan(
  file = "growth_analysis/Models/vbgf5.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 3000,          # number of warmup iterations per chain
  iter = 5000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = list(max_treedepth = 12, adapt_delta = 0.9)
)

saveRDS(fit5, file = "growth_analysis/Models/Fits/vbgf_sim_fit5.rds")

# - Print and plot MCMC
print(fit5, pars=c("mu_linf", "mu_k", "mu_t0", "linf_lineage", "k_lineage", "t0_lineage"), probs=c(.1,.5,.9)) # None of the betas are sig
traceplot(fit5, pars = c("mu_linf", "mu_k", "mu_t0"), inc_warmup = FALSE, nrow = 2)
pairs(fit5, pars = c("mu_linf", "mu_k", "mu_t0"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit5, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

# - Plot fitted model
plot(y = dat$length , x = dat$age, ylab = "Total length (mm)", xlab = "Age (yr)", cex = 2, cex.lab = 1.25, 
     col = cols[factor(round(dat$q))], pch = 16, main = "Model 5")
draws <- as.data.frame(fit5)

# - Plot median curve
lines(1:max(age), apply(draws[,grepl("length_pred\\[1",colnames(draws))], 2, median), col = 1, lty = 1, lwd = 4) # Global

# - SMB
lines(1:max(age), apply(draws[,grepl("length_pred\\[2,",colnames(draws))], 2, median), col = cols[1], lty = 1, lwd = 2)

# - Neosho 
lines(1:max(age), apply(draws[,grepl("length_pred\\[3,",colnames(draws))], 2, median), col = cols[2], lty = 1, lwd = 2)

legend("bottomright", c("Global", "SMB", "Neosho"), col = c(1, cols[1:2]), lty = c(1,1,1), bty = "n", lwd = 2)


#* Test Lineage Parameters #####
# - No sig difference
par(mfrow = c(1,3))
hist(draws$`linf_lineage[1]`-draws$`linf_lineage[2]`, xlab = "Linf diff", main = NA)
hist(draws$`k_lineage[1]`-draws$`k_lineage[2]`, xlab = "K diff", main = "Model 5")
hist(draws$`t0_lineage[1]`-draws$`t0_lineage[2]`, xlab = "t0 diff", main = NA)
par(mfrow = c(1,1))



##### Model 6 #####
# -  VBGF model with  ancestry FIXED and individual level random effects
# - No river/sex effects
fit6 <- stan(
  file = "growth_analysis/Models/vbgf6.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 3000,          # number of warmup iterations per chain
  iter = 5000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = list(max_treedepth = 12, adapt_delta = 0.9)
)

saveRDS(fit6, file = "growth_analysis/Models/Fits/vbgf_sim_fit6.rds")

# - Print and plot MCMC
print(fit6, pars=c("linf_lineage", "k_lineage", "t0_lineage"), probs=c(.1,.5,.9)) # None of the betas are sig
traceplot(fit6, pars = c("linf_lineage", "k_lineage", "t0_lineage"), inc_warmup = FALSE, nrow = 2)
pairs(fit6, pars = c("linf_lineage", "k_lineage", "t0_lineage"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit6, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

# - Plot fitted model
plot(y = dat$length , x = dat$age, ylab = "Total length (mm)", xlab = "Age (yr)", cex = 2, cex.lab = 1.25, 
     col = cols[factor(round(dat$q))], pch = 16, main = "Model 6")
draws <- as.data.frame(fit6)

# - SMB
lines(1:max(age), apply(draws[,grepl("length_pred\\[1,",colnames(draws))], 2, median), col = cols[1], lty = 1, lwd = 2)

# - Neosho 
lines(1:max(age), apply(draws[,grepl("length_pred\\[2,",colnames(draws))], 2, median), col = cols[3], lty = 1, lwd = 2)

legend("bottomright", c("SMB", "Neosho"), col = c(cols[c(1,3)]), lty = c(1,1), bty = "n", lwd = 2)


# * Test Lineage Parameters #####
# - No sig difference
par(mfrow = c(1,3))
hist(draws$`linf_lineage[1]`-draws$`linf_lineage[2]`, xlab = "Linf diff", main = NA)
hist(draws$`k_lineage[1]`-draws$`k_lineage[2]`, xlab = "K diff", main = "Model 6")
hist(draws$`t0_lineage[1]`-draws$`t0_lineage[2]`, xlab = "t0 diff", main = NA)
par(mfrow = c(1,1))

