# Fitting hierarchical growth models using back-calculated length-at-age
# Uses group defined by categorization of q

##### Setup #####
library(rstan)
library(dplyr)
rstan_options(threads_per_chain = 1)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)
cols <- c("#86BBD8", "#F6AE2D", "#F26419","#2F4858") # Colors for the VBGF lines (Females/Males)

##### Assign data to list ##### 
# Real data
load('growth_analysis/data/bc_data/full_bc_data.rda')
full_bc_data <- full_bc_data %>%
  mutate(river_code = factor(river_code),
         river_code = ifelse(river_code == 3, 2, river_code),
         ancestry_group = as.numeric(factor(as.character(ancestry_group))),
         sex = as.numeric(factor(sex)))

# - Define data
age <- full_bc_data$annulus
length <- full_bc_data$bc_tl
group <- as.numeric(as.factor(full_bc_data$ancestry_group)) 
sample.id <- as.numeric(as.factor(as.character(full_bc_data$sample_id)))
Nind <- length(unique(sample.id))

# - Set up model matrix (sex/river) where nrow = number of ind
model_mat <- full_bc_data %>%
  group_by(sample_id) %>%
  slice(n()) %>%
  select(sample_id, ancestry_group, river_code, sex) %>%
  mutate(river_code = factor(river_code),
         sex = factor(sex)) %>%
  arrange(sample_id)

# - Assign to list
dat = list(
  Nobs = length(length),
  Nages = round(max(age)),
  G = length(unique(full_bc_data$ancestry_group)),
  length = length,
  age = age,
  Zero = rep(0, 3),
  
  Nind = Nind,
  Ncoef = 2,
  X = model.matrix(~ sex + river_code, model_mat)[,2:3], # No intercept
  Xhat = matrix(c(0,0,1,0,0,1,1,1), nrow = 4, ncol = 2, byrow = TRUE), # Matrix for prediction
  
  group = full_bc_data$ancestry_group,
  id = sample.id
)

##### Model 3 #####
# -  VBGF model with sex/river effects and ancestry group level random effects
fit_group3 <- stan(
  file = "growth_analysis/Models/vbgf3_group.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 2000,          # number of warmup iterations per chain
  iter = 4000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = list(max_treedepth = 12, adapt_delta = 0.9)
)

saveRDS(fit_group3, file = "growth_analysis/Models/Fits/vbgf_fit_group3.rds")

# - Print and plot MCMC
print(fit_group3, pars=c("mu_linf", "mu_k", "mu_t0", "beta_linf", "beta_k", "beta_t0", "linf_lineage", "k_lineage", "t0_lineage"), probs=c(.1,.5,.9)) # None of the betas are sig
traceplot(fit_group3, pars = c("mu_linf", "mu_k", "mu_t0"), inc_warmup = FALSE, nrow = 2)
pairs(fit_group3, pars = c("mu_linf", "mu_k", "mu_t0"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit_group3, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

# - Plot fitted model
plot(y = length , x = age, ylab = "Total length (mm)", xlab = "Age (yr)", cex = 2, cex.lab = 1.25, 
     col = cols[full_bc_data$sex*2-1 + full_bc_data$river_code-1], pch = c(17, 19)[full_bc_data$sex], main = "Model 3 (Group)")
draws <- as.data.frame(fit_group3)

# - Plot median curve
lines(1:max(age), apply(draws[,grepl("length_pred\\[1,",colnames(draws))], 2, median), col = 1, lty = 1, lwd = 4) # Global

# - Admixed
lines(1:max(age), apply(draws[,grepl("length_pred\\[2,",colnames(draws))], 2, median), col = cols[1], lty = 1, lwd = 2) # Females river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[3,",colnames(draws))], 2, median), col = cols[1], lty = 2, lwd = 2) # Males river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[4,",colnames(draws))], 2, median), col = cols[1], lty = 3, lwd = 2) # Females river 2
lines(1:max(age), apply(draws[,grepl("length_pred\\[5,",colnames(draws))], 2, median), col = cols[1], lty = 4, lwd = 2) # Males river 2

# - Neosho 
lines(1:max(age), apply(draws[,grepl("length_pred\\[6,",colnames(draws))], 2, median), col = cols[2], lty = 1, lwd = 2) # Females river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[7,",colnames(draws))], 2, median), col = cols[2], lty = 2, lwd = 2) # Males river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[8,",colnames(draws))], 2, median), col = cols[2], lty = 3, lwd = 2) # Females river 2
lines(1:max(age), apply(draws[,grepl("length_pred\\[9,",colnames(draws))], 2, median), col = cols[2], lty = 4, lwd = 2) # Males river 2

# - SMB
lines(1:max(age), apply(draws[,grepl("length_pred\\[10,",colnames(draws))], 2, median), col = cols[3], lty = 1, lwd = 2) # Females river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[11,",colnames(draws))], 2, median), col = cols[3], lty = 2, lwd = 2) # Males river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[12,",colnames(draws))], 2, median), col = cols[3], lty = 3, lwd = 2) # Females river 2
lines(1:max(age), apply(draws[,grepl("length_pred\\[13,",colnames(draws))], 2, median), col = cols[3], lty = 4, lwd = 2) # Males river 2

legend("bottomright", c("Admixed", "Neosho","SMB", "Females River 1", "Males River 1", "Females River 2", "Males River 2"), col = c(cols[1:3],1,1,1,1), lty = c(1,1,1,1,2,3,4), bty = "n", lwd = 2)


##### Model 4 #####
# * Final Model ----
# -  VBGF model with sex/river effects and ancestry group and individual level random effects
fit_group4 <- stan(
  file = "growth_analysis/Models/vbgf4_group.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 4000,          # number of warmup iterations per chain
  iter = 7000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = list(max_treedepth = 12, adapt_delta = 0.9)
)

saveRDS(fit_group4, file = "growth_analysis/Models/Fits/vbgf_fit_group4.rds")

# - Print and plot MCMC
print(fit_group4, pars=c("mu_linf", "mu_k", "mu_t0", "beta_linf", "beta_k", "beta_t0", "linf_lineage", "k_lineage", "t0_lineage"), probs=c(.1,.5,.9)) # None of the betas are sig
traceplot(fit_group4, pars = c("mu_linf", "mu_k", "mu_t0"), inc_warmup = FALSE, nrow = 2)
pairs(fit_group4, pars = c("mu_linf", "mu_k", "mu_t0"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit_group4, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

# - Plot fitted model
plot(y = length , x = age, ylab = "Total length (mm)", xlab = "Age (yr)", cex = 2, cex.lab = 1.25, 
     col = cols[full_bc_data$sex*2-1 + full_bc_data$river_code-1], pch = c(17, 19)[full_bc_data$sex], main = "Model 4 (Group)")
draws <- as.data.frame(fit_group4)

# - Plot median curve
lines(1:max(age), apply(draws[,grepl("length_pred\\[1,",colnames(draws))], 2, median), col = 1, lty = 1, lwd = 4) # Global

# - Admixed
lines(1:max(age), apply(draws[,grepl("length_pred\\[2,",colnames(draws))], 2, median), col = cols[1], lty = 1, lwd = 2) # Females river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[3,",colnames(draws))], 2, median), col = cols[1], lty = 2, lwd = 2) # Males river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[4,",colnames(draws))], 2, median), col = cols[1], lty = 3, lwd = 2) # Females river 2
lines(1:max(age), apply(draws[,grepl("length_pred\\[5,",colnames(draws))], 2, median), col = cols[1], lty = 4, lwd = 2) # Males river 2

# - Neosho 
lines(1:max(age), apply(draws[,grepl("length_pred\\[6,",colnames(draws))], 2, median), col = cols[2], lty = 1, lwd = 2) # Females river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[7,",colnames(draws))], 2, median), col = cols[2], lty = 2, lwd = 2) # Males river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[8,",colnames(draws))], 2, median), col = cols[2], lty = 3, lwd = 2) # Females river 2
lines(1:max(age), apply(draws[,grepl("length_pred\\[9,",colnames(draws))], 2, median), col = cols[2], lty = 4, lwd = 2) # Males river 2

# - SMB
lines(1:max(age), apply(draws[,grepl("length_pred\\[10,",colnames(draws))], 2, median), col = cols[3], lty = 1, lwd = 2) # Females river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[11,",colnames(draws))], 2, median), col = cols[3], lty = 2, lwd = 2) # Males river 1
lines(1:max(age), apply(draws[,grepl("length_pred\\[12,",colnames(draws))], 2, median), col = cols[3], lty = 3, lwd = 2) # Females river 2
lines(1:max(age), apply(draws[,grepl("length_pred\\[13,",colnames(draws))], 2, median), col = cols[3], lty = 4, lwd = 2) # Males river 2

legend("bottomright", c("Admixed", "Neosho","SMB", "Females River 1", "Males River 1", "Females River 2", "Males River 2"), col = c(cols[1:3],1,1,1,1), lty = c(1,1,1,1,2,3,4), bty = "n", lwd = 2)


# * Test Lineage Parameters #####
# - No sig difference
par(mfrow = c(3,3))
hist(draws$`linf_lineage[1]`-draws$`linf_lineage[2]`, xlab = "Ad-Neo Linf diff", main = NA)
hist(draws$`k_lineage[1]`-draws$`k_lineage[2]`, xlab = "Ad-Neo K diff", main = "Model 4 Group")
hist(draws$`t0_lineage[1]`-draws$`t0_lineage[2]`, xlab = "Ad-Neo t0 diff", main = NA)


hist(draws$`linf_lineage[1]`-draws$`linf_lineage[3]`, xlab = "Ad-SMB Linf diff", main = NA)
hist(draws$`k_lineage[1]`-draws$`k_lineage[3]`, xlab = "Ad-SMB K diff", main = NA)
hist(draws$`t0_lineage[1]`-draws$`t0_lineage[3]`, xlab = "Ad-SMB t0 diff", main = NA)


hist(draws$`linf_lineage[2]`-draws$`linf_lineage[3]`, xlab = "Neo-SMB Linf diff", main = NA)
hist(draws$`k_lineage[2]`-draws$`k_lineage[3]`, xlab = "Neo-SMB K diff", main = NA)
hist(draws$`t0_lineage[2]`-draws$`t0_lineage[3]`, xlab = "Neo-SMB t0 diff", main = NA)
par(mfrow = c(1,1))



##### Model 5 #####
# -  VBGF model with  ancestry group and individual level random effects
# - No river/sex effects
fit_group5 <- stan(
  file = "growth_analysis/Models/vbgf5_group.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 3000,          # number of warmup iterations per chain
  iter = 5000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = list(max_treedepth = 12, adapt_delta = 0.9)
)

saveRDS(fit_group5, file = "growth_analysis/Models/Fits/vbgf_fit_group5.rds")

# - Print and plot MCMC
print(fit_group5, pars=c("mu_linf", "mu_k", "mu_t0", "linf_lineage", "k_lineage", "t0_lineage"), probs=c(.1,.5,.9)) # None of the betas are sig
traceplot(fit_group5, pars = c("mu_linf", "mu_k", "mu_t0"), inc_warmup = FALSE, nrow = 2)
pairs(fit_group5, pars = c("mu_linf", "mu_k", "mu_t0"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit_group5, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

# - Plot fitted model
plot(y = length , x = age, ylab = "Total length (mm)", xlab = "Age (yr)", cex = 2, cex.lab = 1.25, 
     col = cols[dat$group], pch = 16, main = "Model 5")
draws <- as.data.frame(fit_group5)

# - Plot median curve
lines(1:max(age), apply(draws[,grepl("length_pred\\[1",colnames(draws))], 2, median), col = 1, lty = 1, lwd = 4) # Global

# - Admixed 
lines(1:max(age), apply(draws[,grepl("length_pred\\[2,",colnames(draws))], 2, median), col = cols[2], lty = 1, lwd = 2)

# - Neosho 
lines(1:max(age), apply(draws[,grepl("length_pred\\[3,",colnames(draws))], 2, median), col = cols[2], lty = 1, lwd = 2)

# - SMB
lines(1:max(age), apply(draws[,grepl("length_pred\\[4,",colnames(draws))], 2, median), col = cols[1], lty = 1, lwd = 2)

legend("bottomright", c("Global", "Admixed", "Neosho", "SMB"), col = c(1, cols[1:3]), lty = 1, bty = "n", lwd = 2)


#* Test Lineage Parameters #####
# - No sig difference
par(mfrow = c(3,3))
hist(draws$`linf_lineage[1]`-draws$`linf_lineage[2]`, xlab = "Ad-Neo Linf diff", main = NA)
hist(draws$`k_lineage[1]`-draws$`k_lineage[2]`, xlab = "Ad-Neo K diff", main = "Model 5 Group")
hist(draws$`t0_lineage[1]`-draws$`t0_lineage[2]`, xlab = "Ad-Neo t0 diff", main = NA)


hist(draws$`linf_lineage[1]`-draws$`linf_lineage[3]`, xlab = "Ad-SMB Linf diff", main = NA)
hist(draws$`k_lineage[1]`-draws$`k_lineage[3]`, xlab = "Ad-SMB K diff", main = NA)
hist(draws$`t0_lineage[1]`-draws$`t0_lineage[3]`, xlab = "Ad-SMB t0 diff", main = NA)


hist(draws$`linf_lineage[2]`-draws$`linf_lineage[3]`, xlab = "Neo-SMB Linf diff", main = NA)
hist(draws$`k_lineage[2]`-draws$`k_lineage[3]`, xlab = "Neo-SMB K diff", main = NA)
hist(draws$`t0_lineage[2]`-draws$`t0_lineage[3]`, xlab = "Neo-SMB t0 diff", main = NA)
par(mfrow = c(1,1))


##### Model 6 #####
# -  VBGF model with  ancestry group FIXED and individual level random effects
# - No river/sex effects
fit_group6 <- stan(
  file = "growth_analysis/Models/vbgf6_group.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 3000,          # number of warmup iterations per chain
  iter = 5000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = list(max_treedepth = 12, adapt_delta = 0.9)
)

saveRDS(fit_group6, file = "growth_analysis/Models/Fits/vbgf_fit_group6.rds")

# - Print and plot MCMC
print(fit_group6, pars=c("linf_lineage", "k_lineage", "t0_lineage"), probs=c(.1,.5,.9)) # None of the betas are sig
traceplot(fit_group6, pars = c("linf_lineage", "k_lineage", "t0_lineage"), inc_warmup = FALSE, nrow = 2)
pairs(fit_group6, pars = c("linf_lineage", "k_lineage", "t0_lineage"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit_group6, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

# - Plot fitted model
plot(y = length , x = age, ylab = "Total length (mm)", xlab = "Age (yr)", cex = 2, cex.lab = 1.25, 
     col = cols[dat$group], pch = 16, main = "Model 6")
draws <- as.data.frame(fit_group6)


# - Admixed 
lines(1:max(age), apply(draws[,grepl("length_pred\\[1,",colnames(draws))], 2, median), col = cols[1], lty = 1, lwd = 2)

# - Neosho 
lines(1:max(age), apply(draws[,grepl("length_pred\\[2,",colnames(draws))], 2, median), col = cols[2], lty = 1, lwd = 2)

# - SMB
lines(1:max(age), apply(draws[,grepl("length_pred\\[3,",colnames(draws))], 2, median), col = cols[3], lty = 1, lwd = 2)

legend("bottomright", c("Admixed", "SMB", "Neosho"), col = c(cols[1:3]), lty = 1, bty = "n", lwd = 2)


# * Test Lineage Parameters #####
# - No sig difference
par(mfrow = c(3,3))
hist(draws$`linf_lineage[1]`-draws$`linf_lineage[2]`, xlab = "Ad-Neo Linf diff", main = NA)
hist(draws$`k_lineage[1]`-draws$`k_lineage[2]`, xlab = "Ad-Neo K diff", main = "Model 6 Group")
hist(draws$`t0_lineage[1]`-draws$`t0_lineage[2]`, xlab = "Ad-Neo t0 diff", main = NA)


hist(draws$`linf_lineage[1]`-draws$`linf_lineage[3]`, xlab = "Ad-SMB Linf diff", main = NA)
hist(draws$`k_lineage[1]`-draws$`k_lineage[3]`, xlab = "Ad-SMB K diff", main = NA)
hist(draws$`t0_lineage[1]`-draws$`t0_lineage[3]`, xlab = "Ad-SMB t0 diff", main = NA)


hist(draws$`linf_lineage[2]`-draws$`linf_lineage[3]`, xlab = "Neo-SMB Linf diff", main = NA)
hist(draws$`k_lineage[2]`-draws$`k_lineage[3]`, xlab = "Neo-SMB K diff", main = NA)
hist(draws$`t0_lineage[2]`-draws$`t0_lineage[3]`, xlab = "Neo-SMB t0 diff", main = NA)
par(mfrow = c(1,1))

