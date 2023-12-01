# Fitting hierarchical growth models using back-calculated length-at-age

##### Setup #####
library(rstan)
library(dplyr)
library(bayesplot)
rstan_options(threads_per_chain = 1)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

control = list(adapt_delta=0.999, stepsize=0.001, max_treedepth=18)
control = list()

##### Assign data to list ##### 
# Real data
load('growth_analysis/data/bc_data/full_bc_data.rda')
full_bc_data <- full_bc_data %>%
  mutate(river_code = factor(river_code),
         river_code = ifelse(river_code == 3, 2, river_code),
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
  select(sample_id, smb, river_code, sex) %>%
  mutate(river_code = factor(river_code),
         sex = factor(sex)) %>%
  arrange(sample_id)

# - Assign to list
dat = list(
  Nobs = length(length),
  Nages = round(max(age)),
  length = length,
  age = age,
  Zero = rep(0, 3),
  
  cauchy_scale = 0.5, 
  cholesky_prior = 3,
  beta_scale = 0.25, 
  
  Nind = Nind,
  Ncoef = 2,
  X = model.matrix(~ sex + river_code, model_mat)[,2:3], # No intercept
  Xhat = matrix(c(0,0,1,0,0,1,1,1), nrow = 4, ncol = 2, byrow = TRUE), # Matrix for prediction
  q = model_mat$smb,
  
  id = sample.id
)


##### Model 4 #####
# * Final Model ----
# -  VBGF model with sex/river effects and ancestry and individual level random effects
fit4 <- stan(
  file = "growth_analysis/Models/vbgf4.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 5000,          # number of warmup iterations per chain
  iter = 10000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = control
)

summ <- summary(fit4, probs=c(.1,.5,.9))$summary

traceplot(fit4, pars = c("mu_linf", "mu_k", "mu_t0", "sigma"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4, pars = c("beta_linf", "beta_k", "beta_t0"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4, pars = c("Lcorr_group", "sigma_group"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4, pars = c("Lcorr_ind", "sigma_ind"), inc_warmup = FALSE, nrow = 2)

rhats <- rhat(fit4)
rhats <- data.frame(Parm = names(rhats), Rhat = rhats)

saveRDS(fit4, file = "growth_analysis/Models/Fits/vbgf_fit4.rds")



# Sex only (no river) ----
# - Assign to list
dat_sex = dat
dat_sex$Ncoef = 1
dat_sex$X = as.matrix(model.matrix(~ sex, model_mat)[,2]) # No intercept
dat_sex$Xhat = matrix(c(0,0,1,1), nrow = 4, ncol = 1, byrow = TRUE) # Matrix for prediction

# -  VBGF model with sex/river effects and ancestry and individual level random effects
fit4_sex <- stan(
  file = "growth_analysis/Models/vbgf4.stan",  # Stan program
  data = dat_sex,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 5000,          # number of warmup iterations per chain
  iter = 10000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = control
)

saveRDS(fit4_sex, file = "growth_analysis/Models/Fits/vbgf_fit4_sex.rds")

summ <- summary(fit4_sex, probs=c(.1,.5,.9))$summary

traceplot(fit4_sex, pars = c("mu_linf", "mu_k", "mu_t0", "sigma"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4_sex, pars = c("beta_linf", "beta_k", "beta_t0"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4_sex, pars = c("Lcorr_group", "sigma_group"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4_sex, pars = c("Lcorr_ind", "sigma_ind"), inc_warmup = FALSE, nrow = 2)

rhats <- rhat(fit4_sex)
rhats <- data.frame(Parm = names(rhats), Rhat = rhats)



# River only (no sex) ----
# - Assign to list
dat_river = dat
dat_river$Ncoef = 1
dat_river$X = as.matrix(model.matrix(~ sex, model_mat)[,2]) # No intercept
dat_river$Xhat = matrix(c(0,0,1,1), nrow = 4, ncol = 1, byrow = TRUE) # Matrix for prediction


# -  VBGF model with river effects and ancestry and individual level random effects
fit4_river <- stan(
  file = "growth_analysis/Models/vbgf4.stan",  # Stan program
  data = dat_river,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 5000,          # number of warmup iterations per chain
  iter = 15000,            # total number of iterations per chain
  cores = 4,             # number of cores (could use one per chain)
  control = control
)

saveRDS(fit4_river, file = "growth_analysis/Models/Fits/vbgf_fit4_river.rds")

summ <- summary(fit4_river, probs=c(.1,.5,.9))$summary

traceplot(fit4_river, pars = c("mu_linf", "mu_k", "mu_t0", "sigma"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4_river, pars = c("beta_linf", "beta_k", "beta_t0"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4_river, pars = c("Lcorr_group", "sigma_group"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4_river, pars = c("Lcorr_ind", "sigma_ind"), inc_warmup = FALSE, nrow = 2)

rhats <- rhat(fit4_river)
rhats <- data.frame(Parm = names(rhats), Rhat = rhats)



##### Model 5 #####
# -  VBGF model with  ancestry and individual level random effects
# - No river/sex effects
fit5 <- stan(
  file = "growth_analysis/Models/vbgf5.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 6000,          # number of warmup iterations per chain
  iter = 9000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = control
)

saveRDS(fit5, file = "growth_analysis/Models/Fits/vbgf_fit5.rds")

summ <- summary(fit5, probs=c(.1,.5,.9))$summary

traceplot(fit5, pars = c("mu_linf", "mu_k", "mu_t0", "sigma"), inc_warmup = FALSE, nrow = 2)
traceplot(fit5, pars = c("Lcorr_group", "sigma_group"), inc_warmup = FALSE, nrow = 2)
traceplot(fit5, pars = c("Lcorr_ind", "sigma_ind"), inc_warmup = FALSE, nrow = 2)

rhats <- rhat(fit5)
rhats <- data.frame(Parm = names(rhats), Rhat = rhats)
print(rhats)




