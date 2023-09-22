##### Setup #####
library(rstan)
library(dplyr)
rstan_options(threads_per_chain = 1)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

##### Assign data to list ##### 
# Real data
load('growth_analysis/data/bc_data/full_bc_data.rda')

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
  mutate(river_code = ifelse(river_code == 3, 2, river_code),
         river_code = factor(river_code),
         sex = factor(sex)) %>%
  arrange(sample_id)

# - Assign to list
dat = list(
  Nobs = length(length),
  Nages = round(max(age)),
  length = length,
  age = age,
  
  Nind = Nind,
  Ncoef = 2,
  X = model.matrix(~ sex + river_code, model_mat)[,2:3], # No intercept
  q = model_mat$smb,
  
  id = sample.id
)


##### RUN THE MODEL IN STAN #####
fit1 <- stan(
  file = "growth_analysis/vbgf.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 4000,          # number of warmup iterations per chain
  iter = 6000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = list(max_treedepth = 14, adapt_delta = 0.9)
)

saveRDS(fit1, file = "growth_analysis/vbgf_fit.rds")

##### PRINT AND PLOT #####
print(fit1, pars=c("mu_linf", "mu_k", "mu_t0", "linf_group", "k_group", "t0_group"), probs=c(.1,.5,.9))
traceplot(fit1, pars = c("mu_linf", "mu_k", "mu_t0", "linf_group", "k_group", "t0_group"), inc_warmup = FALSE, nrow = 2)
pairs(fit1, pars = c("mu_linf", "mu_k", "mu_t0"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit1, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)
