# Fitting hierarchical growth models using back-calculated length-at-age

##### Setup #####
library(rstan)
library(dplyr)
library(bayesplot)
rstan_options(threads_per_chain = 1)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

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
  cores = 4              # number of cores (could use one per chain)
  ,control = list(adapt_delta=0.999, stepsize=0.001, max_treedepth=18)
)

summ <- summary(fit4, probs=c(.1,.5,.9))$summary

traceplot(fit4, pars = c("mu_linf", "mu_k", "mu_t0", "sigma"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4, pars = c("beta_linf", "beta_k", "beta_t0"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4, pars = c("Lcorr_group", "sigma_group"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4, pars = c("Lcorr_ind", "sigma_ind"), inc_warmup = FALSE, nrow = 2)

rhats <- rhat(fit4)
rhats <- data.frame(Parm = names(rhats), Rhat = rhats)

saveRDS(fit4, file = "growth_analysis/Models/Fits/vbgf_fit4.rds")


# * Convergence diagnostics ----
# https://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html
posterior_cp <- as.array(fit4)

# - Rhat for chain convergence
rhats <- rhat(fit4)
print(rhats)
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats)

# # - Posterior ll
# lp_cp <- log_posterior(fit4)
# head(lp_cp)
#
# # - Accepted prob
# np_cp <- nuts_params(fit4)
# head(np_cp)
#
# # - Look at par vs divergence
# color_scheme_set("darkgray")
# mcmc_parcoord(posterior_cp, np = np_cp)
#
# # - Look at ll vs acceptance
# color_scheme_set("red")
# mcmc_nuts_divergence(np_cp, lp_cp)

# * Print and plot MCMC ----
print(fit4, pars=c("mu_linf", "mu_k", "mu_t0", "beta_linf", "beta_k", "beta_t0", "linf_lineage", "k_lineage", "t0_lineage"), probs=c(.1,.5,.9)) # None of the betas are sig
traceplot(fit4, pars = c("mu_linf", "mu_k", "mu_t0", "sigma"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4, pars = c("beta_linf", "beta_k", "beta_t0"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4, pars = c("Lcorr_group", "sigma_group"), inc_warmup = FALSE, nrow = 2)
traceplot(fit4, pars = c("Lcorr_ind", "sigma_ind"), inc_warmup = FALSE, nrow = 2)


pairs(fit4, pars = c("mu_linf", "mu_k", "mu_t0"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit4, inc_warmup = FALSE)
summary(do.call(rbind, sampler_params), digits = 2)

# - Plot fitted model
cols <- c("#86BBD8","#2F4858", "#F6AE2D", "#F26419") # Colors for the VBGF lines (Females/Males)
plot(y = length , x = age, ylab = "Total length (mm)", xlab = "Age (yr)", cex = 2, cex.lab = 1.25,
     col = cols[full_bc_data$sex*2-1 + full_bc_data$river_code-1], pch = c(17, 19)[full_bc_data$sex], main = "Model 4")
draws <- as.data.frame(fit4)

# - Plot median curve
lines(1:max(age), apply(draws[,grepl("length_pred\\[1,",colnames(draws))], 2, median), col = 1, lty = 1, lwd = 4) # Global

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

