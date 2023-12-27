# Fitting hierarchical growth models using back-calculated length-at-age

##### Setup #####
library(rstan)
library(dplyr) 
library(bayesplot)
rstan_options(threads_per_chain = 1)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

# - MCMC specifications
control = list(adapt_delta=0.999, stepsize=0.001, max_treedepth=18)
# control = list()
warmup = 6000       # number of warmup iterations per chain
thin = 2
iter = 6000        # final number of iterations per chain

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
  
  eta_scale_prior = 0.5, 
  cholesky_prior = 3,
  beta_scale = 0.5, 
  
  Nind = Nind,
  Ncoef = 2,
  X = model.matrix(~ sex + river_code, model_mat)[,2:3], # No intercept
  Xhat = matrix(c(0,0,1,0,0,1,1,1), nrow = 4, ncol = 2, byrow = TRUE), # Matrix for prediction
  q = model_mat$smb,
  
  id = sample.id
)


##### FIT MODEL #####
# * Final Model ----
# -  VBGF model with sex/river effects and ancestry and individual level random effects
fit_backcalculated <- stan(
  file = "growth_analysis/vbgf_backcalculated.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = warmup,       # number of warmup iterations per chain
  thin = thin,
  iter = iter * thin + warmup,     # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = control
)

saveRDS(fit_backcalculated, file = "growth_analysis/vbgf_fit_backcalculated.rds")
# fit_backcalculated <- readRDS(file = "growth_analysis/vbgf_fit_backcalculated.rds")


# * Convergence diagnostics ----
# https://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html


# # - Posterior ll
# lp_cp <- log_posterior(fit_backcalculated)
# head(lp_cp)
#
# # - Accepted prob
# np_cp <- nuts_params(fit_backcalculated)
# head(np_cp)
#
# # - Look at par vs divergence
# color_scheme_set("darkgray")
# posterior_cp <- as.array(fit_backcalculated)
# mcmc_parcoord(posterior_cp, np = np_cp)
#
# # - Look at ll vs acceptance
# color_scheme_set("red")
# mcmc_nuts_divergence(np_cp, lp_cp)

# * Print and plot MCMC ----
print(fit_backcalculated, pars=c("mu_linf", "mu_k", "mu_t0", "beta_linf", "beta_k", "beta_t0", "linf_lineage", "k_lineage", "t0_lineage"), probs=c(.1,.5,.9)) # None of the betas are sig
traceplot(fit_backcalculated, pars = c("mu_linf", "mu_k", "mu_t0", "sigma"), inc_warmup = FALSE, nrow = 2)
traceplot(fit_backcalculated, pars = c("beta_linf", "beta_k", "beta_t0"), inc_warmup = FALSE, nrow = 2)
traceplot(fit_backcalculated, pars = c("Lcorr", "sigma_group"), inc_warmup = FALSE, nrow = 2)
traceplot(fit_backcalculated, pars = c("Lcorr", "sigma_ind"), inc_warmup = FALSE, nrow = 2)

pairs(fit_backcalculated, pars = c("mu_linf", "mu_k", "mu_t0", "lp__"), las = 1)
pairs(fit_backcalculated, pars = c("linf_lineage", "k_lineage", "t0_lineage", "lp__"), las = 1)
pairs(fit_backcalculated, pars = c("beta_linf", "beta_k", "beta_t0", "lp__"), las = 1)
# pairs(fit_backcalculated, pars = c("sigma_group", "Lcorr", "lp__"), las = 1)
# pairs(fit_backcalculated, pars = c("sigma_ind", "Lcorr", "lp__"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit_backcalculated, inc_warmup = FALSE)
summary(do.call(rbind, sampler_params), digits = 2)

# BFMI low
# pairs_stan <- function(chain, stan_model, pars) {
#   energy <- as.matrix(sapply(get_sampler_params(stan_model, inc_warmup = F), 
#                              function(x) x[,"energy__"]))
#   pars <- extract(stan_model, pars = pars, permuted = F)
#   df <- data.frame(energy[,chain], pars[,chain,])
#   names(df)[1] <- "energy"
#   GGally::ggpairs(df, title = paste0("Chain", chain), 
#                   lower = list(continuous = GGally::wrap("points", alpha = 0.2)))                    
# }
# 
# pairs_stan(1, fit_backcalculated, pars = c("mu_linf", "mu_k", "mu_t0", "lp__"))
# pairs_stan(1, fit_backcalculated, pars = c("linf_lineage", "k_lineage", "t0_lineage", "lp__"))
# pairs_stan(1, fit_backcalculated, pars = c("beta_linf", "beta_k", "beta_t0", "lp__"))
# pairs_stan(1, fit_backcalculated, pars = c("Lcorr", "lp__"))

check_energy(fit_backcalculated)

# - Rhat for chain convergence
rhats <- rhat(fit_backcalculated)
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats)
rhats <- data.frame(Parm = names(rhats), Rhat = rhats)


# * Plot fitted model ----
cols <- c("#86BBD8","#2F4858", "#F6AE2D", "#F26419") # Colors for the VBGF lines (Females/Males)
draws <- as.data.frame(fit_backcalculated)


png( file = "growth_analysis/Figure2_backcalculated.png" , width=5, height = 5, family = "serif", units = "in", res = 300)
par( mar=c(3, 3 , 0.5 , 3) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))

plot(y = length , x = age, ylab = "Back-calculated TL (mm)", xlab = "Annuli", cex = 2, cex.lab = 1.25,
     col = cols[full_bc_data$sex*2-1 + full_bc_data$river_code-1], pch = c(17, 19)[full_bc_data$sex])
# - Plot median curve
lines(1:max(age), apply(draws[,grepl("pred_length\\[1,",colnames(draws))], 2, median), col = 1, lty = 1, lwd = 4) # Global

# - SMB Females
lines(1:max(age), apply(draws[,grepl("pred_length\\[2,",colnames(draws))], 2, median), col = cols[1], lty = 1, lwd = 2) # Females river 1
lines(1:max(age), apply(draws[,grepl("pred_length\\[4,",colnames(draws))], 2, median), col = cols[31], lty = 2, lwd = 2) # Females river 2

# - Neosho Females
lines(1:max(age), apply(draws[,grepl("pred_length\\[6,",colnames(draws))], 2, median), col = cols[2], lty = 1, lwd = 2) # Females river 1
lines(1:max(age), apply(draws[,grepl("pred_length\\[8,",colnames(draws))], 2, median), col = cols[2], lty = 2, lwd = 2) # Females river 2

# - SMB males
lines(1:max(age), apply(draws[,grepl("pred_length\\[3,",colnames(draws))], 2, median), col = cols[3], lty = 1, lwd = 2) # Males river 1
lines(1:max(age), apply(draws[,grepl("pred_length\\[5,",colnames(draws))], 2, median), col = cols[3], lty = 2, lwd = 2) # Males river 2

# - Neosho Males
lines(1:max(age), apply(draws[,grepl("pred_length\\[7,",colnames(draws))], 2, median), col = cols[4], lty = 1, lwd = 2) # Males river 1
lines(1:max(age), apply(draws[,grepl("pred_length\\[9,",colnames(draws))], 2, median), col = cols[4], lty = 2, lwd = 2) # Males river 2

legend("bottomright", c("SMB Females", "Neosho Females","SMB Males", "Neosho Males", "Big Sugar Creek", "Elk River"), col = c(cols,1,1), lty = c(1,1,1,1,1,2), bty = "n", lwd = 2)
dev.off()


# * Test Lineage Parameters #####
# - No sig difference
par(mfrow = c(2,3))
# - Posterior
hist(draws$`linf_lineage[1]`-draws$`linf_lineage[2]`, xlab = "Linf diff", main = NA)
hist(draws$`k_lineage[1]`-draws$`k_lineage[2]`, xlab = "K diff", main = "Back-calculated TL")
hist(draws$`t0_lineage[1]`-draws$`t0_lineage[2]`, xlab = "t0 diff", main = NA)

# - Prior
hist(draws$`prior_linf_lineage[1]`-draws$`prior_linf_lineage[2]`, xlab = "Linf diff", main = NA)
hist(draws$`prior_k_lineage[1]`-draws$`prior_k_lineage[2]`, xlab = "K diff")
hist(draws$`prior_t0_lineage[1]`-draws$`prior_t0_lineage[2]`, xlab = "t0 diff", main = NA)
par(mfrow = c(1,1))


# * Posterior prior plots ----
source("growth_analysis/Plot densities.R")
plot_density(fit_backcalculated, file_name = "growth_analysis/Figure3_backcalculated", probs = c(0.01, 0.99))


# * Summary ----
summ <- summary(fit_backcalculated, probs=c(.1,.5,.9))$summary
summ <- summary(fit_backcalculated, pars=c("mu_linf", "mu_k", "mu_t0", "linf_lineage", "post_linf_diff", "k_lineage", "post_k_diff", "t0_lineage", "post_t0_diff", "beta_linf", "beta_k", "beta_t0"), probs=c(.025, 0.05, .5, 0.95, .975))$summary # None of the betas are sig
summ <- as.data.frame(summ)
summ$parameter <- rownames(summ)
summ <- summ %>%
  select(parameter, mean, `2.5%`, `5%`, `50%`, `95%`, `97.5%`)
write.csv(summ, file = "growth_analysis/Table2_backcalculated_param.csv")


# - SMB - NB Females River 1 difference
length_at_age_diff <- draws[,grepl("pred_length\\[2,",colnames(draws))] - draws[,grepl("pred_length\\[6,",colnames(draws))]
colnames(length_at_age_diff) <- paste0("Difference at ", 1:8, " annuli")
length_at_age_diff <- rbind(
  apply(length_at_age_diff, 2, mean), 
  apply(length_at_age_diff, 2, function(x) quantile(x, probs=c(.025, 0.05, .5, 0.95, .975)))
)
length_at_age_diff <- t(length_at_age_diff)

# - Predicted length-at-age
length_at_age_smb <- draws[,grepl("pred_length\\[2,",colnames(draws))]
colnames(length_at_age_smb) <- paste0("SMB ",  1:8, " annuli")
length_at_age_smb <- rbind(
  apply(length_at_age_smb, 2, mean), 
  apply(length_at_age_smb, 2, function(x) quantile(x, probs=c(.025, 0.05, .5, 0.95, .975)))
)
length_at_age_smb <- t(length_at_age_smb)

length_at_age_n <- draws[,grepl("pred_length\\[6,",colnames(draws))]
colnames(length_at_age_n) <- paste0("NB ", 1:8, " annuli")
length_at_age_n <- rbind(
  apply(length_at_age_n, 2, mean), 
  apply(length_at_age_n, 2, function(x) quantile(x, probs=c(.025, 0.05, .5, 0.95, .975)))
)
length_at_age_n <- t(length_at_age_n)

write.csv(rbind(length_at_age_smb, length_at_age_n, length_at_age_diff), file = "growth_analysis/Table3_backcalculated_length_at_age.csv")


