# Sensitivity analysis fitting hierarchical growth models using observed length-at-age (not back-calculated)

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
warmup = 5000       # number of warmup iterations per chain
thin = 2
iter = 5000        # final number of iterations per chain

##### Assign data to list ##### 
# Real data
load('growth_analysis/data/bc_data/full_bc_data.rda')
full_bc_data <- full_bc_data %>%
  mutate(river_code = factor(river_code),
         river_code = ifelse(river_code == 3, 2, river_code),
         sex = as.numeric(factor(sex)))

# - Set up model matrix (sex/river) and data
model_mat <- full_bc_data %>%
  group_by(sample_id) %>%
  slice(n()) %>%
  select(sample_id, smb, river_code, sex, consensus_age, tl_dead, ancestry_group) %>%
  mutate(
    river_code = factor(river_code),
    sex = factor(sex)) %>%
  arrange(sample_id)

# - Define data
age <- model_mat$consensus_age
length <- model_mat$tl_dead
group <- as.numeric(as.factor(model_mat$ancestry_group)) 
sample.id <- as.numeric(as.factor(as.character(model_mat$sample_id)))
Nind <- length(unique(sample.id))


# - Assign to list
dat = list(
  Nobs = length(length),
  Nages = round(max(age)),
  length = length,
  age = age,
  Zero = rep(0, 3),
  
  eta_scale_prior = 0.5, 
  cholesky_prior = 3,
  beta_scale = 1, 
  
  Nind = Nind,
  Ncoef = 2,
  X = model.matrix(~ sex + river_code, model_mat)[,2:3], # No intercept
  Xhat = matrix(c(0,0,1,0,0,1,1,1), nrow = 4, ncol = 2, byrow = TRUE), # Matrix for prediction
  q = model_mat$smb,
  
  id = sample.id
)


##### FIT MODEL #####
# * Final model ----
# -  VBGF model with sex/river effects and ancestry level random effects
fit_sensitivity <- stan(
  file = "growth_analysis/vbgf_sensitivity.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = warmup,       # number of warmup iterations per chain
  thin = thin,
  iter = iter * thin + warmup,     # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = control
)

saveRDS(fit_sensitivity, file = "growth_analysis/vbgf_fit_sensitivity.rds")
# fit_sensitivity <- readRDS(file = "growth_analysis/vbgf_fit_sensitivity.rds")


# * Convergence diagnostics ----
# https://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html

# # - Posterior ll
# lp_cp <- log_posterior(fit_sensitivity)
# head(lp_cp)
#
# # - Accepted prob
# np_cp <- nuts_params(fit_sensitivity)
# head(np_cp)
#
# # - Look at par vs divergence
# color_scheme_set("darkgray")
# posterior_cp <- as.array(fit_sensitivity)
# mcmc_parcoord(posterior_cp, np = np_cp)
#
# # - Look at ll vs acceptance
# color_scheme_set("red")
# mcmc_nuts_divergence(np_cp, lp_cp)

# * Print and plot MCMC ----
print(fit_sensitivity, pars=c("mu_linf", "mu_k", "mu_t0", "beta_linf", "beta_k", "beta_t0", "linf_lineage", "k_lineage", "t0_lineage"), probs=c(.1,.5,.9)) # None of the betas are sig
traceplot(fit_sensitivity, pars = c("mu_linf", "mu_k", "mu_t0", "sigma"), inc_warmup = FALSE, nrow = 2)
traceplot(fit_sensitivity, pars = c("beta_linf", "beta_k", "beta_t0"), inc_warmup = FALSE, nrow = 2)
traceplot(fit_sensitivity, pars = c("Lcorr", "sigma_group"), inc_warmup = FALSE, nrow = 2)

pairs(fit_sensitivity, pars = c("mu_linf", "mu_k", "mu_t0", "lp__"), las = 1)
pairs(fit_sensitivity, pars = c("linf_lineage", "k_lineage", "t0_lineage", "lp__"), las = 1)
pairs(fit_sensitivity, pars = c("beta_linf", "beta_k", "beta_t0", "lp__"), las = 1)
pairs(fit_sensitivity, pars = c("sigma_group", "Lcorr", "lp__"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit_sensitivity, inc_warmup = FALSE)
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
# pairs_stan(1, fit_sensitivity, pars = c("mu_linf", "mu_k", "mu_t0", "lp__"))
# pairs_stan(1, fit_sensitivity, pars = c("linf_lineage", "k_lineage", "t0_lineage", "lp__"))
# pairs_stan(1, fit_sensitivity, pars = c("beta_linf", "beta_k", "beta_t0", "lp__"))
# pairs_stan(1, fit_sensitivity, pars = c("Lcorr", "lp__"))

check_energy(fit_sensitivity)

# - Rhat for chain convergence
rhats <- rhat(fit_sensitivity)
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats)
rhats <- data.frame(Parm = names(rhats), Rhat = rhats)


# * Plot fitted model ----
cols <- c("#86BBD8","#2F4858", "#F6AE2D", "#F26419") # Colors for the VBGF lines (Females/Males)
draws <- as.data.frame(fit_sensitivity)

png( file = "growth_analysis/Figure2_sensitivity.png" , width=5, height = 5, family = "serif", units = "in", res = 300)
par( mar=c(3, 3 , 0.5 , 3) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))

plot(y = dat$length , x = dat$age, ylab = "Total length (mm)", xlab = "Age (yr)", cex = 2, cex.lab = 1.25, 
     col = cols[as.numeric(model_mat$sex)*2-1 + as.numeric(model_mat$river_code)-1], pch = c(17, 19)[as.numeric(model_mat$sex)])

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
dev.off()


# * Test Lineage Parameters #####
par(mfrow = c(2,3))
# - Posterior
hist(draws$`linf_lineage[1]`-draws$`linf_lineage[2]`, xlab = "Linf diff", main = NA)
hist(draws$`k_lineage[1]`-draws$`k_lineage[2]`, xlab = "K diff", main = "Observed TL")
hist(draws$`t0_lineage[1]`-draws$`t0_lineage[2]`, xlab = "t0 diff", main = NA)

# - Prior
hist(draws$`prior_linf_lineage[1]`-draws$`prior_linf_lineage[2]`, xlab = "Linf diff", main = NA, 
     xlim = range(draws$`linf_lineage[1]`-draws$`linf_lineage[2]`))
hist(draws$`prior_k_lineage[1]`-draws$`prior_k_lineage[2]`, xlab = "K diff",
     xlim = range(draws$`k_lineage[1]`-draws$`k_lineage[2]`))
hist(draws$`prior_t0_lineage[1]`-draws$`prior_t0_lineage[2]`, xlab = "t0 diff", main = NA,
     xlim = range(draws$`t0_lineage[1]`-draws$`t0_lineage[2]`))
par(mfrow = c(1,1))


# * Summary ----
summ <- summary(fit_sensitivity, probs=c(.1,.5,.9))$summary
summ <- summary(fit_sensitivity, pars=c("mu_linf", "mu_k", "mu_t0", "linf_lineage", "post_linf_diff", "k_lineage", "post_k_diff", "t0_lineage", "post_t0_diff", "beta_linf", "beta_k", "beta_t0"), probs=c(.0125, 0.05, .5, 0.95, .975))$summary # None of the betas are sig
summ <- as.data.frame(summ)
summ$parameter <- rownames(summ)
summ <- summ %>%
  select(parameter, mean, `1.25%`, `5%`, `50%`, `95%`, `97.5%`)
write.csv(summ, file = "growth_analysis/Table2_sensitivity_param.csv")


# - SMB - NB Females River 1 difference
length_at_age_diff <- draws[,grepl("pred_length\\[2,",colnames(draws))] - draws[,grepl("pred_length\\[6,",colnames(draws))]
colnames(length_at_age_diff) <- paste0("Difference at age ", 1:8)
length_at_age_diff <- rbind(
  apply(length_at_age_diff, 2, mean), 
  apply(length_at_age_diff, 2, function(x) quantile(x, probs=c(.0125, 0.05, .5, 0.95, .975)))
)
length_at_age_diff <- t(length_at_age_diff)

# - Predicted length-at-age
length_at_age_smb <- draws[,grepl("pred_length\\[2,",colnames(draws))]
colnames(length_at_age_smb) <- paste0("SMB age ", 1:8)
length_at_age_smb <- rbind(
  apply(length_at_age_smb, 2, mean), 
  apply(length_at_age_smb, 2, function(x) quantile(x, probs=c(.0125, 0.05, .5, 0.95, .975)))
)
length_at_age_smb <- t(length_at_age_smb)

length_at_age_n <- draws[,grepl("pred_length\\[6,",colnames(draws))]
colnames(length_at_age_n) <- paste0("NB age ", 1:8)
length_at_age_n <- rbind(
  apply(length_at_age_n, 2, mean), 
  apply(length_at_age_n, 2, function(x) quantile(x, probs=c(.0125, 0.05, .5, 0.95, .975)))
)
length_at_age_n <- t(length_at_age_n)

write.csv(rbind(length_at_age_smb, length_at_age_n, length_at_age_diff), file = "growth_analysis/Table3_sensitivity_length_at_age.csv")

