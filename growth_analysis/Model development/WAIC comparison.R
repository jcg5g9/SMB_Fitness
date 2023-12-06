library(rstan)
library(loo)

# Back-calculated length-at-age
fit_4 <- readRDS(file = "growth_analysis/Models/Fits/vbgf_fit4.rds")
fit_group4 <- readRDS(file = "growth_analysis/Models/Fits/vbgf_fit_group4.rds")

waic(extract_log_lik(fit_4))
waic(extract_log_lik(fit_group4))

print(fit_4, pars=c( "beta_linf", "beta_k", "beta_t0"), probs=c(.1,.5,.9)) # None of the betas are sig
print(fit_group4, pars=c( "beta_linf", "beta_k", "beta_t0"), probs=c(.05,.5,.95)) # None of the betas are sig

# Sensitivity (true length-at-age)
sen_fit_3 <- readRDS(file = "growth_analysis/Models/Fits/vbgf_sen_fit3.rds")
sen_fit_group3 <- readRDS(file = "growth_analysis/Models/Fits/vbgf_sen_fit_group3.rds")

waic(extract_log_lik(sen_fit_3))
waic(extract_log_lik(sen_fit_group3))

print(sen_fit_3, pars=c( "beta_linf", "beta_k", "beta_t0"), probs=c(.1,.5,.9)) # None of the betas are sig
print(sen_fit_group3, pars=c( "beta_linf", "beta_k", "beta_t0"), probs=c(.1,.5,.9)) # None of the betas are sig
