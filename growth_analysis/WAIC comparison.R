library(rstan)
library(loo)

# Back-calculated length-at-age
fit_4 <- readRDS(file = "growth_analysis/Models/Fits/vbgf_fit4.rds")
fit_group4 <- readRDS(file = "growth_analysis/Models/Fits/vbgf_fit_group4.rds")

waic(extract_log_lik(fit_4))
waic(extract_log_lik(fit_group4))

# Sensitivity (true length-at-age)
sen_fit_3 <- readRDS(file = "growth_analysis/Models/Fits/vbgf_sen_fit3.rds")
sen_fit_group3 <- readRDS(file = "growth_analysis/Models/Fits/vbgf_sen_fit_group3.rds")

waic(extract_log_lik(sen_fit_3))
waic(extract_log_lik(sen_fit_group3))