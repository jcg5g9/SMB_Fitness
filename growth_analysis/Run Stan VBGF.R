library(rstan)
library(dplyr)

##### Assign data to list ##### 
# Real data
load('growth_analysis/data/bc_data/full_bc_data.rda')

# - Define data
age <- full_bc_data$annulus
length <- full_bc_data$bc_tl
group <- as.numeric(as.factor(full_bc_data$ancestry_group)) 
sample.id <- as.numeric(as.factor(full_bc_data$sample_id))
Nind <- length(unique(full_bc_data$sample_id))

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
  length = length,
  age = age,
  
  Nind = Nind,
  Ncoef = 2,
  X = model.matrix(~-1 + river_code + sex, model_mat),
  q = model_mat$smb,
  
  id = sample.id
)


##### RUN THE MODEL IN STAN #####
fit1 <- stan(
  file = "growth_analysis/vbgf.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 1,              # number of cores (could use one per chain)
  refresh = 0             # no progress shown
)
