##### READ AND PREPARE DATA #####
library(matrixStats)
# Load data
load('growth_analysis/data/bc_data/full_bc_data.rda')

# - Define data
age <- full_bc_data$annulus
length <- full_bc_data$bc_tl
group <- as.numeric(as.factor(full_bc_data$ancestry_group)) 

cols <- c("#007FFF","#FF7F00") # Colors for the VBGF lines
cols.points <- c("#80bfff", "#ffc180")


##### Get Linf, K, and t0 ##### 
draws <- as.data.frame(model)


##### Plot #####
tiff(file="growth_analysis/Figures/VBGF_curves.tiff" , height= 121.4, pointsize=11,  width=85 , res=300  , units = "mm", family = "serif")
par(mar=c(0.25 , 0.25 , 0 , 0) +0.1, tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0, 0 , 0), cex = 0.75)

# 1) Plot data
plot(y = length , x = age, col = cols.points[group], ylab = "Total length (mm)", xlab = "Age (yr)", cex = .5, cex.lab = 1.25, pch = 20)

# 2) Hyper parameters
mu_linf <- draws[,"mu_linf"] # Subset the mcmc chains
mu_k <- draws[,"mu_k"] # Subset the mcmc chains
mu_t0 <- draws[,"mu_t0"] # Subset the mcmc chains

# - Predict overall curve
ages <- seq(from = -5, to = 50, by = .2)
lengths.mat <- matrix(NA, ncol = length(ages), nrow = length(t0.sub))

for (m in 1: length(ages)){
  lengths.mat[,m] <- mu_linf * (1 - exp(-mu_k * (ages[m] - mu_t0)))
}

# - Plot median curve
median.lengths <- apply(lengths.mat, 2, median)
lines(ages, median.lengths, col = 1, lty = 1, lwd = 2)

# 3) Plot neosho vs smb
for (i in 1:2){ 
  # - Get parameters
  linf_group <- draws[,paste0("linf_group_",i)] # Subset the mcmc chains
  k_group <- draws[,paste0("k_group_",i)] # Subset the mcmc chains
  t0_group <- draws[,paste0("t0_group_",i)] # Subset the mcmc chains
  
  # - Predict curve
  for (m in 1: length(ages)){
    lengths.mat[,m] <- linf_group * (1 - exp(-k_group * (ages[m] - t0_group)))
  }
  
  # - Plot median curve
  median.lengths <- apply(lengths.mat, 2, median)
  lines(ages, median.lengths, col = cols[i], lty = 1, lwd = 2)
}

dev.off()