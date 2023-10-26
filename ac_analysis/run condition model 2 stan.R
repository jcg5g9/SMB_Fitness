# Ancestry-Condition correlation
# SMB Fitness paper
# SJC
# Oct 2023

library(rstan)

# check data
# load condition data
head(condition_data)
cd <- condition_data

# Basic model
# condition ~ smb + river + sex

con <- cd$condition *100
smb <- as.numeric(scale(cd$smb))
n <- length(c)
riv <- as.numeric(as.factor(cd$river))
riv[which(riv==3)] <- 2 # recode bc 1 and 3 was causing issues with jags
sex <- as.numeric(as.factor(cd$sex))
# make riv and sex 0 or 1 for stan
riv1 <- riv-1
sex1 <- sex-1


dat <- list(con=con, smb=smb, N=n, riv=riv1, sex=sex1)

fit2 <- stan(
  file = "ac_analysis/condition_model_2.stan",  
  data = dat,    
  chains = 4,             
  warmup = 2000,          
  iter = 4000,            
  cores = 4,              
)

print(fit2,digits=3)

rstan::traceplot(fit2)
#save (fit2, file='ac_analysis/condition_model_2_out.rds')
load('ac_analysis/condition_model_2_out.rds')
#
sampler_params <- get_sampler_params(fit2, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

draws <- as.data.frame(fit2)

# boxplot of effect sizes
# really just to look at (keeping in mind smb is scaled, and river and sex were binary)
boxplot(draws$beta1, draws$beta2, draws$beta3, xaxt='n', cex.axis=1.5)
mtext('% SMB', side=1, at=1, line=1, cex=1.5)
mtext('River\n (Elk)', side=1, at=2, line=2, cex=1.5)
mtext('Sex\n (Male)', side=1, at=3, line=2, cex=1.5)
abline(h=0, lty='dashed', lwd=2, col='red')


# prediction plot - smb
par(mar=c(4,4,1,1))

b1 <- draws$beta1
b2 <- draws$beta2
b3 <- draws$beta3
alpha <- draws$alpha
min(smb)
max(smb)
x <- seq(-1.3, 2.3, 0.1)
x_us <- ((x*sd(smb))+mean(smb))

pred <- matrix(NA, length(b1), length(x)) # empty matrix for prediction plot data
for (i in 1:length(x)){ # fill empty matrix with predicted values using the model
  pred[,i] <- alpha + b1*x[i] + b2*1 + b3*1
}

meanb1 <- upperb1 <- lowerb1 <- numeric()
for (i in 1:length(x)){
  meanb1[i] <- mean(pred[,i]) # mean at each x
  upperb1[i] <- quantile(pred[,i], probs=0.975) # 97.5% quantile (top of 95% CRI) at each x
  lowerb1[i] <- quantile(pred[,i], probs=0.025)
}

blue1 <- rgb(10,40,100,85,maxColorValue=255)

plot(x=smb, y=con, type='p',xlab='', ylab='', col='gray80',
     pch=20, cex.axis=1.3, main="", yaxt='n')
axis(side=2, at=seq(0.8, 1.8, 0.1), las=2, cex.axis=1.3)
polygon(x=c(x_us,rev(x_us)), y=c(lowerb1,rev(upperb1)), col=blue1, border=NA) # polygon for 95% CRI 
lines(meanb1 ~ x_us, type='l', col='deepskyblue4', lty='solid', pch=16, cex=2, lwd=2.5) # line for mean
mtext("Body Condition Index", side=2, line=2.5, cex=1.4)
mtext("% SMB", side=1, line=2.5, cex=1.4)

