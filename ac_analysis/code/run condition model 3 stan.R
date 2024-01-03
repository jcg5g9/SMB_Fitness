# Ancestry-Condition correlation
# Run Condition Model 3 (with interactions)
# SMB Fitness paper
# SJC
# DEC 2023

library(rstan)

load("~/GitHub/SMB_Fitness/ac_analysis/data/condition_data.Rda")

# check data
# load condition data
head(condition_data)
cd <- condition_data

# Basic model
# condition ~ smb + river + sex

con <- cd$condition *100
smb <- as.numeric(scale(cd$smb))
n <- length(con)
riv <- as.numeric(as.factor(cd$river))
riv[which(riv==3)] <- 2 # recode to 1 and 2 instead of 1 and 3
sex <- as.numeric(as.factor(cd$sex))
# make riv and sex 0 or 1 for stan
riv1 <- riv-1
sex1 <- sex-1


dat <- list(con=con, smb=smb, N=n, riv=riv1, sex=sex1)

fit2 <- stan(
  file = "ac_analysis/condition_model_3_TEST.stan",  
  data = dat,    
  chains = 4,             
  warmup = 2000,          
  iter = 4000,            
  cores = 4,              
)



rstan::traceplot(fit2)
#save (fit2, file='ac_analysis/condition_model_3_out.rds')
load('ac_analysis/condition_model_3_out.rds')

print(fit2,digits=3)

sampler_params <- get_sampler_params(fit2, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)
draws <- as.data.frame(fit2)

# boxplot of effect sizes
# really just to look at (keeping in mind smb is scaled, and river and sex were binary)
# make sure to put in caption that big sugar and female are both represented by the dotted red line
boxplot(draws$beta1, draws$beta2, draws$beta3, draws$beta4, cex.axis=1.5, outline=F, xaxt="n",
        col='gray50',border='gray50', las=2
        ,boxwex=0.3, whisklty=1, staplelwd=0, staplewex=0, whisklwd=4, boxlwd=2, medcol="white")
mtext('% SMB', side=1, at=1, line=1, cex=1.5)
mtext('Stream\n (Elk)', side=1, at=2, line=2, cex=1.5)
mtext('Sex\n (Male)', side=1, at=3, line=2, cex=1.5)
mtext('SMB x\n Stream', side=1, at=4, line=2, cex=1.5)
abline(h=0, lty='dashed', lwd=2, col='black')
mtext('Effect Sizes', side=2, line=4, cex=1.5)


# prediction plot - smb
par(mar=c(4,6,1,1))

b1 <- draws$beta1
b2 <- draws$beta2
b3 <- draws$beta3
b41 <- draws$beta4
b42

# find p 
length(which(b1<0))/length(b1) # 0.963
length(which(b2>0))/length(b2) # 0.986
length(which(b3<0))/length(b3) # 0.727
length(which(b4<0))/length(b4) # 0.595

alpha <- draws$alpha
min(smb)
max(smb)
x <- seq(-1.3, 2.3, 0.1)
x_us <- x*sd(cd$smb)+mean(cd$smb)

pred <- matrix(NA, length(b1), length(x)) # empty matrix for prediction plot data
for (i in 1:length(x)){ # fill empty matrix with predicted values using the model
  pred[,i] <- alpha + b1*x[i] + b2*0 + b3*0 + b4*0
}

meanb1 <- upperb1 <- lowerb1 <- numeric()
for (i in 1:length(x)){
  meanb1[i] <- mean(pred[,i])/1000 # mean at each x
  upperb1[i] <- quantile(pred[,i], probs=0.975)/1000 # 97.5% quantile (top of 95% CRI) at each x
  lowerb1[i] <- quantile(pred[,i], probs=0.025)/1000
}

blue1 <- rgb(10, 40, 100, 85,maxColorValue=255)
red1 <- rgb(100, 40, 10, 85, maxColorValue=255)

plot(x=smb*sd(cd$smb)+mean(cd$smb), y=con/1000, type='p',xlab='', ylab='',
     pch=16, cex=1.3, col='black', cex.axis=1.3, main="", yaxt='n', xaxt='n')
points(x=smb*sd(cd$smb)+mean(cd$smb), y=con/1000, pch=16, cex=1, col='gray80')
axis(side=2, at=seq(0.8, 1.8, 0.1)/1000, las=2, cex.axis=1.2)
axis(side=1, at=seq(0, 1, 0.1), cex.axis=1.3)
polygon(x=c(x_us,rev(x_us)), y=c(lowerb1,rev(upperb1)), col=blue1, border=NA) # polygon for 95% CRI 
lines(meanb1 ~ x_us, type='l', col='deepskyblue4', lty='solid', pch=16, cex=2, lwd=2.5) # line for mean
mtext(expression(paste("Body Condition Index (g/mm"^"3",")")), side=2, line=4.2, cex=1.4)
mtext("% SMB", side=1, line=2.5, cex=1.4)
