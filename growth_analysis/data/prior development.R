# Code to get the lognormal prior from estimates in:
# Starks, T. A., & Rodger, A. W. (2020). Otolith and scale‐based growth standards for lotic Smallmouth Bass. North American Journal of Fisheries Management, 40(4), 986-994.

# Get mean and sd ----
# - Linf 578 (555-605) # One of these is mis-reported
mulinf <- 578
# - 578 + 1.92 * X = 605
sdlinf <- (605 - 578)/1.92
# - 578 - 1.92 * X = 555
(555 - 578)/-1.92

# - K 0.125 (0.114-0.136)
muk <- 0.125
sdk <- (0.136 - 0.125)/1.92
(0.114 - 0.125)/-1.92

# - T0 -1.67 (-1.79- -1.55)
mut0 <- -1.67
sdt0 <- (1.79 - 1.67)/1.92
(1.55 - 1.67)/-1.92


# Get CV -----
cvlinf <- sdlinf/mulinf
cvk <- sdk/muk

# Convert CV to logSD ----
# https://stats.stackexchange.com/questions/366470/coefficient-of-variation-cv-of-log-transformed-data
# CV = sqrt(exp(logSD^2)-1)
logSDlinf <- sqrt(log((cvlinf^2)+1))
logSDk <- sqrt(log((cvk^2)+1))


curve(dcauchy(x, 0, 0.5), from = 0, to = 1)
