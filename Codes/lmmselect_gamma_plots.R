## Sample analysis
## Adult income data from UCI repository
rm(list=ls())
#setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/Rain")
setwd('D:/Study/My projects/Climate-indian-monsoon/Codes')
source('misc_functions.R')

library(rrcov)
library(fda.usc)
library(lme4)
library(MuMIn)
library(doSNOW)
library(parallel)

# read in data
rainsmall = read.csv("../Data/rainsmall.csv")
long.list = rainsmall$LONGITUDE[1:36]
lat.list = rainsmall$LATITUDE[1:36]

# check full model
rainsmall[-(1:3)] = scale(rainsmall[-(1:3)])
varnames = names(rainsmall)[-(1:3)]
fixed.full = paste(varnames, collapse="+")
random_terms = "+ (1|year)"
form.full = as.formula(paste("log(PRCP+1) ~", fixed.full, random_terms))
mod.full = lmer(form.full, data=rainsmall)
summary(mod.full)
anova(mod.full)
r.squaredGLMM(mod.full)

# set up residuals
y = getME(mod.full, 'y')
fixed = getME(mod.full, 'X') %*% fixef(mod.full)
eta = unlist(ranef(mod.full))
Z = t(as.matrix(getME(mod.full,'Zt')))
random = Z %*% eta
r = y - fixed - random

# bootstrap parameters
n = nrow(rainsmall)
nr = length(eta)
p = ncol(rainsmall)-3

load('bestmodel.rda')
best.index = which.min(lapply(outputs, function(x) x[[3]]))
Cn.frame = outputs[[best.index]][[1]]
mod.final = outputs[[best.index]][[2]]
summary(mod.final)
anova(mod.final)
r.squaredGLMM(mod.final)
anova(mod.final, mod.full)

# future prediction
testyrs = 2003:2012
ntest = length(testyrs)
pred.mat.MSE = matrix(0, ncol=2, nrow=ntest)
pred.mat.bias = pred.mat.MSE

for (i in 1:ntest){
iyr = testyrs[i]
itrain = which(rainsmall$year >= iyr-25 & rainsmall$year < iyr)
itest = which(rainsmall$year == iyr)
testX = rainsmall[itest,-c(1:3)]

mod.full.train = update(mod.full, subset=itrain)
mod.final.train = update(mod.final, subset=itrain)

ytest = log(rainsmall$PRCP[itest]+1)
yhat.full = getME(mod.full,"X")[itest,] %*% as.numeric(fixef(mod.full.train))
diff.full = ytest - yhat.full
pred.mat.bias[i,1] = mean(diff.full)
pred.mat.MSE[i,1] = mean(diff.full^2)

yhat.final = getME(mod.final,"X")[itest,] %*% as.numeric(fixef(mod.final.train))
diff.final = ytest - yhat.final
pred.mat.bias[i,2] = mean(diff.final)
pred.mat.MSE[i,2] = mean(diff.final^2)
}

pdf('rolling_predMSE_full_vs_reduced_gamma.pdf',5,5)
plot(pred.mat.MSE[,1]~testyrs, type="b", ylim=c(0,ceiling(max(pred.mat.MSE[,1]))), lwd=2,
	xlab="year", ylab="MSE")
lines(pred.mat.MSE[,2]~testyrs, type="b", lty=2, lwd=2)
legend("topleft", c("Full model", "Reduced model"), lty=1:2, lwd=2)
dev.off()

pdf('rolling_predbias_full_vs_reduced_gamma.pdf',5,5)
plot(pred.mat.bias[,1]~testyrs, type="b", ylim=c(-3,3), lwd=2,
	xlab="year", ylab="Bias")
lines(pred.mat.bias[,2]~testyrs, type="b", lty=2, lwd=2)
abline(h=0, lty=2)
legend("topleft", c("Full model", "Reduced model"), lty=1:2, lwd=2)
dev.off()

pdf('rolling_density2012_full_vs_reduced_gamma.pdf',5,5)
plot(density(ytest), xlim=c(-2,10), ylim=c(0,.5), lwd=2,
	xlab="log(PRCP+1)", ylab="density", main="Year 2012")
lines(density(yhat.full), col='red', lwd=2)
lines(density(yhat.final), col='blue', lwd=2)
legend("topright", c("Truth", "Full model pred", "Reduced model pred"),
	col=c('black','red','blue'), lty=1, lwd=2)
dev.off()

# plot on map
library(maps)
library(maptools)

r.final = ytest - yhat.final
q.final = quantile(r.final, c(.05,.25,.75,.95))
pdf('rolling_map2012_full_vs_reduced_gamma.pdf',5,5)
defaultPar = par()
par(mar=rep(0,4))
India <- map("world", ylim=c(8,38), xlim=c(68,98))
title("2012")
points(long.list, lat.list, pch=ifelse(r.final<0, 19, 17),
	cex=ifelse(r.final < q.final[1] | r.final > q.final[4], 2,
		ifelse(r.final < q.final[2] | r.final > q.final[3], 1, 0.75)),
	col=ifelse(r.final<0, "black", "red"))
legend('topright', c("Positive resid", "negative resid"),
	pch=c(17,19), col=c("red", "black"), bty="n")
par(defaultPar)
dev.off()

