## Sample analysis
## Adult income data from UCI repository
rm(list=ls())
#setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/Rain")
setwd('C:/Study/My projects/Climate-indian monsoon/Codes')
source('misc_functions.R')

library(rrcov)
library(fda.usc)
library(lme4)
library(MuMIn)
library(doSNOW)
library(parallel)

# read in data
rainsmall = read.csv("rainsmall.txt")

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
sdn = n^(.05)
sdr = nr^(.05)

## loop to get drop 1 bootstrap estimates
loopfun = function(i){
  require(lme4)
  require(fda.usc)
  
  iyb = fixed + (Z %*% (rnorm(nr)*sdr*eta) + rnorm(n)*sdn*r)/sqrt(1 - p/n)
    
	SSPvec = rep(0,p)
  # loop for all variables
	for(j in 1:p){
		# build model
		jformula = paste(varnames[-j], collapse="+")
		jformula = as.formula(paste("iyb ~", jformula, random_terms))
		ijmod = lmer(jformula, data=rainsmall)

		# extract coef and extend to lengthp by appending 0's
		ibeta.j = as.numeric(coef(ijmod)$year[1,-1])
    		ibeta.j.full = rep(0,p)
		ibeta.j.full[-j] = ibeta.j

    # get depth-based criterion
		SSPvec[j] = mdepth.RP(ibeta.j.full, beta.mat)$dep
  	}
  SSPvec
}

# loop for full model bootstrap
# to get bootstrap distribution of parameters
nboot = 1e3
SSPmat = matrix(0, nrow=nboot, ncol=p+1)
SSPmat.d = SSPmat
beta.mat = matrix(0, nrow=nboot, ncol=p)
set.seed(07152015)

# create progress bar
pb <- txtProgressBar(min = 0, max = nboot, style = 3)

for(i in 1:nboot){
  iyb = fixed + (Z %*% (rnorm(nr)*sdr*eta) + rnorm(n)*sdn*r)/sqrt(1 - p/n)
  iformula = paste(varnames, collapse="+")
  iformula = as.formula(paste("iyb ~", iformula, random_terms))
  imod = lmer(iformula, data=rainsmall)
  beta.mat[i,] = as.numeric(coef(imod)$year[1,-1])
  setTxtProgressBar(pb, i)
}
close(pb)
SSPmat.d[,p+1] = mdepth.RP(beta.mat, beta.mat)$dep

# run function in parallel
cl <- makeCluster(detectCores()-2)
registerDoSNOW(cl)
system.time(SSPtab <- foreach(i=1:nboot) %dopar% loopfun(i))
SSPmat.d[,1:p] = matrix(unlist(SSPtab), ncol=p, byrow=T)
stopCluster(cl)

# get p-values
pVal = rep(1, p+1)
for(i in 1:p){
  pVal[i] = t.test(SSPmat.d[,i], SSPmat.d[,p+1], alternative="less")$p.value
}

Cn.frame = data.frame(DroppedVar = c(paste("-", names(data.frame(getME(mod.full, 'X')))[-1]), "<none>"),
                      Cn = apply(SSPmat.d, 2, mean),
                      pValue = pVal)
Cn.frame = Cn.frame[with(Cn.frame, order(Cn)),]
row.names(Cn.frame) = NULL
Cn.frame


write.csv(Cn.frame, 'wild_ntothepoint1.csv')

# parameter distributions
pairs(beta.mat[,1:5], pch=19, cex=.2)

# final model
noneCn = Cn.frame$Cn[which(Cn.frame$DroppedVar == "<none>")]
which.final = which(apply(SSPmat.d, 2, mean) < noneCn & pVal < 0.05)
fixed.final = paste(varnames[which.final], collapse="+")
form.final = as.formula(paste("log(PRCP+1) ~", fixed.final, random_terms))
mod.final = lmer(form.final, data=rainsmall)
summary(mod.final)
anova(mod.final)
r.squaredGLMM(mod.final)
anova(mod.final, mod.full)

# out-of-sample prediction
set.seed(07222015)

pred.mat = matrix(0,1e2,2)
system.time(for(i in 1:1e2){
test = sample(1:n, ceiling(.1*n), replace=F)
mod.full.train = update(mod.full, subset=-test)
mod.final.train = update(mod.final, subset=-test)

ytest = log(rainsmall$PRCP[test]+1)
yhat.full = predict(mod.full.train, newdata=rainsmall[test,])
pred.mat[i,1] = mean((ytest - yhat.full)^2)

yhat.final = predict(mod.final.train, newdata=rainsmall[test,])
pred.mat[i,2] = mean((ytest - yhat.final)^2)
})

exp(apply(pred.mat, 2, median))-1
apply(pred.mat, 2, sd)

# future prediction
testyrs = 2003:2012
ntest = length(testyrs)
pred.mat = matrix(0, ncol=2, nrow=ntest)

for (i in 1:ntest){
iyr = testyrs[i]
itrain = which(rainsmall$year >= iyr-25 & rainsmall$year < iyr)
itest = which(rainsmall$year == iyr)
testX = rainsmall[itest,-c(1:3)]

#mod.full.train = lm(as.formula(paste("log(PRCP+1) ~", fixed.full)), data=rainsmall, subset=itrain)
#mod.final.train = lm(as.formula(paste("log(PRCP+1) ~", fixed.final)), data=rainsmall, subset=itrain)

mod.full.train = lmer(form.full, data=rainsmall, subset=itrain)
mod.final.train = lmer(form.final, data=rainsmall, subset=itrain)

ytest = log(rainsmall$PRCP[itest]+1)
yhat.full = as.matrix(cbind(1, testX)) %*% as.numeric(fixef(mod.full.train))
pred.mat[i,1] = mean((ytest - yhat.full)^2)

yhat.final = as.matrix(cbind(1, testX[,which.final])) %*% as.numeric(fixef(mod.final.train))
pred.mat[i,2] = mean((ytest - yhat.final)^2)
}

plot(pred.mat[,1]~testyrs, type="b", ylim=c(0,ceiling(max(pred.mat[,1]))), lwd=2,
	xlab="year", ylab="MSE of BLUE", main="25 year rolling prediction of next year's median rainfall")
lines(pred.mat[,2]~testyrs, type="b", lty=2, lwd=2)
legend("topleft", c("Full model", "Reduced model"), lty=1:2, lwd=2)

plot(density(ytest), xlim=c(-2,10), ylim=c(0,.5), lwd=2,
	xlab="log(PRCP+1)", ylab="density", main="Year 2012")
lines(density(yhat.full), col='red', lwd=2)
lines(density(yhat.final), col='blue', lwd=2)
legend("topright", c("Truth", "Full model pred", "Reduced model pred"),
	col=c('black','red','blue'), lty=1, lwd=2)

