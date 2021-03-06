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
rainsmall = read.csv("../data/rainsmall.csv")

# check full model
rainsmall[-(1:3)] = scale(rainsmall[-(1:3)])
varnames = names(rainsmall)[-(1:3)]
formula = paste(varnames, collapse="+")
random_terms = "+ (1|year)"
formula = as.formula(paste("log(PRCP+1) ~", formula, random_terms))
mod.full = lmer(formula, data=rainsmall)
summary(mod.full)
anova(mod.full)
r.squaredGLMM(mod.full)

# set up residuals
y = getME(mod.full, 'y')
x = getME(mod.full, 'X')
fixed = x %*% fixef(mod.full)
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
  pVal[i] = t.test(SSPmat.d[,i], SSPmat.d[,i+1], paired=TRUE)$p.value
}

Cn.frame = data.frame(DroppedVar = c(paste("-", names(data.frame(x))[-1]), "<none>"),
                      Cn = apply(SSPmat.d, 2, mean),
                      pValue = pVal)
Cn.frame = Cn.frame[with(Cn.frame, order(Cn)),]
row.names(Cn.frame) = NULL
Cn.frame

write.table(Cn.frame, 'wild_ntothepoint1.txt')

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