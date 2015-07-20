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

# read in data
rainsmall = read.csv("../data/rainsmall.csv")

# check full model
rainsmall[-(1:3)] = scale(rainsmall[-(1:3)])
varnames = names(rainsmall)[-(1:3)]
formula = paste(varnames, collapse="+")
formula = as.formula(paste("log(PRCP+1) ~", formula, "+ (1|year)"))
mod.full = lmer(formula, data=rainsmall)
summary(mod.full)
anova(mod.full)
r.squaredGLMM(mod.full)

# bootstrap function
n = nrow(rainsmall)
nyr = length(unique(rainsmall$year))
p = ncol(rainsmall)-3
nboot = 1e3
sdn = n^.2
sdnyr = nyr^.2

# set up residuals
y = getME(mod.full, 'y')
fixed = getME(mod.full, 'X') %*% fixef(mod.full)
eta = unlist(ranef(mod.full))
Z = t(as.matrix(getME(mod.full,'Zt')))
random = Z %*% eta
r = y - fixed - random

SSPmat = matrix(0, nrow=nboot, ncol=p+1)
SSPmat.d = SSPmat
beta.mat = matrix(0, nrow=nboot, ncol=p)

# loop for full model bootstrap
# to get bootstrap distribution of parameters
set.seed(07152015)
# create progress bar
pb <- txtProgressBar(min = 0, max = nboot, style = 3)
for(i in 1:nboot){
  iyb = fixed + (Z %*% (rnorm(nyr)*sdnyr*eta) + rnorm(n)*sd*r)/sqrt(1 - p/n)
	iformula = paste(varnames, collapse="+")
	iformula = as.formula(paste("iyb ~", iformula, "+ (1|year)"))
	imod = lmer(iformula, data=rainsmall)
	beta.mat[i,] = as.numeric(coef(imod)$year[1,-1])
  setTxtProgressBar(pb, i)
}
close(pb)
SSPmat.d[,p+1] = mdepth.RP(beta.mat, beta.mat)$dep
  
## loop to get drop 1 bootstrap estimates
# create progress bar
pb <- txtProgressBar(min = 0, max = nboot, style = 3)
for(i in 1:nboot){
  iyb = fixed + (Z %*% (rnorm(nyr)*sd*eta) + rnorm(n)*sd*r)/sqrt(1 - p/n)
    
	# loop for all variables
	for(j in 1:p){
		# build model
		jformula = paste(varnames[-j], collapse="+")
		jformula = as.formula(paste("iyb ~", jformula, "+ (1|year)"))
		ijmod = lmer(jformula, data=rainsmall)

		# extract coef and extend to lengthp by appending 0's
		ibeta.j = as.numeric(coef(ijmod)$year[1,-1])
    		ibeta.j.full = rep(0,p)
		ibeta.j.full[-j] = ibeta.j

    		# get depth-based criterion
		#SSPmat.d[i,j] = 1/(1 + t(ibeta.j.full - beta) %*% XtX.inv %*% (ibeta.j.full - beta))
		SSPmat.d[i,j] = mdepth.RP(ibeta.j.full, beta.mat)$dep
  	}
  setTxtProgressBar(pb, i)
}
close(pb)

t = data.frame(apply(SSPmat.d, 2, mean))
rownames(t) = c(varnames,"full")
t

write.table(t, 'allmedianMJOtele.txt')

# parameter distributions
pairs(beta.mat[,1:8], pch=19, cex=.2)

# final model
formula = paste(varnames[which(as.numeric(unlist(t)) < t[nrow(t),1])], collapse="+")
formula = as.formula(paste("log(PRCP+1) ~", formula, "+ (1|year)"))
mod.final = lmer(formula, data=rainsmall)
summary(mod.final)
anova(mod.final)
r.squaredGLMM(mod.final)
anova(mod.final, mod.full)

# robust lmm
######## doesn't fucking work
library(robustlmm)
rmod.full = rlmer(formula, data=rainsmall)
summary(mod.full)
anova(mod.full)
r.squaredGLMM(mod.full)
