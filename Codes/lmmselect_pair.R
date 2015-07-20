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
m = ceiling(n^(.6))
p = ncol(rainsmall)-3
nboot = 1e3

SSPmat = matrix(0, nrow=nboot, ncol=p+1)
SSPmat.d = SSPmat
beta.mat = matrix(0, nrow=nboot, ncol=p)

# loop for full model bootstrap
# to get bootstrap distribution of parameters
set.seed(07152015)
for(i in 1:nboot){
	iind = sample(1:n, m, replace=T)
	imod = update(mod.full, subset=iind)	
	beta.mat[i,] = as.numeric(coef(imod)$year[1,-1])
}
SSPmat.d[,p+1] = mdepth.RP(beta.mat, beta.mat)$dep
  
## loop to get drop 1 bootstrap estimates
for(i in 1:nboot){
	iind = sample(1:n, m, replace=T)
    
	# loop for all variables
	for(j in 1:p){
		# build model
		jformula = paste(varnames[-j], collapse="+")
		jformula = as.formula(paste("log(PRCP+1) ~", jformula, "+ (1|year)"))
		ijmod = lmer(jformula, data=rainsmall, subset=iind)

		# extract coef and extend to lengthp by appending 0's
		ibeta.j = as.numeric(coef(ijmod)$year[1,-1])
    		ibeta.j.full = rep(0,p)
		ibeta.j.full[-j] = ibeta.j

    		# get depth-based criterion
		#SSPmat.d[i,j] = 1/(1 + t(ibeta.j.full - beta) %*% XtX.inv %*% (ibeta.j.full - beta))
		SSPmat.d[i,j] = mdepth.RP(ibeta.j.full, beta.mat)$dep
  	}
}

t = data.frame(apply(SSPmat.d, 2, mean))
rownames(t) = c(varnames,"full")
t

write.table(t, 'allmedian.txt')

# parameter distributions
pairs(beta.mat[,1:8], pch=19, cex=.2)

# final model
formula = paste(varnames[c(3,4,7,16:19,23:24)], collapse="+")
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
