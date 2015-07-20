## Sample analysis
## Adult income data from UCI repository
rm(list=ls())
setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/Rain")
#source('misc_functions.R')

## depth-based model selection
library(plyr)
library(mgcv)
library(car)
library(fda.usc)

### Read in and process data
bmd = read.table("bmd_bl.txt")
# bmd = within(bmd, {
#   treat = factor(treat)
#   gender = mapvalues(gender, from=c(1,2), to=c("M","F"))
#   vitD = mapvalues(vitD, from=c(1,2,3), to=c("low","medlow","high"))
#   steroids = mapvalues(steroids, from=c(0,1), to=c("no","yes"))
# })
bmd.hip = bmd[,-c(1,2,10,12)]
bmd.spine = bmd[,-c(1,2,10,11)]

d = read.table("diabetes.txt", header=T)
X = as.matrix(d[,-11])
y = as.matrix(d[,11])

X = as.matrix(bmd.spine[,-8])
y = as.numeric(bmd.spine[,8])

# bootstrap function
n = nrow(X)
sdn = n^.05
#m = ceiling(n^(.5))
p = ncol(X)
nboot = 1e3

SSPmat = matrix(0, nrow=nboot, ncol=p+1)
SSPmat.d = SSPmat
beta.mat = matrix(0, nrow=nboot, ncol=p)

# make full model
XtX.inv = solve(t(X) %*% X)
Const.full = XtX.inv %*% t(X)
beta = Const.full %*% y
Xbeta = X %*% beta
r = (y - Xbeta)/sqrt(1-p/n)

# loop for full model bootstrap
# to get bootstrap distribution of parameters
for(i in 1:nboot){
  iyb = Xbeta + sdn*rnorm(n)*r
  ibeta = Const.full %*% iyb
  beta.mat[i,] = ibeta
}

## initialize some quantities
Const = list()
for(j in 1:p){
	Xj = X[,-j]
	Const[[j]] = solve(t(Xj) %*% Xj) %*% t(Xj)
}

## loop to get drop 1 bootstrap estimates
for(i in 1:nboot){
  iyb = Xbeta + sdn*rnorm(n)*r
  
  for(j in 1:p){
    ibeta.j = Const[[j]] %*% iyb
    
    # get depth-based criterion
    ibeta.j.full = rep(0,p)
    ibeta.j.full[-j] = ibeta.j
    SSPmat.d[i,j] = mdepth.RP(ibeta.j.full, beta.mat)$dep
    #SSPmat.d[i,p+1] = 1/(1 + t(as.numeric(ibeta.j.full) - beta) %*% XtX.inv %*% (as.numeric(ibeta.j.full) - beta))
  }
  
  # for full model
  ibeta = Const.full %*% iyb
  SSPmat.d[i,p+1] = mdepth.RP(as.numeric(ibeta), beta.mat)$dep
  #SSPmat.d[i,p+1] = 1/(1 + t(as.numeric(ibeta) - beta) %*% XtX.inv %*% (as.numeric(ibeta) - beta))
}


# get p-values
pVal = rep(1, p+1)
for(i in 1:p){
	pVal[i] = t.test(SSPmat.d[,i], SSPmat.d[,p+1], paired=TRUE)$p.value
}

mod = glm.fit(x=X, y=y)
Cn.frame = data.frame(DroppedVar = c(paste("-", names(coef(mod))), "<none>"),
	Cn = apply(SSPmat.d, 2, mean),
	pValue = pVal)
Cn.frame = Cn.frame[with(Cn.frame, order(Cn)),]
row.names(Cn.frame) = NULL
Cn.frame
