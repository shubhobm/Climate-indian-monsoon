## depth-based model selection
library(plyr)
library(mgcv)
library(car)

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

X = as.matrix(bmd.spine[,-8])
y = as.numeric(bmd.spine[,8])

# bootstrap function
n = nrow(X)
m = ceiling(n^(.7))
p = ncol(X)
nboot = 1e3

SSPmat = matrix(0, nrow=nboot, ncol=p+1)
SSPmat.d = SSPmat
beta.mat = matrix(0, nrow=nboot, ncol=p)


# loop for full model bootstrap
# to get bootstrap distribution of parameters
for(i in 1:nboot){
  iind = sample(1:n, m, replace=T)
  iX = X[iind,]
  iyb = y[iind]
  #	ibeta = glm.fit(x=iX, y=iyb, family=binomial())$coef
  XtX.inv = solve(t(iX) %*% as.matrix(iX))
  ibeta = XtX.inv %*% t(iX) %*% iyb
  beta.mat[i,] = ibeta
}
#SSPmat.d[,p+1] = 1/(1 + diag((beta.mat - beta) %*% XtX.inv %*% t(beta.mat - beta)))

## loop to get drop 1 bootstrap estimates
for(i in 1:nboot){
  iind = sample(1:n, m, replace=T)
  iX = X[iind,]
  iyb = y[iind]
  XtX.inv = solve(t(iX)%*%iX)
  
  for(j in 1:p){
    iXj = iX[,-j]
    ibeta.j = solve(t(iXj) %*% iXj) %*% t(iXj) %*% iyb
    
    # get depth-based criterion
    ibeta.j.full = rep(0,p)
    ibeta.j.full[-j] = ibeta.j
    #    SSPmat.d[i,j] = 1/(1 + t(ibeta.j.full - beta) %*% XtX.inv %*% (ibeta.j.full - beta))
    SSPmat.d[i,j] = mdepth.RP(ibeta.j.full, beta.mat)$dep
  }
  
  # for full model
  ibeta = XtX.inv %*% t(iX) %*% iyb
  SSPmat.d[i,p+1] = mdepth.RP(as.numeric(ibeta), beta.mat)$dep
  
}

apply(SSPmat.d, 2, mean)
print(paste('done',yr.vec[yr-1979]))
