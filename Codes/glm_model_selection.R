## Sample analysis
## Adult income data from UCI repository
rm(list=ls())
#setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/Rain")
setwd('C:/Study/My projects/Climate-indian monsoon')
source('misc_functions.R')

library(rrcov)
library(fda.usc)

# readin data
rain = read.csv("rain_1973_2013_test.csv")
rainsmall = with(rain, data.frame(cbind(PRCP,
	LATITUDE, LONGITUDE, log(ELEVATION), TMAX, TMIN,
	del_TT_Deg_Celsius, DMI, Nino34,
	Air_200_Temp, Air_600_Temp,
	u_wind_200, u_wind_600, u_wind_850,
	v_wind_200, v_wind_600, v_wind_850)))

rainsmall = aggregate(cbind(PRCP,
	LATITUDE, LONGITUDE, ELEVATION, TMAX, TMIN,
	del_TT_Deg_Celsius, DMI, Nino34) ~ year+STATION_NAME,
	data=rain, FUN=mean)

yr.vec = 1980:2013
Cn.mat = matrix(0, ncol=9, nrow=length(yr.vec))

nboot = 1e2

for(yr in yr.vec){
  rainsmall_yr = rainsmall[which(rainsmall$year == yr),]

  X = as.matrix(rainsmall_yr[,-(1:3)])
  y = log(rainsmall_yr$PRCP+1)
  beta = as.numeric(solve(t(X) %*% X) %*% t(X) %*% y)
  #hist(y)

  # bootstrap function
  n = nrow(X)
  m = ceiling(n^(.9))
  p = ncol(X)

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
	XtX.inv = solve(t(iX)%*%iX)
	ibeta = XtX.inv %*% t(iX) %*% iyb
	beta.mat[i,] = ibeta
  }
  #SSPmat.d[,p+1] = 1/(1 + diag((beta.mat - beta) %*% XtX.inv %*% t(beta.mat - beta)))
  SSPmat.d[,p+1] = mdepth.RP(beta.mat, beta.mat)$dep
  
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
  }

  Cn.mat[yr-1979,] = apply(SSPmat.d, 2, mean)
  print(paste('done',yr.vec[yr-1979]))
}

plot(Cn.mat[,9], type='l', lwd=2, ylim=c(0,1))
for(i in 1:8){
  lines(Cn.mat[,i], lwd=2, lty=2, col=i+1)
}
legend('topright', names(rainsmall_yr)[-(1:3)], col=2:9, lty=1, lwd=2)

write.table(Cn.mat, file='Cnmat.txt')
