## Sample analysis
## Adult income data from UCI repository
rm(list=ls())
setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/Rain")
source('misc_functions.R')

library(rrcov)
library(fda.usc)

# readin data
rain = read.csv("rain_1973_2013_test.csv")
rainsmall = with(rain, data.frame(cbind(year, month,
	LATITUDE, LONGITUDE, STATION_NAME, log(ELEVATION+1),
	DATE, PRCP, TMAX, TMIN, PRCP_mm,
	del_TT_Deg_Celsius, DMI, Nino34,
	Air_200_Temp, Air_600_Temp,
	u_wind_200, u_wind_600, u_wind_850,
	v_wind_200, v_wind_600, v_wind_850)))

rain76X = rainsmall[which(rainsmall$year==2001),-c(1,2,5,7,8,11)]

pcmod = PcaClassic(rain76X)
pcmod.D = PcaRank(rain76X)

#pcmod@loadings
#summary(pcmod)

#pcmod.D@loadings
#summary(pcmod.D)

p = pcmod@eigenvalues / sum(pcmod@eigenvalues)
pD = pcmod.D@eigenvalues / sum(pcmod.D@eigenvalues)
plot(p, type='b', col="red", lwd=2, ylim=c(0,1))
lines(pD, type='b', col="blue", lwd=2)

# orthogonal distances
k=2
pca2 = PcaClassic(rain76X, k=k)
pca2S = PcaLocantore(rain76X, k=k)
pca2D = PcaRank(rain76X, k=k)

dist <- pca2@od^2
distS <- pca2S@od^2
distD <- pca2D@od^2

qclass  <- round(quantile(dist, probs = seq(0, 1, 0.1)[-c(1,11)]), 1)
qclassS  <- round(quantile(distS, probs = seq(0, 1, 0.1)[-c(1,11)]), 1)
qclassD  <- round(quantile(distD, probs = seq(0, 1, 0.1)[-c(1,11)]), 1)
rbind(qclass, qclassS, qclassD)

###### check how pc loadings evolve
### Depth PCA
yr.vec = 1976:2013
loadmat = matrix(0, ncol=ncol(rain76X), nrow=length(yr.vec))

load.list.D = list(loadmat,loadmat,loadmat)
for(yr in yr.vec){
	rainyr = rainsmall[which(rainsmall$year==yr),-c(1,2,5,7,11)]
	pcmod.D = PcaRank(rainyr[which(rainyr$PRCP < 644),-4], k=k)
	load.list.D[[1]][yr-1975,] = pcmod.D@loadings[,1]
	
	pcmod.D = PcaRank(rainyr[which(rainyr$PRCP>=644 & rainyr$PRCP<1244),-4], k=k)
	load.list.D[[2]][yr-1975,] = pcmod.D@loadings[,1]

	pcmod.D = PcaRank(rainyr[which(rainyr$PRCP>=1244),-4], k=k)
	load.list.D[[3]][yr-1975,] = pcmod.D@loadings[,1]
}

### Classical PCA
load.list = list(loadmat,loadmat,loadmat)
for(yr in yr.vec){
	rainyr = rainsmall[which(rainsmall$year==yr),-c(1,2,5,7,11)]
	pcmod = PcaClassic(rainyr[which(rainyr$PRCP < 644),-4], k=k)
	load.list[[1]][yr-1975,] = pcmod@loadings[,1]
	
	pcmod = PcaClassic(rainyr[which(rainyr$PRCP>=644 & rainyr$PRCP<1244),-4], k=k)
	load.list[[2]][yr-1975,] = pcmod@loadings[,1]

	pcmod = PcaClassic(rainyr[which(rainyr$PRCP>=1244),-4], k=k)
	load.list[[3]][yr-1975,] = pcmod@loadings[,1]
}

## Plot everything
defaultPar = par()
colOpts = rep(1:8, 2)
lineOpts = rep(1:2, rep(8,2))

### Depth PCA plot
par(mfrow=c(1,3), oma=c(5,0,4,0))
plot(load.list[[1]][,1]~yr.vec, type='l', ylim=c(-1,1),
	main="Light rainfall", xlab="Year", ylab="Loadings")
for(i in 2:ncol(rain76)){
	lines(load.list.D[[1]][,i]~yr.vec, type='l', colOpts[i], lty=lineOpts[i])
}

plot(load.list[[2]][,1]~yr.vec, type='l', ylim=c(-1,1),
	main="Medium rainfall", xlab="Year", ylab="Loadings")
for(i in 2:ncol(rain76)){
	lines(load.list.D[[2]][,i]~yr.vec, type='l', colOpts[i], lty=lineOpts[i])
}

plot(load.list[[3]][,1]~yr.vec, type='l', ylim=c(-1,1),
	main="Heavy rainfall", xlab="Year", ylab="Loadings")
for(i in 2:ncol(rain76)){
	lines(load.list.D[[3]][,i]~yr.vec, type='l', colOpts[i], lty=lineOpts[i])
}
mtext("Depth PCA", 3, 1, outer=TRUE, cex=1.5)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend('bottom', names(rain76), col=colOpts, lty=lineOpts, xpd=T, ncol=4, bty="n")
par(defaultPar)

### Classical PCA plot
par(mfrow=c(1,3), oma=c(5,0,4,0))
plot(load.list[[1]][,1]~yr.vec, type='l', ylim=c(-1,1),
	main="Light rainfall", xlab="Year", ylab="Loadings")
for(i in 2:ncol(rain76)){
	lines(load.list[[1]][,i]~yr.vec, type='l', colOpts[i], lty=lineOpts[i])
}

plot(load.list[[2]][,1]~yr.vec, type='l', ylim=c(-1,1),
	main="Medium rainfall", xlab="Year", ylab="Loadings")
for(i in 2:ncol(rain76)){
	lines(load.list[[2]][,i]~yr.vec, type='l', colOpts[i], lty=lineOpts[i])
}

plot(load.list[[3]][,1]~yr.vec, type='l', ylim=c(-1,1),
	main="Heavy rainfall", xlab="Year", ylab="Loadings")
for(i in 2:ncol(rain76)){
	lines(load.list[[3]][,i]~yr.vec, type='l', colOpts[i], lty=lineOpts[i])
}
mtext("Classical PCA", 3, 1, outer=TRUE, cex=1.5)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend('bottom', names(rain76), col=colOpts, lty=lineOpts, xpd=T, ncol=4, bty="n")
par(defaultPar)

## contour plot of top depth-PC
dat = rainyr[which(rainyr$PRCP<644),-4]
pcmod = PcaRank(dat, k=2)
datgrid = data.frame(as.numeric(pcmod@scores[,1]), dat$LONGITUDE ,dat$LATITUDE)
names(datgrid) = c("score","lon","lat")
datgrid1 = aggregate(score~lon+lat, data=datgrid, FUN=median)

grid = expand.grid(datgrid1$lon, datgrid1$lat); names(grid) = c("lon", "lat")
grid1 = merge(datgrid1, grid, by=c("lon","lat"), all=T)

filled.contour(x=sort(datgrid1$lon), y=unique(grid1$lat), z=matrix(grid1$score, ncol=36, byrow=T),
	plot.axes=map('world', ylim=c(8,38), xlim=c(68,98), add=T))