## Sample analysis
## Adult income data from UCI repository
rm(list=ls())
setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/Rain")
source('misc_functions.R')

library(rrcov)
library(fda.usc)

# readin data
rain = read.csv("../data/rain_1973_2013_test.csv")
rainsmall = rain[,-(12:28)]
rainsmall$DMI = as.numeric(paste(rainsmall$DMI))


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
yr.vec = 1978:2010
p = ncol(rainsmall)-8
loadmat = matrix(0, ncol=p, nrow=length(yr.vec))

load.list.D = list(loadmat,loadmat,loadmat)
for(yr in yr.vec){
	rainyr = rainsmall[which(rainsmall$year==yr),-c(1,2,5,7,11,20,21)]
	pcmod.D = PcaRank(rainyr[which(rainyr$PRCP < 644),-4], k=k)
	load.list.D[[1]][yr-1977,] = pcmod.D@loadings[,1]
	
	pcmod.D = PcaRank(rainyr[which(rainyr$PRCP>=644 & rainyr$PRCP<1244),-4], k=k)
	load.list.D[[2]][yr-1977,] = pcmod.D@loadings[,1]

	pcmod.D = PcaRank(rainyr[which(rainyr$PRCP>=1244),-4], k=k)
	load.list.D[[3]][yr-1977,] = pcmod.D@loadings[,1]
}

### Classical PCA
load.list = list(loadmat,loadmat,loadmat)
for(yr in yr.vec){
	rainyr = rainsmall[which(rainsmall$year==yr),-c(1,2,5,7,11,20,21)]
	pcmod = PcaClassic(rainyr[which(rainyr$PRCP < 644),-4], k=k)
	load.list[[1]][yr-1977,] = pcmod@loadings[,1]
	
	pcmod = PcaClassic(rainyr[which(rainyr$PRCP>=644 & rainyr$PRCP<1244),-4], k=k)
	load.list[[2]][yr-1977,] = pcmod@loadings[,1]

	pcmod = PcaClassic(rainyr[which(rainyr$PRCP>=1244),-4], k=k)
	load.list[[3]][yr-1977,] = pcmod@loadings[,1]
}

## Plot everything
defaultPar = par()
colOpts = rep(1:(p/2), 2)
lineOpts = rep(1:2, rep(p/2,2))

### Depth PCA plot
par(mfrow=c(1,3), oma=c(5,0,4,0))
plot(load.list[[1]][,1]~yr.vec, type='l', ylim=c(-1,1),
	main="Light rainfall", xlab="Year", ylab="Loadings")
for(i in 2:p){
	lines(load.list.D[[1]][,i]~yr.vec, type='l', col=colOpts[i], lty=lineOpts[i])
}

plot(load.list[[2]][,1]~yr.vec, type='l', ylim=c(-1,1),
	main="Medium rainfall", xlab="Year", ylab="Loadings")
for(i in 2:p){
	lines(load.list.D[[2]][,i]~yr.vec, type='l', col=colOpts[i], lty=lineOpts[i])
}

plot(load.list[[3]][,1]~yr.vec, type='l', ylim=c(-1,1),
	main="Heavy rainfall", xlab="Year", ylab="Loadings")
for(i in 2:p){
	lines(load.list.D[[3]][,i]~yr.vec, type='l', col=colOpts[i], lty=lineOpts[i])
}
mtext("Depth PCA", 3, 1, outer=TRUE, cex=1.5)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend('bottom', names(rainyr[,-4]), col=colOpts, lty=lineOpts, xpd=T,
       ncol=4, bty="n")
par(defaultPar)

### Classical PCA plot
par(mfrow=c(1,3), oma=c(5,0,4,0))
plot(load.list[[1]][,1]~yr.vec, type='l', ylim=c(-1,1),
	main="Light rainfall", xlab="Year", ylab="Loadings")
for(i in 2:p){
	lines(load.list[[1]][,i]~yr.vec, type='l', col=colOpts[i], lty=lineOpts[i])
}

plot(load.list[[2]][,1]~yr.vec, type='l', ylim=c(-1,1),
	main="Medium rainfall", xlab="Year", ylab="Loadings")
for(i in 2:p){
	lines(load.list[[2]][,i]~yr.vec, type='l', col=colOpts[i], lty=lineOpts[i])
}

plot(load.list[[3]][,1]~yr.vec, type='l', ylim=c(-1,1),
	main="Heavy rainfall", xlab="Year", ylab="Loadings")
for(i in 2:p){
	lines(load.list[[3]][,i]~yr.vec, type='l', col=colOpts[i], lty=lineOpts[i])
}
mtext("Classical PCA", 3, 1, outer=TRUE, cex=1.5)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend('bottom', names(rainyr[,-4]), col=colOpts, lty=lineOpts, xpd=T,
        ncol=4, bty="n")
par(defaultPar)

