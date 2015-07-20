## Sample analysis
rm(list=ls())
#setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/Rain")
setwd('C:/Study/My projects/Climate-indian monsoon/Codes')

# read in data
rain = read.csv("../data/rain_1973_2013_test.csv")

rainsmall1 = aggregate(PRCP ~ year+STATION_NAME, data=rain, FUN=median)
rainsmall2 = aggregate(cbind(LATITUDE, LONGITUDE, ELEVATION, TMAX, TMIN,
                             del_TT_Deg_Celsius, DMI, Nino34,
                             u_wind_200, u_wind_600, u_wind_850,
                             v_wind_200, v_wind_600, v_wind_850) ~ year+STATION_NAME,
                       data=rain, FUN=median)
rainsmall = cbind(rainsmall1, rainsmall2)
rainsmall[,4:5] = list(NULL)
rm(rainsmall1, rainsmall2)

# Madden-Julien Oscillation data
MJO = read.table("../data/madden_julien_1978_2015.txt", header=T)
MJO$year = as.numeric(substr(paste(MJO$PENTAD),1,4))
MJO = MJO[-which(MJO$year > 2013),]
for(i in 1:ncol(MJO)){
  MJO[,i] = as.numeric(paste(MJO[,i]))
}
MJOsmall = aggregate(cbind(X20E, X70E, X80E, X100E,
                           X120E, X140E, X160E,
                           X120W, X40W, X10W)~year, data=MJO, FUN=median)
rainsmall = merge(rainsmall, MJOsmall)

# teleconnections data
tele = read.table("../data/teleconnections_1950_2015.txt", header=T)
telesmall = aggregate(cbind(NAO, EA, WP, EPNP, PNA,
                            EAWR, SCA, TNH, POL)~yyyy,
                      data=tele, FUN=median)
names(telesmall)[1] = 'year'
rainsmall = merge(rainsmall, telesmall)

# solar flux data
solar = read.table('../data/solar_flux_1948_2015.txt', header=F)
solarsmall = data.frame(year = solar[,1], SolarFlux = apply(solar[,-1], 1, median))
rainsmall = merge(rainsmall, solarsmall)

# Temperature anomaly data
temp = read.table('../data/temp_anomaly_index_1948_2012.txt', header=F)
tempsmall = data.frame(year = temp[,1], TempAnomaly = apply(temp[,-1], 1, median))
rainsmall = merge(rainsmall, tempsmall)

# save data
write.csv(rainsmall, '../Data/rainsmall.csv', row.names=FALSE)
