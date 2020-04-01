# YWAR per-population BBS analysis
# last edited 01.11.2018, J. Saracco
# 09.16.19 R. Bay

################################################################################
# BBS data
################################################################################

library(rgdal)
library(maptools)
library(sp)
library(rgeos)
library(raster)
library(reshape2)
library(foreign)


setwd("~/Documents/MigratoryBirds/RAD/YWAR/BBS/DataFiles/")

Aou <- "6520" #This is the identifying for YWAR

# raw data files downloaded here: https://www.pwrc.usgs.gov/BBS/RawData/
# bbs data files are now provided as individual files for each State/province
files <- list.files(path="States")
states <- gsub("(.*).csv","\\1",files)
bbs.df <- data.frame(countrynum=NULL, statenum=NULL, Route=NULL, RPID=NULL, 
                     Year=NULL, Aou=NULL, count10=NULL, count20=NULL, count30=NULL, count40=NULL, 
                     count50=NULL, StopTotal=NULL, SpeciesTotal=NULL)

for (i in 1:length(files)){
#  print(files[i])
#  unzip(zipfile=paste("States/", files[i], sep=""),exdir="States/") # Only do this once!
  x <- read.csv(paste("States/", states[i], ".csv",sep=""))
  x <- x[x$AOU %in% Aou,] # 
  bbs.df <- rbind(bbs.df, x)
}
# select only regular bbs data; i.e., no double observers, replication
bbs.df$Rte <- do.call(paste, list(bbs.df$CountryNum, bbs.df$StateNum, bbs.df$Route, sep="."))
bbs.df <- bbs.df[bbs.df$RPID%in%"101",] #No replication
obs.df <- read.csv("weather.csv") 
obs.df$Rte <- do.call(paste, list(obs.df$CountryNum, obs.df$StateNum, obs.df$Route, sep="."))
obs.df <- obs.df[obs.df$RPID%in%"101" & obs.df$RunType%in%"1" & obs.df$Rte %in% bbs.df$Rte,]# 
obs.df <- obs.df[,c(1:4,6,9,23)] #Which columns are we trying to get here?
obs.df$obs <- do.call(paste, list(obs.df$ObsN, obs.df$Rte, sep="."))
obs.df$fy <- ifelse(duplicated(obs.df$obs), 0, 1) #Is this a duplicate?
bbs.df <- merge(bbs.df, obs.df, all.y=T) # this for matrix format; see below
bbs.df$SpeciesTotal[is.na(bbs.df$SpeciesTotal)] <- 0
routes = read.csv("routes.csv")
routes$Rte <- do.call(paste, list(routes$CountryNum, routes$StateNum, routes$Route, sep="."))
bbs.df <- merge(bbs.df, routes)
mig <- read.csv("MigrantNonBreeder/MigrantSummary.csv")
mig <- mig[mig$AOU%in%Aou,]
mig$Rte <- do.call(paste, list(mig$CountryNum, mig$StateNum, mig$Route, sep="."))
bbs.df <- bbs.df[!(bbs.df$Rte%in%mig$Rte),] #Get rid of nonbreeders (found in migrant file)

bbs.mn <- aggregate(SpeciesTotal~Rte, bbs.df, mean)

yrdetec <- aggregate(SpeciesTotal~Rte, bbs.df, length)

names(yrdetec)[2] <- "no.yrs"
sumtab <- merge(bbs.mn, yrdetec)

zeromean.rtes <- sumtab$Rte[sumtab$SpeciesTotal==0]

bbs.df <- bbs.df[!(bbs.df$Rte%in% zeromean.rtes),] #Get rid of routes with 0 individuals

hist(bbs.df$SpeciesTotal)
bbs.locs <- bbs.df[,c(22,21)]

# Read in regions shapefile - This is for K=5
ywar.regions <- readRDS("K5_shapefile.rds")
ywar.proj<- proj4string(ywar.regions)
ywar.proj

# # calculate region areas
ywar.reg.df <- data.frame(Reg=c("ywar1","ywar2","ywar3","ywar4","ywar5"), Area=gArea(ywar.regions, byid=T))
ywar.regions <- SpatialPolygonsDataFrame(ywar.regions, ywar.reg.df)

# Read route locs to spatial points
bbs.locs.sp <- SpatialPoints(bbs.locs, proj4string=CRS(
  "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
bbs.locs.sp <- spTransform(bbs.locs.sp, CRS(ywar.proj))

w <- seq(10000, 50000, 10000)

bbs.strat <- over(bbs.locs.sp, ywar.regions)  

for (i in w){
  bbs.locs.sp.b <- gBuffer(bbs.locs.sp, byid=TRUE, width=i)
  # overlay point and region data. 
  # Some coastal locs end up as NA because just outside polygon boundary. So
  # have to enlarge points to polygons by buffering to get strata overlay. Here I 
  # do this by 10 km increments adding classified points at each radius.
  bbs.strat.b <- over(bbs.locs.sp.b, ywar.regions)
  bbs.strat$Reg[is.na(bbs.strat$Reg)] <- bbs.strat.b$Reg[is.na(bbs.strat$Reg)]
  bbs.strat$Area[is.na(bbs.strat$Area)] <- bbs.strat.b$Area[is.na(bbs.strat$Area)]
}

bbs.df <- cbind(bbs.df, bbs.strat)

bbs.df <- bbs.df[!is.na(bbs.df$Reg),]
rm(bbs.strat, bbs.strat.b)

bbs.df$Reg.BCR <- interaction(bbs.df$Reg, bbs.df$BCR)


bbs.df <- bbs.df[bbs.df$Year>1967,]

# summarize no. routes/region
rte.sum <- unique(subset(bbs.df, select=c(Reg, Rte)))
rte.sum <- aggregate(Rte~Reg, rte.sum, length)

rte.sum

# ywar1 = arctic
# ywar2 = central
# ywar3 = west US
# ywar4 = east
# ywar5 = west canada

# data inputs needed for JAGS

count <- bbs.df$SpeciesTotal
ncounts <- length(count)
obser <- bbs.df$obs
obser <- as.numeric(factor(bbs.df$obs))
nobservers <- length(unique(obser))
firstyr <- bbs.df$fy
year <- as.numeric(factor(bbs.df$Year))
nyears <- length(unique(year))
strat <- as.numeric(factor(bbs.df$Reg))
nstrata <- length(unique(strat))
aw <- unique(subset(bbs.df, select = c(Reg, Area)))
aw <- aw[order(aw$Reg),]
areaweight <- aw$Area
totareaweight <- sum(aw$Area)

# calculate z weights 
rte.all.locs <- routes[,7:6]
rte.all.locs.sp <- SpatialPoints(rte.all.locs, 
                                 proj4string=
                                   CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
rte.all.locs.sp <- spTransform(rte.all.locs.sp, CRS(ywar.proj))


rte.strat <- over(rte.all.locs.sp, ywar.regions) 
for (i in w){
  rte.all.locs.sp.b <- gBuffer(rte.all.locs.sp, byid=TRUE, width=i)
  # overlay point and region data. 
  rte.strat.b <- over(rte.all.locs.sp.b, ywar.regions)
  rte.strat$Reg[is.na(rte.strat$Reg)] <- rte.strat.b$Reg[is.na(rte.strat$Reg)]
}

rte.all <- cbind(routes, rte.strat)
rte.all <- rte.all[!is.na(rte.all$Reg),]

rte.sum <- aggregate(Rte~Reg, rte.all, length) 
names(rte.sum)[2] <- "tot.rtes"
spec.rte.sum <- aggregate(Rte~Reg, unique(subset(bbs.df, select=c(Rte, Reg))), length)
names(spec.rte.sum)[2] <- "detec.rtes"
wts <- merge(spec.rte.sum, rte.sum)
wts$nonzeroweight <- wts$detec.rtes/wts$tot.rtes
nonzeroweight <- wts$nonzeroweight

jags.data <- list(count=count, year = year, obser=obser,  nyears=nyears, areaweight=areaweight, nonzeroweight=nonzeroweight,  firstyr=firstyr, ncounts=ncounts, strat=strat,  
                  nobservers=nobservers, nstrata=nstrata, fixedyear=25) 

saveRDS(jags.data,"jags_data_wRte_09.20.19.rds")