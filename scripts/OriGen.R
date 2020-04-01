library(OriGen)
library(dplyr)
library(rgdal)

####Run OriGen with Breeding birds as training set wintering birds as unknowns
Real <- ConvertUnknownPEDData("Breeding_12.19.18","Breeding_12.19.18.loc","Unknowns_12.19.18")
Origen <- FitOriGenModelFindUnknowns(Real$DataArray,Real$SampleCoordinates,Real$UnknownData,MaxGridLength=70,RhoParameter=10)
PlotUnknownHeatMap(Origen,UnknownNumber=2)
unkMeta <- read.delim("Unknowns_12.19.18.loc",header=T)

#####Clip output grid using range shapefile
range <- readOGR("../data/shapefile/",layer="YWAR") # This is the shapefile retrieved from NatureServe
breeding <- subset(range,SEASONAL==2)
grid <- expand.grid(Origen$GridCoordinates[1,],Origen$GridCoordinates[2,])
grid <- grid[grid$Var2!=0,]
names(grid) <- c("x","y")
pts <- SpatialPointsDataFrame(grid[,c(1,2)],data=grid[,c(2,1)],proj4string=attributes(breeding)$proj4string)
overlap <- sp::over(pts,as(breeding,"SpatialPolygons"))
subgrid <- grid[!is.na(overlap),]
subgrid$xind <- match(subgrid$x,Origen$GridCoordinates[1,])
subgrid$yind <- match(subgrid$y,Origen$GridCoordinates[2,])


##write out probability surfaces
for (i in 1:nrow(Real$UnknownData)) {
  myGrid <- Origen$UnknownGrids[,,i]
  for (j in 1:nrow(subgrid)) {
    subgrid[j,i+4] <- myGrid[subgrid$xind[j],subgrid$yind[j]]
  }
}

winterind <- c(1,2,which(unkMeta$Stage=="Wintering")+4)
wintersurf <- subgrid[,winterind]
names(wintersurf) <- c("x","y",as.character(unkMeta$ID[unkMeta$Stage=="Wintering"]))
saveRDS(wintersurf,"WinterProbSurfaces.rds")



