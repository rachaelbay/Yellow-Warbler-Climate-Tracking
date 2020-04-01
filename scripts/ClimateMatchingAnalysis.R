rm(list=ls())

###########################
# librariers and functions#
###########################
library(beeswarm)
library(raster)
library(lme4)
library(scales)
library(RColorBrewer)
library(geosphere)
source('~/Google Drive/YWAR project/Code/myfunctions.R', chdir = TRUE)

mround<-function(x,base){
	base*floor(x/base)
}

###############
# Read in Data#
###############
ProbSurface=readRDS("~/Google Drive/YWAR project/Data/WinterProbSurfaces.rds") # Probability surfaces across the breeding range for each sample
WinterData=readRDS("~/Google Drive/YWAR project/Data/WinteringBirdData.rds") # Data on where wintering birds were captured 
cru=read.csv("~/Google Drive/YWAR project/Data/CRUclim_adjusted.csv",header=TRUE) # Read in CRU climate data

################
# Redfine sites# # Two sites (Breeding site: -160.4055 66.43942; and Winter site: -85.67 10.78) are located over the ocean. Need to move slightly so can download bioclim variables
################
ProbSurface[which(paste(round(ProbSurface$x,4), round(ProbSurface$y,5))=="-160.4055 66.43942"),]$y=66.35 # Move breeding site slightly down so it is on land
WinterData [which(paste(WinterData$longitude, WinterData$latitude)=="-85.67 10.78"),]$longitude= -85.665 # Move winter site slightly west so it is on land 

#######################
# Create Site Database#
#######################
BreedSiteIDs=paste("B",1:length(ProbSurface[,1]),sep="_") # Names of breeding sites
AllSites=cbind(BreedSiteIDs,ProbSurface$x, ProbSurface$y)

UniqueWinSites=unique(paste(WinterData$longitude,WinterData$latitude)) # Find all Unique Winter Sites
WinSites=t(array(unlist(strsplit(UniqueWinSites," ")),dim=c(2, length(UniqueWinSites)))) # Create matrix of Winter Sites
WinSiteIDs=paste("W",1:length(WinSites[,1]),sep="_") #Names of wintering sites
WinterData$WinSiteIds=WinSiteIDs [match( paste(WinterData$longitude,WinterData$latitude), UniqueWinSites)]

AllSites=as.data.frame(rbind(AllSites,cbind(WinSiteIDs,WinSites))) # This dataset will contain the X and Y coordinates of every wintering and breeding site, along with a unique identifier
colnames(AllSites)=c("ID","X","Y")
AllSites$ID=as.character(AllSites$ID)
AllSites$X=as.numeric(as.character(AllSites$X))
AllSites$Y=as.numeric(as.character(AllSites$Y))

#####################################
# Extract Climate Data for Each Site#
#####################################
# Define breeding season as june/july (months 6,7)
# Define nonbreeding season as nov/dec/jan/feb (months 1, 2, 11, 12)
cru_analyze=cru[,c(1:8,17:20)] # Subset CRU data to just the coloumns that will be used. 

WinterData=cbind(WinterData,cru_analyze[match(paste(WinterData $longitude, WinterData $latitude),paste(cru_analyze $longitude, cru_analyze $latitude)),5:12]) # add climate to wintering sites 
BreedSiteValues= cru_analyze[which(substr(cru_analyze $ID,1,1)=="B"),] # Matrix of just the breeding sites
for (i in 5:12){
	BreedSiteValues[,i]=as.numeric(as.character(BreedSiteValues[,i]))
}
for (i in 8:15){
	WinterData[,i]=as.numeric(as.character(WinterData[,i]))
}

#############################################################
# Loop over all the sites, calculatingn all relevant metrics#
#############################################################
ProbSurface=cbind(BreedSiteIDs,ProbSurface) # Add breeding sample codes to prob surface 
ProbSurfaceEqualAbundance =rep(1,length(BreedSiteIDs))/length(BreedSiteIDs) # Make a matrix that gives equal probabilities to each breeding location 

# Arrays to store values 
IndivdualPredictedClimate=array(dim=c(0,8))
colnames(IndivdualPredictedClimate) =colnames(BreedSiteValues)[5:12]

DistanceWeightedNullClimate=array(dim=c(0,8))
colnames(DistanceWeightedNullClimate) = colnames(BreedSiteValues)[5:12]

SampleWeightedDistances=array(dim=c(0,8))
colnames(SampleWeightedDistances)= colnames(BreedSiteValues)[5:12]

DistanceWeightedDistances =array(dim=c(0,8))
colnames(DistanceWeightedDistances)= colnames(BreedSiteValues)[5:12]

for (i in 1:length(WinterData[,1])){ # loop over all the winter sites and calculate all the climate variables and distances. build resulting matrices
	ProbValues= ProbSurface[,match(WinterData[i,1], colnames(ProbSurface))] # Choose the probabiliy surface that corresponds to the wintering sample
	RescaledProbSurfaceValues = ProbValues/sum(ProbValues) # Make sure prob surface sums to 1

	# First calculate predicted individual climate values at breeding ranges
	multiplied=sweep(BreedSiteValues[,5:12],1, RescaledProbSurfaceValues,"*")  
	PredictValues=apply(multiplied,2,sum)	
	IndivdualPredictedClimate =rbind(IndivdualPredictedClimate, PredictValues)

	# Next calculated null climate values at breeding ranges, correcting for distance. 		
	SampleSiteName=WinterData[i,]$WinSiteIds # Identify the name of the site the sample correspond to 
	GeoDist=distHaversine (cbind(rep(WinterData[i,]$longitude,length(BreedSiteValues$longitude)),rep(WinterData[i,]$latitude,length(BreedSiteValues$latitude))), cbind(BreedSiteValues$longitude, BreedSiteValues$latitude))/1000 # Calculate distance of site to all Breeding sites
	RoundedDistanceProbs=cbind(mround(GeoDist,500), RescaledProbSurfaceValues)
	Ag_RoundedDistanceProbs =aggregate(RoundedDistanceProbs[,2],list(RoundedDistanceProbs[,1]),sum)
	MatchedProbs= Ag_RoundedDistanceProbs[match(mround(GeoDist,500), Ag_RoundedDistanceProbs[,1]),2]
	MatchedProbsAdjust = MatchedProbs/sum(MatchedProbs)
	Multiplied=sweep(BreedSiteValues[,5:12],1, MatchedProbsAdjust,"*") # Multiply the probability values by the distance values
	ProbWeightedDistance=apply(Multiplied,2,sum)
	DistanceWeightedNullClimate =rbind(DistanceWeightedNullClimate, ProbWeightedDistance)

	# Now calculate the climate distances, for individual predicted birds
	ClimDistances=abs(sweep(BreedSiteValues[,5:12],2,as.numeric(WinterData[i,8:15]),"-")) # Take clim distances. USE ABSOLUTE VALUE SO THIS IS DISTANCE FROM BREEDING SITE TO WINTERING SITE 
	Multiplied=sweep(ClimDistances,1, RescaledProbSurfaceValues,"*") # Multiply the probability values by the distance values
	ProbWeightedDistance=apply(Multiplied,2,sum)
	SampleWeightedDistances =rbind(SampleWeightedDistances, ProbWeightedDistance)

	# Finally, calculate the average distances for the average bird assuming equal abundances across breeding range but weight by geographic distance likelihood
	Multiplied=abs(sweep(ClimDistances,1, MatchedProbsAdjust,"*")) # Multiply the probability values by the distance values. USE ABSOLUTE VALUE SO THIS IS DISTANCE FROM BREEDING SITE TO WINTERING SITE
	ProbWeightedDistance=apply(Multiplied,2,sum)
	DistanceWeightedDistances =rbind(DistanceWeightedDistances, ProbWeightedDistance)
}

colnames(WinterData)=paste("Winter",colnames(WinterData),sep="_")
IndivdualPredictedClimate =as.data.frame(cbind(WinterData , IndivdualPredictedClimate))
DistanceWeightedNullClimate =as.data.frame(cbind(WinterData ,DistanceWeightedNullClimate))
SampleWeightedDistances =as.data.frame(cbind(WinterData , SampleWeightedDistances))
DistanceWeightedDistances =as.data.frame(cbind(WinterData ,DistanceWeightedDistances))

DistanceWeightedNullClimate $Type="Null"
IndivdualPredictedClimate $Type="Predicted"
DifferenceMatrix=DistanceWeightedDistances[16:23]-SampleWeightedDistances[,16:23]
head(DifferenceMatrix)
DifferenceMatrix =cbind(WinterData , DifferenceMatrix)

ActualPrecipValues=rbind(DistanceWeightedNullClimate, IndivdualPredictedClimate) # Make final dataset with all the information needed for modeling and graph making 
ActualPrecipValues$Winter_Monthly_Precip2= ActualPrecipValues$Winter_Monthly_Precip^2
ActualPrecipValues$Winter_Monthly_Precip_scale= scale(ActualPrecipValues$Winter_Monthly_Precip)
ActualPrecipValues$Winter_Monthly_Precip_scale2= ActualPrecipValues$Winter_Monthly_Precip_scale^2

#####################################################################################
# Make Figure: Monthly precipitation niche mapping and conduct corresponding modeling #
###################################################################################

layout(t(c(1,2,3)),widths=c(1,1,.6))
###########
# Panel 1 #
###########
par(mar=c(5,5,1,1))
M1=lmer(Monthly_Precip ~Type* Winter_Monthly_Precip_scale +Type: Winter_Monthly_Precip_scale2 + Winter_Monthly_Precip_scale2 +(1| Winter_WinSiteIds/Winter_ID),data=ActualPrecipValues, REML=FALSE)
M2=update(M1,~.-Type: Winter_Monthly_Precip_scale2) # Different 'types' have different quadratic terms
anova(M1,M2)
M2=update(M1,~.-Type: Winter_Monthly_Precip_scale) # Different 'types' have different linear terms as well. Can report the p values 
anova(M1,M2)
summary(M1)
#qqnorm(resid(M1))
#hist(resid(M1))
#plot(resid(M1)~fitted(M1)) 

Winter_Monthly_Precip_scale =seq(min( ActualPrecipValues$Winter_Monthly_Precip_scale),max(ActualPrecipValues$Winter_Monthly_Precip_scale),.01)
data=expand.grid(Winter_Monthly_Precip_scale = Winter_Monthly_Precip_scale, Type =c("Null","Predicted"), Monthly_Precip =0)
data$Winter_Monthly_Precip_scale2= Winter_Monthly_Precip_scale^2
mm=model.matrix(terms(M1), data)
data$Monthly_Precip = mm %*% fixef(M1) 
pvar1 <- diag(mm %*% vcov(M1) %*% t(mm))
data$y=data $Monthly_Precip
data$plow=data $Monthly_Precip-2*sqrt(pvar1)
data$phigh=data $Monthly_Precip +2*sqrt(pvar1)
Null=data[which(data$Type=="Null"),]
Predicted=data[which(data$Type=="Predicted"),]

BrewerCols=brewer.pal(3, "BrBG")[c(1,3)]
MyCol=ActualPrecipValues$Type
MyCol[which(MyCol=="Null")]= alpha(BrewerCols[1],.7)
MyCol[which(MyCol=="Predicted")]= alpha(BrewerCols[2],.4)

plot(ActualPrecipValues$Monthly_Precip ~ ActualPrecipValues$Winter_Monthly_Precip,col= MyCol,pch=19,ylab="Breeding Site Precipitation (mm/month)",xlab="Winter Site Precipitation (mm/month)",xlim=c(0,150),ylim=c(0,110),axes=FALSE,frame.plot=TRUE,cex.lab=1.4)

Winter_Monthly_Precip_unscale=Winter_Monthly_Precip_scale*sd(ActualPrecipValues$Winter_Monthly_Precip)+mean(ActualPrecipValues$Winter_Monthly_Precip)
polygon(x=c(Winter_Monthly_Precip_unscale, rev(Winter_Monthly_Precip_unscale)),y=c(Null $phigh, rev(Null $plow)), col= alpha(BrewerCols[1],.4),border= alpha(BrewerCols[1],.4))
polygon(c(Winter_Monthly_Precip_unscale, rev(Winter_Monthly_Precip_unscale)),c(Predicted $phigh, rev(Predicted $plow)), col= alpha(BrewerCols[2],.4),border= alpha(BrewerCols[2],.4))
lines(Winter_Monthly_Precip_unscale, Null$y,col= BrewerCols[1],lwd=4)
lines(Winter_Monthly_Precip_unscale, Predicted$y,col= BrewerCols[2],lwd=4)
legend(20,17,c("Predicted Precipitation","Null Expected Precipitation"),bty="n",cex=1.3,fill=BrewerCols[2:1])

axis(1,seq(0,150,50),cex.axis=1.3)
axis(2,seq(0,150,50),cex.axis=1.3,las=2)
text(2,107,"A.",cex=2)

###########
# Panel 2 #
###########
par(mar=c(5,5,1,1))
DifferenceMatrix$Winter_Monthly_Precip_scale=scale(DifferenceMatrix$Winter_Monthly_Precip)
DifferenceMatrix$Winter_Monthly_Precip_scale2= DifferenceMatrix$Winter_Monthly_Precip_scale^2
M1=lmer(Monthly_Precip ~ Winter_Monthly_Precip_scale + Winter_Monthly_Precip_scale2  +(1| Winter_WinSiteIds),data= DifferenceMatrix, REML=FALSE)
summary(M1)
M2=update(M1,~.-Winter_Monthly_Precip_scale2)
anova(M1,M2) 
#hist(resid(M1))
#qqnorm(resid(M1))
#plot(fitted(M1)~resid(M1))

Winter_Monthly_Precip_scale =seq(min( ActualPrecipValues$Winter_Monthly_Precip_scale),max(ActualPrecipValues$Winter_Monthly_Precip_scale),.1)
data=expand.grid(Winter_Monthly_Precip_scale = Winter_Monthly_Precip_scale, Monthly_Precip =0)
data$Winter_Monthly_Precip_scale2= Winter_Monthly_Precip_scale ^2
mm=model.matrix(terms(M1), data)
data$Monthly_Precip = mm %*% fixef(M1) 
pvar1 <- diag(mm %*% vcov(M1) %*% t(mm))
data$y=data $Monthly_Precip
data$plow=data $Monthly_Precip-2*sqrt(pvar1)
data$phigh=data $Monthly_Precip +2*sqrt(pvar1)

mycol = brewer.pal(11, "BrBG")[10]
plot(DifferenceMatrix$Monthly_Precip ~DifferenceMatrix$Winter_Monthly_Precip,xlab="Winter Site Precipitation (mm/month)",ylab="Climate Matching Index (mm/month)",pch=19,  ylim=c(-50,50),xlim=c(0,150),col= alpha(mycol,.4),axes=FALSE,frame.plot=TRUE,cex.lab=1.4)
abline(h=0,lty=2)

Winter_Monthly_Precip_unscale=Winter_Monthly_Precip_scale*sd(ActualPrecipValues$Winter_Monthly_Precip)+mean(ActualPrecipValues$Winter_Monthly_Precip)
polygon(x=c(Winter_Monthly_Precip_unscale, rev(Winter_Monthly_Precip_unscale)),y=c(data$phigh,rev(data$plow)), col= alpha(mycol,.4),border= alpha(mycol,.4))

lines(Winter_Monthly_Precip_unscale, data$y,col= mycol,lwd=4)

axis(1,seq(0,150,50),cex.axis=1.3)
axis(2,seq(-50,50,25),cex.axis=1.3,las=2)
text(2,48,"B.",cex=2)

###########
# Panel 3 #
###########
par(mar=c(5,0,1,1))
M1=lmer(Monthly_Precip ~  +(1| Winter_WinSiteIds),data= DifferenceMatrix)
#hist(resid(M1))
#qqnorm(resid(M1))
#plot(fitted(M1)~resid(M1))

beeswarm(DifferenceMatrix$Monthly_Precip,ylab="",xlab="",axes=FALSE,pch=19,col= alpha(mycol,.4),frame.plot=TRUE,ylim=c(-50,50))
abline(h=0,lty=2)
points(1,fixef(M1),cex=2,pch=19,col= mycol)
lines(c(1,1),confint(M1)[3,],lwd=3,col=mycol)
text(.55,48,"C.",cex=2)

#######################################################################################
# Now analyze the other climate variables and make graph of their (negligible) effects#
#######################################################################################

layout(rbind(c(1,2,3,4,5),c(6,7,8,9,10)))
###########
# Panel 1 #
###########
par(mar=c(2,6,1,1))

M1=lmer(Monthly_Tave ~  +(1| Winter_WinSiteIds),data= DifferenceMatrix)
#hist(resid(M1))
#qqnorm(resid(M1))
#plot(fitted(M1)~resid(M1))
summary(M1)

mycol=brewer.pal(11, "PiYG")[1]
beeswarm(DifferenceMatrix$Monthly_Tave,ylab= "",xlab="",axes=FALSE,pch=19,col= alpha(mycol,.4),frame.plot=TRUE,ylim=c(-6,6),cex=.95)
abline(h=0,lty=2)
axis(2,at=seq(-6,6,2),las=2,cex.axis=1.2)

mtext("Climate Matching Index",2,4.5)
mtext(expression(paste("Ave. Monthly Temp (",degree,"C)")),2,3,cex=.8)
text(.55,5.8,"A.",cex=2)

###########
# Panel 2 #
###########
par(mar=c(2,6,1,1))
M1=lmer(Monthly_Tmax ~  +(1| Winter_WinSiteIds),data= DifferenceMatrix)
#hist(resid(M1))
#qqnorm(resid(M1))
#plot(fitted(M1)~resid(M1))
summary(M1)

mycol=brewer.pal(11, "PiYG")[1]
beeswarm(DifferenceMatrix$Monthly_Tmax,ylab= "",xlab="",axes=FALSE,pch=19,col= alpha(mycol,.4),frame.plot=TRUE,ylim=c(-6,6),cex=.95)
abline(h=0,lty=2)
axis(2,at=seq(-6,6,2),las=2,cex.axis=1.2)

mtext("Climate Matching Index",2,4.5)
mtext(expression(paste("Max. Monthly Temp (",degree,"C)")),2,3,cex=.8)
text(.55,5.8,"B.",cex=2)

###########
# Panel 3 #
###########
par(mar=c(2,6,1,1))

M1=lmer(Monthly_Tmin ~  +(1| Winter_WinSiteIds),data= DifferenceMatrix)
#hist(resid(M1))
#qqnorm(resid(M1))
#plot(fitted(M1)~resid(M1))
summary(M1)

mycol=brewer.pal(11, "PiYG")[1]
beeswarm(DifferenceMatrix$Monthly_Tmin,ylab= "",xlab="",axes=FALSE,pch=19,col= alpha(mycol,.4),frame.plot=TRUE,ylim=c(-4.2,4.2),cex=.95)
abline(h=0,lty=2)
axis(2,at=seq(-6,6,2),las=2,cex.axis=1.2)
mtext("Climate Matching Index",2,4.5)
mtext(expression(paste("Min. Monthly Temp (",degree,"C)")),2,3,cex=.8)
text(.55,4,"C.",cex=2)

###########
# Panel 4 #
###########
par(mar=c(2,6,1,1))
M1=lmer(DifferenceMatrix$Annual_Tave ~  +(1| Winter_WinSiteIds),data= DifferenceMatrix)
#hist(resid(M1)) #Residuals completely ok
#qqnorm(resid(M1))
#shapiro.test(resid(M1))
#plot(fitted(M1)~resid(M1))
summary(M1)

mycol=brewer.pal(11, "PiYG")[1]
beeswarm(DifferenceMatrix$Annual_Tave,ylab= "",xlab="",axes=FALSE,pch=19,col= alpha(mycol,.4),frame.plot=TRUE,ylim=c(-10,10),cex=.95)
abline(h=0,lty=2)
axis(2,at=seq(-10,10,2),las=2,cex.axis=1.2)
mtext("Climate Matching Index",2,4.5)
mtext(expression(paste("Ave. Annual Temp (",degree,"C)")),2,3,cex=.8)
text(.55,9.7,"D.",cex=2)

###########
# Panel 5 #
###########
par(mar=c(2,6,1,1))

M1=lmer(DifferenceMatrix$AnnualPrecip ~  +(1| Winter_WinSiteIds),data= DifferenceMatrix)
#hist(resid(M1)) #Pretty bad residuals but not much to be done about it :(
#qqnorm(resid(M1))
#shapiro.test(resid(M1))
#plot(fitted(M1)~resid(M1))
summary(M1)

mycol=brewer.pal(11, "BrBG")[10]
beeswarm(DifferenceMatrix$AnnualPrecip,ylab= "",xlab="",axes=FALSE,pch=19,col= alpha(mycol,.4),frame.plot=TRUE,ylim=c(-800,800),cex=.85)
abline(h=0,lty=2)
axis(2,at=seq(-800,800,200),las=2,cex.axis=1.2)
mtext("Climate Matching Index",2,4.5)
mtext("Total Annual Precipitation (mm)",2,3,cex=.8)
text(.55,775,"E.",cex=2)




############################################
# Now Repeat everything with Worldclim data#
############################################
rm(list=ls())
###############
# Read stuff in again#
###############
ProbSurface=readRDS("~/Google Drive/YWAR project/Data/WinterProbSurfaces.rds") # Probability surfaces across the breeding range for each sample
WinterData=readRDS("~/Google Drive/YWAR project/Data/WinteringBirdData.rds") # Data on where wintering birds were captured 
mround<-function(x,base){
	base*floor(x/base)
}

################
# Redfine sites# # Two sites (Breeding site: -160.4055 66.43942; and Winter site: -85.67 10.78) are located over the ocean. Need to move slightly so can download bioclim variables
################
ProbSurface[which(paste(round(ProbSurface$x,4), round(ProbSurface$y,5))=="-160.4055 66.43942"),]$y=66.35 # Move breeding site slightly down so it is on land
WinterData [which(paste(WinterData$longitude, WinterData$latitude)=="-85.67 10.78"),]$longitude= -85.665 # Move winter site slightly west so it is on land 

#######################
# Create Site Database#
#######################
BreedSiteIDs=paste("B",1:length(ProbSurface[,1]),sep="_") # Names of breeding sites
AllSites=cbind(BreedSiteIDs,ProbSurface$x, ProbSurface$y)

UniqueWinSites=unique(paste(WinterData$longitude,WinterData$latitude)) # Find all Unique Winter Sites
WinSites=t(array(unlist(strsplit(UniqueWinSites," ")),dim=c(2, length(UniqueWinSites)))) # Create matrix of Winter Sites
WinSiteIDs=paste("W",1:length(WinSites[,1]),sep="_") #Names of wintering sites
WinterData$WinSiteIds=WinSiteIDs [match( paste(WinterData$longitude,WinterData$latitude), UniqueWinSites)]

AllSites=as.data.frame(rbind(AllSites,cbind(WinSiteIDs,WinSites))) # This dataset will contain the X and Y coordinates of every wintering and breeding site, along with a unique identifier
colnames(AllSites)=c("ID","X","Y")
AllSites$ID=as.character(AllSites$ID)
AllSites$X=as.numeric(as.character(AllSites$X))
AllSites$Y=as.numeric(as.character(AllSites$Y))

#####################################
# Extract Climate Data for Each Site#
#####################################
# Define breeding season as june/july (months 6,7)
# Define nonbreeding season as nov/dec/jan/feb (months 1, 2, 11, 12)

points=SpatialPointsDataFrame(coords= cbind(AllSites$X,AllSites$Y),data=AllSites,proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
climDATA=getData("worldclim",var="bio",res="2.5")
climVALUES=extract(climDATA,points)
AllSites=cbind(AllSites,climVALUES)
Breed_Index=which(substr(AllSites[,1],1,1)=="B")
Winter_Index=which(substr(AllSites[,1],1,1)=="W")

#####################
# Build monthly data#
#####################
climGenerator=function(clim_1,clim_2,clim_6,clim_7,clim_11,clim_12){
	clim1_VALUES=extract(clim_1,points)
	clim2_VALUES=extract(clim_2,points)
	clim6_VALUES=extract(clim_6,points)
	clim7_VALUES=extract(clim_7,points)
	clim11_VALUES=extract(clim_11,points)
	clim12_VALUES=extract(clim_12,points)
	summer_clim=apply(cbind(clim6_VALUES, clim7_VALUES),1,mean)
	winter_clim=apply(cbind(clim1_VALUES, clim2_VALUES, clim11_VALUES, clim12_VALUES),1,mean)
	Combine=rbind(cbind(AllSites[Breed_Index,1], summer_clim[Breed_Index]),cbind(AllSites[Winter_Index,1], winter_clim[Winter_Index]))
	as.numeric(Combine[match(Combine[,1],AllSites[,1]),2])
}

clim_1 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_prec_01.tif")
clim_2 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_prec_02.tif")
clim_6 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_prec_06.tif")
clim_7 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_prec_07.tif")
clim_11 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_prec_11.tif")
clim_12 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_prec_12.tif")
AllSites$MonthlyPre=climGenerator(clim_1,clim_2,clim_6,clim_7,clim_11,clim_12)

clim_1 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tavg_01.tif")
clim_2 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tavg_02.tif")
clim_6 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tavg_06.tif")
clim_7 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tavg_07.tif")
clim_11 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tavg_11.tif")
clim_12 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tavg_12.tif")
AllSites$MonthlyTAvg=climGenerator(clim_1,clim_2,clim_6,clim_7,clim_11,clim_12)

clim_1 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tmax_01.tif")
clim_2 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tmax_02.tif")
clim_6 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tmax_06.tif")
clim_7 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tmax_07.tif")
clim_11 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tmax_11.tif")
clim_12 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tmax_12.tif")
AllSites$MonthlyTMax=climGenerator(clim_1,clim_2,clim_6,clim_7,clim_11,clim_12)

clim_1 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tmin_01.tif")
clim_2 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tmin_02.tif")
clim_6 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tmin_06.tif")
clim_7 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tmin_07.tif")
clim_11 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tmin_11.tif")
clim_12 <- raster("~/Google Drive/YWAR project/Data/WorldClim/wc2.0_2.5m_tmin_12.tif")
AllSites$MonthlyTMin=climGenerator(clim_1,clim_2,clim_6,clim_7,clim_11,clim_12)


#############################################################
# Loop over all the sites, calculating all relevant metrics#
#############################################################

ProbSurface=cbind(BreedSiteIDs,ProbSurface) # Add breeding sample codes to prob surface 
ProbSurfaceEqualAbundance =rep(1,length(BreedSiteIDs))/length(BreedSiteIDs) # Make a matrix that gives equal probabilities to each breeding location 
BreedSiteValues=AllSites[which(substr(AllSites$ID,1,1)=="B"),] # Matrix of just the breeding sites


# Arrays to store values 
IndivdualPredictedClimate=array(dim=c(0,23))
colnames(IndivdualPredictedClimate) =paste("breed",c(paste("bio",1:19,sep="_"),"MonthlyPre","MonthlyTAvg","MonthlyTMax","MonthlyTMin"),sep="_")

DistanceWeightedNullClimate=array(dim=c(0,23))
colnames(DistanceWeightedNullClimate) =paste("breed",c(paste("bio",1:19,sep="_"),"MonthlyPre","MonthlyTAvg","MonthlyTMax","MonthlyTMin"),sep="_")

SampleWeightedDistances=array(dim=c(0,23))
colnames(SampleWeightedDistances)=paste("breeddist",colnames(DistanceWeightedNullClimate),sep="_")

DistanceWeightedDistances =array(dim=c(0,23))
colnames(DistanceWeightedDistances)=paste("breeddist",colnames(DistanceWeightedNullClimate),sep="_")

for (i in 1:length(WinterData[,1])){ # loop over all the winter sites and calculate all the climate variables and distances. build resulting matrices
	ProbValues= ProbSurface[,match(WinterData[i,1], colnames(ProbSurface))]
	RescaledProbSurfaceValues = ProbValues/sum(ProbValues) # Make sure prob surface sums to 1
	
	# First calculate predicted individual climate values at breeding ranges
	multiplied=sweep(BreedSiteValues[,4:26],1, RescaledProbSurfaceValues,"*")  
	PredictValues=apply(multiplied,2,sum)	
	IndivdualPredictedClimate =rbind(IndivdualPredictedClimate, PredictValues)

	# Next calculated null climate values at breeding ranges, assuming equal abundances and correct for distances. 		
	SampleSiteName=WinterData[i,]$WinSiteIds # Identify the name of the site the sample correspond to 
	GeoDist=distHaversine (cbind(rep(WinterData[i,]$longitude,length(BreedSiteValues$X)),rep(WinterData[i,]$latitude,length(BreedSiteValues$Y))), cbind(BreedSiteValues$X, BreedSiteValues$Y))/1000 # Calculate distance of site to all Breeding sites
	RoundedDistanceProbs=cbind(mround(GeoDist,500), RescaledProbSurfaceValues)
	Ag_RoundedDistanceProbs =aggregate(RoundedDistanceProbs[,2],list(RoundedDistanceProbs[,1]),sum)
	MatchedProbs= Ag_RoundedDistanceProbs[match(mround(GeoDist,500), Ag_RoundedDistanceProbs[,1]),2]
	MatchedProbsAdjust = MatchedProbs/sum(MatchedProbs)
	Multiplied=sweep(BreedSiteValues[,4:26],1, MatchedProbsAdjust,"*") # Multiply the probability values by the distance values
	ProbWeightedDistance=apply(Multiplied,2,sum)
	DistanceWeightedNullClimate =rbind(DistanceWeightedNullClimate, ProbWeightedDistance)

	# Now calculate the climate distances, for individual predicted birds
	WinSiteValues=AllSites[which(AllSites$ID== WinterData[i,]$WinSiteIds),]	
	ClimDistances=abs(sweep(BreedSiteValues[,4:26],2,as.numeric(WinSiteValues[,4:26]),"-")) # Take clim distances.USE ABSOLUTE VALUE SO THIS IS DISTANCE FROM BREEDING SITE TO WINTERING SITE 
	Multiplied=sweep(ClimDistances,1, RescaledProbSurfaceValues,"*") # Multiply the probability values by the distance values
	ProbWeightedDistance=apply(Multiplied,2,sum)
	SampleWeightedDistances =rbind(SampleWeightedDistances, ProbWeightedDistance)

	# Finally, calculate the average distances for the average bird assuming equal abundances across breeding range but weight by geographic distance likelihood
	Multiplied=abs(sweep(ClimDistances,1, MatchedProbsAdjust,"*")) # Multiply the probability values by the distance values. USE ABSOLUTE VALUE SO THIS IS DISTANCE FROM BREEDING SITE TO WINTERING SITE
	ProbWeightedDistance=apply(Multiplied,2,sum)
	DistanceWeightedDistances =rbind(DistanceWeightedDistances, ProbWeightedDistance)
}

IndivdualPredictedClimate =as.data.frame(cbind(WinterData ,IndivdualPredictedClimate))
DistanceWeightedNullClimate =as.data.frame(cbind(WinterData ,DistanceWeightedNullClimate))
SampleWeightedDistances =as.data.frame(cbind(WinterData , SampleWeightedDistances))
DistanceWeightedDistances =as.data.frame(cbind(WinterData ,DistanceWeightedDistances))

WinterSiteValues=AllSites[which(substr(AllSites$ID,1,1)=="W"),] # Matrix of just the winter sites
colnames(WinterSiteValues)[4:26]=paste("WintSiteValue",colnames(WinterSiteValues)[4:26],sep="_")
colnames(WinterSiteValues)[1]=c("WintSiteID")
IndivdualPredictedClimate =cbind(IndivdualPredictedClimate ,WinterSiteValues[match( IndivdualPredictedClimate$WinSiteIds, WinterSiteValues$WintSiteID),])
DistanceWeightedNullClimate =cbind(DistanceWeightedNullClimate ,WinterSiteValues[match( DistanceWeightedNullClimate$WinSiteIds, WinterSiteValues$WintSiteID),])

DistanceWeightedNullClimate $Type="Null"
IndivdualPredictedClimate $Type="Predicted"

DifferenceMatrix=DistanceWeightedDistances[8:30]-SampleWeightedDistances[,8:30]
DifferenceMatrix$SiteID=IndivdualPredictedClimate$WinSiteIds
DifferenceMatrix =cbind(DifferenceMatrix ,WinterSiteValues[match( DifferenceMatrix $SiteID, WinterSiteValues$WintSiteID),])

######################################
# NOW COME BACK TO ADDING TO THE PLOT#
######################################

ActualPrecipValues=rbind(DistanceWeightedNullClimate, IndivdualPredictedClimate)

###########
# Panel 6 #
###########
par(mar=c(1,6,1,1))

M1=lmer(breeddist_breed_MonthlyTAvg ~  +(1| WintSiteID),data= DifferenceMatrix)
#hist(resid(M1))  
#qqnorm(resid(M1))
#plot(fitted(M1)~resid(M1))
summary(M1)

mycol=brewer.pal(11, "PiYG")[1]
beeswarm(DifferenceMatrix$breeddist_breed_MonthlyTAvg,ylab= "",xlab="",axes=FALSE,pch=19,col= alpha(mycol,.4),frame.plot=TRUE,ylim=c(-5,5),cex=.95)
abline(h=0,lty=2)
axis(2,at=seq(-6,6,2),las=2,cex.axis=1.2)

mtext("Climate Matching Index",2,4.5)
mtext(expression(paste("Ave. Monthly Temp (",degree,"C)")),2,3,cex=.8)
text(.55,4.8,"F.",cex=2)


###########
# Panel 7 #
###########
par(mar=c(1,6,1,1))
M1=lmer(breeddist_breed_MonthlyTMax ~  +(1| WintSiteID),data= DifferenceMatrix)
#hist(resid(M1))  
#qqnorm(resid(M1))
#plot(fitted(M1)~resid(M1))
summary(M1)

mycol=brewer.pal(11, "PiYG")[1]
beeswarm(DifferenceMatrix$breeddist_breed_MonthlyTMax,ylab= "",xlab="",axes=FALSE,pch=19,col= alpha(mycol,.4),frame.plot=TRUE,ylim=c(-6,6),cex=.95)
abline(h=0,lty=2)
axis(2,at=seq(-6,6,2),las=2,cex.axis=1.2)

mtext("Climate Matching Index",2,4.5)
mtext(expression(paste("Max. Monthly Temp (",degree,"C)")),2,3,cex=.8)
text(.55,5.8,"G.",cex=2)

###########
# Panel 8 #
###########
par(mar=c(1,6,1,1))

M1=lmer(breeddist_breed_MonthlyTMin ~  +(1| WintSiteID),data= DifferenceMatrix)
#hist(resid(M1))
#qqnorm(resid(M1))
#plot(fitted(M1)~resid(M1))
summary(M1)

mycol=brewer.pal(11, "PiYG")[1]
beeswarm(DifferenceMatrix$breeddist_breed_MonthlyTMin,ylab= "",xlab="",axes=FALSE,pch=19,col= alpha(mycol,.4),frame.plot=TRUE,ylim=c(-4,4),cex=.95)
abline(h=0,lty=2)
axis(2,at=seq(-6,6,2),las=2,cex.axis=1.2)
mtext("Climate Matching Index",2,4.5)
mtext(expression(paste("Min. Monthly Temp (",degree,"C)")),2,3,cex=.8)
text(.55,3.85,"H.",cex=2)

###########
# Panel 9 #
###########
par(mar=c(1,6,1,1))

M1=lmer(DifferenceMatrix$breeddist_breed_bio_1/10 ~  +(1| WintSiteID),data= DifferenceMatrix)
#hist(resid(M1)) 
#qqnorm(resid(M1))
#plot(fitted(M1)~resid(M1))
summary(M1)

mycol=brewer.pal(11, "PiYG")[1]
beeswarm(DifferenceMatrix$breeddist_breed_bio_1/10,ylab= "",xlab="",axes=FALSE,pch=19,col= alpha(mycol,.4),frame.plot=TRUE,ylim=c(-10,10),cex=.95)
abline(h=0,lty=2)
axis(2,at=seq(-10,10,2),las=2,cex.axis=1.2)
mtext("Climate Matching Index",2,4.5)
mtext(expression(paste("Ave. Annual Temp (",degree,"C)")),2,3,cex=.8)
text(.55,9.75,"I.",cex=2)

###########
# Panel 10 #
###########
par(mar=c(1,6,1,1))

M1=lmer(DifferenceMatrix$breeddist_breed_bio_12 ~  +(1| WintSiteID),data= DifferenceMatrix)
#hist(resid(M1)) 
#qqnorm(resid(M1))
#plot(fitted(M1)~resid(M1))
summary(M1)

mycol=brewer.pal(11, "BrBG")[10]
beeswarm(DifferenceMatrix$breeddist_breed_bio_12,ylab= "",xlab="",axes=FALSE,pch=19,col= alpha(mycol,.4),frame.plot=TRUE,ylim=c(-900,900),cex=.85)
abline(h=0,lty=2)
axis(2,at=seq(-800,800,200),las=2,cex.axis=1.2)
mtext("Climate Matching Index",2,4.5)
mtext("Total Annual Precipitation (mm)",2,3,cex=.8)
text(.55,875,"J.",cex=2)


################################################################
# Now create monthly precipitation graphs, using worldclim data#
################################################################

ActualPrecipValues=rbind(DistanceWeightedNullClimate, IndivdualPredictedClimate)
ActualPrecipValues$WintSiteValue_MonthlyPre2= ActualPrecipValues$WintSiteValue_MonthlyPre ^2
ActualPrecipValues$WintSiteValue_MonthlyPre_scale= scale(ActualPrecipValues$WintSiteValue_MonthlyPre)
ActualPrecipValues$WintSiteValue_MonthlyPre_scale2= ActualPrecipValues$WintSiteValue_MonthlyPre_scale ^2

layout(t(c(1,2,3)),widths=c(1,1,.6))
###########
# Panel 1 #
###########
par(mar=c(5,5,1,1))
M1=lmer(breed_MonthlyPre~Type* WintSiteValue_MonthlyPre_scale +Type:WintSiteValue_MonthlyPre_scale2+ WintSiteValue_MonthlyPre_scale2+(1| WintSiteID/ID),data=ActualPrecipValues, REML=FALSE)
M2=update(M1,~.-Type: WintSiteValue_MonthlyPre_scale2) # Different 'types' have different quadratic terms
anova(M1,M2)
M2=update(M1,~.-Type: WintSiteValue_MonthlyPre_scale) # Different 'types' have different linear terms as well. Can report the p values 
anova(M1,M2)
summary(M1)
#qqnorm(resid(M1))
#plot(resid(M1)~fitted(M1)) 

WintSiteValue_MonthlyPre_scale =seq(min( ActualPrecipValues$WintSiteValue_MonthlyPre_scale),max(ActualPrecipValues$WintSiteValue_MonthlyPre_scale),.01)
data=expand.grid(WintSiteValue_MonthlyPre_scale = WintSiteValue_MonthlyPre_scale, Type =c("Null","Predicted"), breed_MonthlyPre =0)
data$WintSiteValue_MonthlyPre_scale2= WintSiteValue_MonthlyPre_scale ^2
mm=model.matrix(terms(M1), data)
data$breed_MonthlyPre = mm %*% fixef(M1) 
pvar1 <- diag(mm %*% vcov(M1) %*% t(mm))
data$y=data $breed_MonthlyPre
data$plow=data $breed_MonthlyPre-2*sqrt(pvar1)
data$phigh=data $breed_MonthlyPre +2*sqrt(pvar1)
Null=data[which(data$Type=="Null"),]
Predicted=data[which(data$Type=="Predicted"),]

BrewerCols=brewer.pal(3, "BrBG")[c(1,3)]
MyCol=ActualPrecipValues$Type
MyCol[which(MyCol=="Null")]= alpha(BrewerCols[1],.7)
MyCol[which(MyCol=="Predicted")]= alpha(BrewerCols[2],.4)

plot(ActualPrecipValues$breed_MonthlyPre~ ActualPrecipValues$WintSiteValue_MonthlyPre,col= MyCol,pch=19,ylab="Breeding Site Precipitation (mm/month)",xlab="Winter Site Precipitation (mm/month)",xlim=c(0,150),ylim=c(0,110),axes=FALSE,frame.plot=TRUE,cex.lab=1.4)

Winter_Monthly_Precip_unscale= WintSiteValue_MonthlyPre_scale*sd(ActualPrecipValues$WintSiteValue_MonthlyPre)+mean(ActualPrecipValues$WintSiteValue_MonthlyPre)
polygon(x=c(Winter_Monthly_Precip_unscale, rev(Winter_Monthly_Precip_unscale)),y=c(Null $phigh, rev(Null $plow)), col= alpha(BrewerCols[1],.4),border= alpha(BrewerCols[1],.4))
polygon(c(Winter_Monthly_Precip_unscale, rev(Winter_Monthly_Precip_unscale)),c(Predicted $phigh, rev(Predicted $plow)), col= alpha(BrewerCols[2],.4),border= alpha(BrewerCols[2],.4))
lines(Winter_Monthly_Precip_unscale, Null$y,col= BrewerCols[1],lwd=4)
lines(Winter_Monthly_Precip_unscale, Predicted$y,col= BrewerCols[2],lwd=4)
legend(20,17,c("Predicted Precipitation","Null Expected Precipitation"),bty="n",cex=1.3,fill=BrewerCols[2:1])

axis(1,seq(0,150,50),cex.axis=1.3)
axis(2,seq(0,150,50),cex.axis=1.3,las=2)
text(2,107,"A.",cex=2)

###########
# Panel 2 #
###########

par(mar=c(5,5,1,1))

DifferenceMatrix$WintSiteValue_MonthlyPre2= DifferenceMatrix$WintSiteValue_MonthlyPre^2
DifferenceMatrix$WintSiteValue_MonthlyPre_scale=scale(DifferenceMatrix$WintSiteValue_MonthlyPre)
DifferenceMatrix$WintSiteValue_MonthlyPre_scale2= DifferenceMatrix$WintSiteValue_MonthlyPre_scale^2

M1=lmer(breeddist_breed_MonthlyPre ~ WintSiteValue_MonthlyPre  +(1| WintSiteID),data= DifferenceMatrix, REML=FALSE) # Model fails to converge if you include quadratic effect, so we have to ignore it. Especially because looking at the data there do not seem to be quadratic effects 
M2=update(M1,~.-WintSiteValue_MonthlyPre)
anova(M1,M2)
#hist(resid(M1)) 
#qqnorm(resid(M1))
#plot(fitted(M1)~resid(M1))
summary(M1)

WintSiteValue_MonthlyPre =seq(min( ActualPrecipValues$WintSiteValue_MonthlyPre),max(ActualPrecipValues$WintSiteValue_MonthlyPre),.1)
data=expand.grid(WintSiteValue_MonthlyPre = WintSiteValue_MonthlyPre, breeddist_breed_MonthlyPre =0)
mm=model.matrix(terms(M1), data)
data$breeddist_breed_MonthlyPre = mm %*% fixef(M1) 
pvar1 <- diag(mm %*% vcov(M1) %*% t(mm))
data$y=data $breeddist_breed_MonthlyPre
data$plow=data $breeddist_breed_MonthlyPre-2*sqrt(pvar1)
data$phigh=data $breeddist_breed_MonthlyPre +2*sqrt(pvar1)

mycol = brewer.pal(11, "BrBG")[10]
plot(DifferenceMatrix$breeddist_breed_MonthlyPre~DifferenceMatrix$WintSiteValue_MonthlyPre,xlab="Winter Site Precipitation (mm/month)",ylab="Climate Matching Index (mm/month)",pch=19,  ylim=c(-50,50),xlim=c(0,150),col= alpha(mycol,.4),axes=FALSE,frame.plot=TRUE,cex.lab=1.4)
abline(h=0,lty=2)
polygon(x=c(WintSiteValue_MonthlyPre, rev(WintSiteValue_MonthlyPre)),y=c(data$phigh,rev(data$plow)), col= alpha(mycol,.4),border= alpha(mycol,.4))

lines(WintSiteValue_MonthlyPre, data$y,col= mycol,lwd=4)

axis(1,seq(0,150,50),cex.axis=1.3)
axis(2,seq(-50,50,25),cex.axis=1.3,las=2)
text(2,48,"B.",cex=2)

###########
# Panel 3 #
###########
par(mar=c(5,0,1,1))
M1=lmer(breeddist_breed_MonthlyPre ~  +(1| WintSiteID),data= DifferenceMatrix)
#hist(resid(M1)) 
#qqnorm(resid(M1))
#plot(fitted(M1)~resid(M1))

beeswarm(DifferenceMatrix$breeddist_breed_MonthlyPre,ylab="",xlab="",axes=FALSE,pch=19,col= alpha(mycol,.4),frame.plot=TRUE,ylim=c(-50,50))
abline(h=0,lty=2)
points(1,fixef(M1),cex=2,pch=19,col= mycol)
lines(c(1,1),confint(M1)[3,],lwd=3,col=mycol)
text(.55,48,"C.",cex=2)

