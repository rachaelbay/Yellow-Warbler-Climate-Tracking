######JUST START HERE
setwd("~/Google Drive/YWAR project/Data/")
library(Pstat)
library(nlme)
library(tidyverse)

breed <- read.csv("Morphology.csv")
breed <- breed[!is.na(breed$K5cur) & !is.na(breed$TARSUS),]
breed <- breed[-57,] # remove bird with unknown precip
morph <- breed[,13:21] # just the morphology data


##Get wintering bird climates
clim <- read.csv("CRUclim_adjusted.csv")
names(clim)
orderclim <- clim[match(breed$ID1,clim$ID),]
clim <- as.numeric(as.character(orderclim$Monthly_Precip))

##Scale by tarsus
fixmorph <- Res(data=cbind(breed$K5cur,morph),r="TARSUS")

##Color scale
my_cols <- c('#984ea3','#ff7f00','#e41a1c','#FF69B4','#377eb8','#ffff33')
ordercol <- my_cols[c(2,3,1,5,4)]

###gls with xy correlation structure to account for spatial autocorrelation.
par(mfrow=c(1,2),cex.lab=1.4)
pvals <- c()
names <- c(NA,"Scaled Bill Length",
           "Scaled Bill Width",
           "Scaled Bill Depth",NA,
           "Scaled Sixth Primary",
           "Scaled Ninth Primary",
           "Scaled Tail Length")
for (k in c(2,3,4,6,7,8)) {
frame <- data.frame(cbind(trait=fixmorph[,k],clim,x=breed$LONG,y=breed$LAT))
frame <- na.omit(frame)
M1 <- gls(trait~clim,correlation=corExp(form=~x+jitter(y,amount=0.01),nugget=TRUE),data=frame)
p <- coefficients(summary(M1))[2,4]
pvals[k] <- p
x <- seq(0,max(clim),by=0.5)
plotframe <- data.frame(cbind(clim=x,trait=0))
mm=model.matrix(terms(M1),plotframe)
plotframe$clim <- mm %*% coef(M1)
pvar <- diag(mm %*% vcov(M1) %*% t(mm))
plow <- plotframe$clim-2*sqrt(pvar)
phigh <- plotframe$clim+2*sqrt(pvar)

###Plot
plot(fixmorph[,k]~clim,col=ordercol[as.numeric(breed$K5cur)],
     pch=19,xlab="Monthly Precipitation (mm)",ylab=names[k])
if (k%in%c(2,4)) { # These are the significant ones
lines(x,plotframe$clim)
polygon(c(rev(x),x),c(rev(phigh),plow),col=alpha("black",0.2),border=NA)
}
}
p.adj <- p.adjust(pvals,method="fdr")
p.adj




