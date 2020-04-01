library(WGCNA)
library(data.table)
library(RColorBrewer)
library(gdsfmt)
library(SNPRelate)

####Read in the following files:
## meta: metadata - make sure sample names match sample names from barcode file used in demultiplexing
## data: SNP table output from pipeline (vcftools -012 format)
## snps: SNP positions (this is also output from vcftools -012)
## inds: sample names (also vcftools -012 output)
meta <- read.delim("YWAR_RAD_meta.txt",header=T)
data <- data.frame(fread("YWAR_final.012",na.strings="-1",header=F))
inds <- read.delim("YWAR_final.012.indv",header=F)
rownames(data) <- inds[,1]
data <- data[,-1]
snps <- read.delim("YWAR_final.012.pos",header=F)
ordermeta <- meta[match(rownames(data),meta$Field_Number),]

oneper <- read.delim("longest_align.txt") # chromosome assignments for scaffolds based on alignment to zebra finch
snps[,3] <- oneper[match(snps[,1],oneper$SCAF),12]

##Order states by longitude. This helps make the graphs pretty
statelong <- aggregate(ordermeta$Long,list(ordermeta$State),mean)
statelong <- statelong[order(statelong$x),]
statelong$order <- 1:nrow(statelong)
ordermeta$order <- statelong$order[match(ordermeta$State,statelong$Group.1)]


##remove outlier individuals
bad <- c("6022","3510")
badInd <- which(ordermeta$Field_Number%in%bad)
filtmeta <- ordermeta[-badInd,]
data <- data[-badInd,]
inds <- inds[-badInd,]

###Create GDS object for SNPrelate
nonZ <- data[,which(snps[,3]!="Z")]
snpgdsCreateGeno("YWAR_auto.gds",genmat=as.matrix(nonZ),sample.id=filtmeta$Field_Number,
	snp.chromosome=snps[which(snps[,3]!="Z"),1],snp.position=snps[which(snps[,3]!="Z"),2],snpfirstdim=F)

gfile <- snpgdsOpen("YWAR_auto.gds")

###Calculate PCA
pca <- snpgdsPCA(gfile, num.thread=2,autosome.only=F)
pc.perc <- pca$varprop*100
snpgdsClose(gfile)

cols <- brewer.pal(11,"Spectral")[c(1:5,7:11)] ##You will need a different color pallette if you need >12 colors
colors <- colorRampPalette(cols)(length(unique(filtmeta$State)))

###Assayed SNPs
pan96 <- readRDS("RAD_96snps.rds") # pulled out the 96 assayed SNPs from RAD data
panmeta <- ordermeta[match(rownames(pan96),ordermeta$Field_Number),]

snpgdsCreateGeno("YWAR_rad96.gds",genmat=as.matrix(pan96),sample.id=rownames(pan157),
	snpfirstdim=F)
file96 <- snpgdsOpen("YWAR_rad96.gds")

pc96 <- snpgdsPCA(file96, num.thread=2,autosome.only=F)
pc96$varprop*100
snpgdsClose(file96)


##Plot
par(mar=c(5,4,1,1),cex.lab=1.3,mgp=c(1.5,0,0),tck=0.01)
layout(matrix(c(1,2,3,4,5,3),nrow=2,byrow=T))

#All SNPs
plot(pca$eigenvect[,2]~pca$eigenvect[,1],col=labels2colors(filtmeta$order,colorSeq=colors),
     cex=1.3,cex.lab=1.7,xlab="PC1 (26.15%)",ylab="PC2 (3.72%)",pch=19,axes=T)
plot(pca$eigenvect[,1]~filtmeta$Long,pch=19,cex=1.3,cex.lab=1.7,col=labels2colors(filtmeta$order,colorSeq=colors),xlab="PC1",ylab="Longitude")

#Legend
plot(1:10,1:10,type="n",axes=F,xlab="",ylab="")
legend("left",bty="n",cex=1.4,legend=unique(ordermeta$State[order(ordermeta$order)]),fill=labels2colors(unique(ordermeta$order[order(ordermeta$order)]),colorSeq=colors))
text(pca$eigenvect[,1],pca$eigenvect[,2],labels=ordermeta$Field_Number)

#96 SNPs
plot(-pc96$eigenvect[,2]~pc96$eigenvect[,1],col=labels2colors(panmeta$order,colorSeq=colors),
     cex=1.3,cex.lab=1.7,xlab="PC1 (34.32%)",ylab="PC2 (4.08%)",lwd=2,pch=19,axes=T)
plot(pc96$eigenvect[,1]~panmeta$Long,pch=19,cex=1.3,cex.lab=1.7,col=labels2colors(panmeta$order,colorSeq=colors),xlab="PC1",ylab="Longitude")







plot(pca$eigenvect[,1]~filtmeta$Lat,pch=19,cex=1.5,cex.lab=1.7,col=labels2colors(filtmeta$order,colorSeq=colors))
summary(lm(pca$eigenvect[,4]~filtmeta$Lat))
summary(lm(pc96$eigenvect[,2]~panmeta$Lat))

ordermeta$Group <- factor(ordermeta$Group)
snpgdsFst(gfile,sample.id=ordermeta$Field_Number,population=as.factor(ordermeta$Group),
	method="W&C84",autosome.only=F)
	
####spatial interpolation?
library(maptools)
library(sp)
library(gstat)
library(rworldmap)
randpoints <- read.delim("climate/YWARcurrent100k.txt")
Error in file(file, "rt") : cannot open the connection
In addition: Warning message:
In file(file, "rt") :
  cannot open file 'climate/YWARcurrent100k.txt': No such file or directory
PCpts <- SpatialPointsDataFrame(ordermeta[,c("Long","Lat")],data=data.frame(pca$eigenvect),proj4string=CRS("+proj=longlat"))
pts <- SpatialPointsDataFrame(randpoints[,c(1,2)],data=randpoints[,3:9],proj4string=CRS("+proj=longlat"))

int1 <- gstat::idw(X1~1,PCpts,newdata=pts)
int2 <- gstat::idw(X2~1,PCpts,newdata=pts)
int3 <- gstat::idw(X3~1,PCpts,newdata=pts)

par(mfrow=c(2,2),mar=c(0,0,3,0))
###PC1
col1 <- colorRampPalette(brewer.pal(9,"YlOrRd"))(20)[as.numeric(cut(int1$var1.pred,20))]
plot(pts,col=col1,pch=".",cex=2,main="PC1")
newmap <- getMap(resolution="low")
plot(newmap,lwd=0.5,add=T,border="grey40")

###PC2
col2 <- colorRampPalette(brewer.pal(9,"YlGnBu"))(20)[as.numeric(cut(int2$var1.pred,20))]
plot(pts,col=col2,pch=".",cex=2,main="PC2")
plot(newmap,lwd=0.5,add=T,border="grey40")

###PC1
col3 <- colorRampPalette(brewer.pal(9,"BuPu"))(20)[as.numeric(cut(int3$var1.pred,20))]
plot(pts,col=col3,pch=".",cex=2,main="PC3")
plot(newmap,lwd=0.5,add=T,border="grey40")

###Combine 3 pcs
R <- int1$var1.pred + int2$var1.pred
G <- -int2$var1.pred
B <- int3$var1.pred + int2$var1.pred - int1$var1.pred
#R <- int1$var1.pred
#G <- -int2.var1.pred
#B <- int3$var1.pred
R <- (R-min(R))/(max(R)-min(R)) * 255
G <- (G-min(G))/(max(G)-min(G)) * 255
B <- (B-min(B))/(max(B)-min(B)) * 255
plot(pts,col=rgb(R,G,B,max=255),pch=".",cex=2,main="Combined")
plot(newmap,lwd=0.5,add=T,border="grey40")

source("~/Documents/scripts/scalebar.R",int1$var1.pred)
image.scale(int1$var1.pred,col=colorRampPalette(brewer.pal(9,"YlOrRd"))(20))

