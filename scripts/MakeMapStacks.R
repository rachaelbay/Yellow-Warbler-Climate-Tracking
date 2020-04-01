library(RColorBrewer)
library(rgdal)
library(tess3r) #This is Eric's version of Tess for mapping rasters
library(tidyverse)
library(raster)
library(ggmap)
library(ggspatial)


####Read in qmatrices
str <- read.delim("Breeding_157snps.str",header=F) # This is the structure input file (get the names from here)
inds <- str[,1]
K5 <- read.delim("K5.txt",row.names=1,header=F) # This is the q matrix
rownames(K5) <- inds

###Read and order meta
raw <- read.delim("Breeding_meta_05.29.17.txt")
meta <- raw[match(rownames(K5),raw$Field_Number),]
meta$Pop <- as.numeric(meta$Town)
  
popcoords <- aggregate(meta[,c("Long","Lat")],list(meta$Pop),mean)
ordercoords <- popcoords[order(popcoords$Long),]
ordercoords$Fakepop <- 1:nrow(ordercoords)
meta$Fakepop <-ordercoords$Fake[match(meta$Pop,ordercoords$Group.1)]
K5order <- K5[order(meta$Fakepop),]
ordermeta <- meta[order(meta$Fakepop),]

####YWAR shapefiles for the range
range <- readOGR("../data/shapefile/",layer="YWAR") # Shapefile from NatureServe
breeding <- subset(range,SEASONAL==2)
breeding <- crop(breeding,extent(-170,-55,30,70))
palette <- CreatePalette(brewer.pal(8,"Set2"))

###Do some interpolation - for YWAR the tps interpolation worked better
plotK <- K5order
stack <- tess3Q_map_rasters(as.matrix(plotK),as.matrix(ordermeta[,c("Long","Lat")]),method="map.max", interpol = FieldsTpsModel(),  
     main = "Ancestry coefficients",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4,
     col.palette = palette, map.polygon = breeding)
#plot(stack)
names(stack) <- paste0("grp",1:ncol(plotK))
maxes <- max(stack)
keep <- stack>=maxes
stack[keep<=0] <- NA

saveRDS(stack,"K5breeding.rds")

