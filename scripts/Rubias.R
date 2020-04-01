library(devtools)
library(rubias)
library(tidyverse)
library(rworldmap)

str <- read.delim("Breeding_157snps.str",header=F)
inds <- str[,1]
K5 <- read.delim("K5.txt",row.names=1,header=F) # This is the Q matrix from structure
rownames(K5) <- inds

###Read and order meta
raw <- read.delim("Breeding_meta_05.29.17.txt")
meta <- raw[match(rownames(K5),raw$Field_Number),]
meta$Pop <- as.numeric(meta$Town)

###Find mode for each population
agg <- aggregate(K5,list(meta$Pop),mean)
agg$maxPop <- apply(agg[,2:6],1,function(x) which(x==max(x))) #CHANGE THIS FOR DIFFERENT K VALUES
meta$maxPop <- agg$maxPop[match(meta$Pop,agg$Group.1)]

###Double check assignments by plotting
popcoords <- aggregate(meta[,c("Long","Lat")],by=list(meta$Pop),mean)
newmap <- getMap(resolution="low")
plot(newmap,lwd=1.5,border="gray30",xlim=c(-150,-50),ylim=c(30,70),axes=T)
text(popcoords[,2],popcoords[,3],agg$maxPop,col="red",cex=0.8)
#text(popcoords[,2],popcoords[,3],agg$Group.1,col="blue",cex=0.5)

###Read in genotypes
breed <- read.delim("Breeding_96snps.ped",header=F,row.names=1)
breed <- breed[,6:ncol(breed)]
snps <- read.delim("Breeding_96snps.map",header=F)[,2]
meta <- meta[match(rownames(breed),meta$Field_Number),]

### Do some formatting
breed[breed==0] <- NA
colnames(breed) <- rep(snps,each=2)
#breedFrame <- data.frame(cbind(sample_type=rep("reference",nrow(breeding)),repunit=meta$maxPop,collection=meta$Source,indiv=meta$Field_Number,breeding))
breedFrame <- data.frame(cbind(sample_type=rep("reference",nrow(breed)),repunit=meta$maxPop,collection=meta$maxPop,source=meta$Source,indiv=meta$Field_Number,breed))
breedFrame$repunit <- as.factor(breedFrame$repunit)
breedFrame$collection <- as.character(breedFrame$collection)
breedFrame$indiv <- as.character(breedFrame$indiv)
#breedFrame$sample<- as.character(breedFrame$sample)
breedFrame <- breedFrame[!is.na(breedFrame$repunit),]

breedFrame <- lapply(breedFrame,as.character) %>%
	as_data_frame()

###Do Self Assignment
self <- self_assign(breedFrame,gen_start_col=6)

max <- self %>%
	group_by(indiv) %>%
	filter(near(log_likelihood,max(log_likelihood)))
	
###Add meta data
frame <- max %>%
	left_join(., breedFrame %>% select(indiv,source))

###
tally <- frame %>%
	filter(source=="Fluidigm") %>%
	group_by(collection,inferred_collection) %>%
	tally() %>%
	rename(from=collection) %>%
	spread(key=inferred_collection,value=n,fill=0)
	
###Filter
filt <- frame %>% 
#	filter(source=="Fluidigm") %>%
	filter(scaled_likelihood>0.9) %>%
	group_by(collection,inferred_collection) %>%
	tally() %>%
	rename(from=collection) %>%
	spread(key=inferred_collection,value=n,fill=0)
	
###That looks pretty good!
###Now lets get together the migrants and wintering birds
###Read in genotypes
unk <- read.delim("Unknowns_96snps.ped",header=F,row.names=1)
unk <- unk[,6:ncol(unk)]
unk[unk==0] <- NA
unkMeta <- read.delim("Unknowns_96snps.loc",header=T)


### Format wintering data for rubias
colnames(unk) <- rep(snps,each=2)
unkFrame <- data.frame(cbind(sample_type=rep("mixture",nrow(unk)),
	repunit=rep(NA,nrow(unk)),collection=rep("mixture",nrow(unk)),
	source=rep("Fluidigm",nrow(unk)),indiv=rownames(unk),unk))
unkFrame$repunit <- as.factor(unkFrame$repunit)
unkFrame$collection <- as.character(unkFrame$collection)
unkFrame$indiv <- as.character(unkFrame$indiv)

unkFrame <- lapply(unkFrame,as.character) %>%
	as_data_frame()

### Run mixture model to assign wintering and migrating birds
mix <- infer_mixture(reference=breedFrame,mixture=unkFrame,gen_start_col=6,method="MCMC")

###
mixMax <- mix$indiv_posteriors %>%
	group_by(indiv) %>%
	filter(near(PofZ,max(PofZ)))

mixMax %>%
	group_by(collection) %>%
	tally()


###Add meta data
frame <- cbind(unkMeta,collection=mixMax$collection,PofZ=mixMax$PofZ)
saveRDS(frame,"Unknown_assignments_K5.rds")
