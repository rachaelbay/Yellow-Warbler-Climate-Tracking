library(data.table)
library(tidyverse)

###Read in "twocol" format from Fluidigm
p1 <- read.delim("Panel1_TwoCol_bases.txt",row.names=1,na.strings=0,stringsAsFactors=F,colClasses="character")
p2 <- read.delim("Panel2_TwoCol_bases.txt",row.names=1,na.strings=0,stringsAsFactors=F,colClasses="character")
p1p2 <- data.frame(cbind(p1,p2))
p3 <- read.delim("Panel3_TwoCol_bases_08.04.17.txt",row.names=1,na.strings=0,stringsAsFactors=F,colClasses="character")


###Merge three SNP panels into a single dataframe
p1p2 <- apply(p1p2,2,function(x) factor(x,levels=c("A","T","C","G","0")))
p3 <- apply(p3,2,function(x) factor(x,levels=c("A","T","C","G","0")))
match <- match(colnames(p1p2),colnames(p3))
p3order <- p3[,match]
names(p3order) <- names(p1p2)
p1p2 <- p1p2[-which(rownames(p1p2)%in%rownames(p3)),] #This removes bad quality samples from p1p2 that we resequenced in p3
all <- rbind(p1p2,p3order)
row.names(all) <- c(rownames(p1p2),rownames(p3))
all <- data.frame(all)

###Filter down to quality assays. Quality scores are based on scoring of Fluidigm plots
index <- seq(1,ncol(all),by=2)
assays <- names(all)[index]
qual <- read.delim("AssayQuality.txt")
qualord <- qual[match(assays,qual$Assay),]
mean <- apply(qualord[,4:7],1,mean)
good <- sort(c(which(mean>3)*2,(which(mean>3)*2)-1))
snpfilt <- all[,good]


##metadata
meta <- read.delim("meta_07.06.17.txt",header=T)
ordermeta <- meta[match(rownames(snpfilt),meta$Field.Number),]
ordermeta$Source <- "Fluidigm"

###Get RAD data
meta <- read.delim("../../../YWAR_RAD_meta.txt",header=T)
data <- data.frame(fread("../../../YWAR_final.012",na.strings="-1",header=F),row.names=1)
inds <- read.delim("../../../YWAR_final.012.indv",header=F)
rownames(data) <- inds[,1]
snps <- read.delim("../../../YWAR_final.012.pos",header=F)
radmeta <- meta[match(rownames(data),meta$Field_Number),]
radmeta$Breedmap <- 1
radmeta$Source <- "RAD"
radmeta$Stage <- "Breeding"
badones <- c("96N2316","96N2314") #These are actually migrants!
radmeta <- radmeta[-match(badones,radmeta$Field_Number),]
data <- data[-match(badones,rownames(data)),]

###Pull assayed snps out of RAD data
info <- read.delim("../../YWAR_assay_200snps.txt")
snpnames <- names(snpfilt)[seq(1,ncol(snpfilt),by=2)]
info$scaf <- factor(info$scaf,levels=levels(snps[,1]))
assInd <- which(outer(snps[,1],info$scaf[match(snpnames,info$Name)],"==") & outer(snps[,2],info$pos[match(snpnames,info$Name)],"=="),arr.ind=T)
radgens <- data[,assInd[,1]]
names(radgens) <- snpnames

###Break radgens into two columns and convert to bases
twocolRAD <- matrix(data=NA,nrow=nrow(radgens))
for(i in 1:ncol(radgens)) {
	h1 <- radgens[,i]
	h2 <- radgens[,i]
	h1[h1==1] <- 0
	h2[h1==2] <- 1
	h1[h1==2] <- 1
	name <- snpnames[i]
	ref <- info$ref[info$Name==name]
	alt <- info$alt[info$Name==name]
	h1[h1==0] <- as.character(ref)
	h1[h1==1] <- as.character(alt)
	h2[h2==0] <- as.character(ref)
	h2[h2==1] <- as.character(alt)
	twocolRAD <- cbind(twocolRAD,h1,h2)
}
twocolRAD <- data.frame(twocolRAD[,-1])
names(twocolRAD) <- names(snpfilt)
rownames(twocolRAD) <- rownames(radgens)

###Combine feather and RAD data
allgens <- rbind(snpfilt,twocolRAD)
rownames(allgens) <- c(rownames(snpfilt),rownames(twocolRAD))
featherMeta <- ordermeta[,c("Field.Number","State","Town","Lat","Long","Stage","Breedmap","Source")]
names(featherMeta) <- c("Field_Number","State","Town","Lat","Long","Stage","Breedmap","Source")
radMeta <- radmeta[,c("Field_Number","State","Near_Town","Lat","Long","Stage","Breedmap","Source")]
names(radMeta) <- names(featherMeta)
allmeta <- rbind(featherMeta,radMeta)
dim(allgens)
dim(allmeta)

###This is for strucuture: pull out breeding samples and filter
base2structure <- function(frame,info) {
	new <- frame
	new <- apply(new,2,as.character)
	starts <- seq(1,ncol(frame),by=2)
	for (k in starts) {
		name <- names(frame)[k]
		ref <- info$ref[info$Name==name]
		alt <- info$alt[info$Name==name]
		new[,k][new[,k]==ref] <- 0
		new[,k][new[,k]==alt] <- 1
		new[,k+1][new[,k+1]==ref] <- 0
		new[,k+1][new[,k+1]==alt] <- 1
	}
	return(data.frame(new))
}
goodBreeding <- allgens[which(allmeta$Breedmap==1),]
cov <- apply(goodBreeding,1,function(x) length(which(is.na(x)))/length(is.na(x)))
filtBreed <- goodBreeding[which(cov<0.2),]
str <- base2structure(filtBreed,info)
rownames(str) <- rownames(filtBreed)
strMeta <- allmeta[match(rownames(str),allmeta$Field_Number),]
str <- data.frame(cbind(pop=as.numeric(strMeta$Town),str))
str <- apply(str,2,as.numeric)
str[is.na(str)] <- -9
rownames(str) <- strMeta$Field_Number
#write.table(str,"Breeding_05.29.17.str",quote=F,sep="\t",col.names=F,row.names=T)
#write.table(strMeta,"Breeding_meta_05.29.17.txt",col.names=T,quote=F,row.names=F,sep="\t")

###96 snps only
panelsnps <- allgens[,colnames(allgens)%in%colnames(p3)]
panCov <- apply(panelsnps,1,function(x) length(which(is.na(x)))/length(x))
panelFilt <- panelsnps[panCov<0.2,]
panelFilt <- apply(panelFilt,2,function(x) factor(x,levels=c("A","G","C","T","0")))
panelFilt[is.na(panelFilt)] <- "0"
panelMeta <- allmeta[match(rownames(panelFilt),allmeta$Field_Number),]

###PED files for OriGen - these are a decent starting place for rubias too
breedPanel <- panelFilt[which(panelMeta$Breedmap==1),]
breedMeta <- panelMeta[which(panelMeta$Breedmap==1),]
breedFake <- rep(0,nrow(breedPanel))
breedPED <- data.frame(cbind(rownames(breedPanel),rownames(breedPanel),
	breedFake,breedFake,breedFake,breedFake,breedPanel))
breedLocs <- breedMeta[,c(1,2,3,6,5,4)]
names(breedLocs) <- c("ID","State","Town","Stage","longitude","latitude")
write.table(breedPED,"Breeding_12.19.18.ped",row.names=F,col.names=F,sep="\t")
write.table(breedLocs,"Breeding_12.19.18.loc",row.names=F,col.names=T,sep="\t")

unkPanel <- panelFilt[panelMeta$Stage%in%c("Migrant","Wintering"),]
unkMeta <- panelMeta[panelMeta$Stage%in%c("Migrant","Wintering"),]
unkFake <- rep(0,nrow(unkPanel))
unkPED <- data.frame(cbind(rownames(unkPanel),rownames(unkPanel),
	unkFake,unkFake,unkFake,unkFake,unkPanel))
unkLocs <- unkMeta[,c(1,2,3,6,5,4)]
names(unkLocs) <- c("ID","State","Town","Stage","longitude","latitude")
write.table(unkPED,"Unknowns_12.19.18.ped",row.names=F,col.names=F,sep="\t")
write.table(unkLocs,"Unknowns_12.19.18.loc",row.names=F,col.names=T,sep="\t")

panelnames <- colnames(panelsnps)[seq(1,ncol(panelsnps),by=2)]
panelInfo <- info[match(panelnames,info$Name),]
panelMap <- panelInfo[,c("scaf","Name","ref","pos")] #ref is a fake!
panelMap$ref <- rep(0,nrow(panelMap))
write.table(panelMap,"Breeding.map",row.names=F,col.names=F,sep="\t")
write.table(panelMap,"Unknowns.map",row.names=F,col.names=F,sep="\t")