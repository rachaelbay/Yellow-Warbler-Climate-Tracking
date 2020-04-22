
#---------------------------------------------------------------------------------
# load packages
#---------------------------------------------------------------------------------
setwd("~/Google Drive/YWAR_BBS")
library(plyr)
library(jagsUI)
#---------------------------------------------------------------------------------
# read in bird data
#---------------------------------------------------------------------------------
rm(list=ls())

bbs.dat <- readRDS("data/jags_data_wRte_09.20.19.rds")
source('src/make.jags.data.R')
#---------------------------------------------------------------------------------
# read in weather data
#---------------------------------------------------------------------------------
jags.data<-make.jags.data(lag2=T)

inits <- function(){
  list()
}  

parameters <- c("eta", "n", "N", "sigma.noise", "sigma.obs", "sdyear", "CompIndex", "Bbar", "diffgof", "posdiff", "beta1", "beta2", "beta3", "r", "strata")

            
# MCMC settings
nsamp<-1000
ni <- 50000
na <- 10000
nb <- 20000
nc <- 4
nt <- round((ni-nb)/(nsamp/nc))

# Calling JAGS
mod <- jags(jags.data, inits, parameters, 
							"src/models/bbs_model_final.txt", 
							n.adapt = na, n.chains = nc, n.thin = nt, 
							n.iter = ni, n.burnin = nb, parallel=TRUE)

saveRDS(mod,"BBS_strata_Vfinal.rds")

mod<-readRDS("BBS_strata_Vfinal.rds")

cols<-c('mean', '2.5%', '97.5%', "Rhat", "n.eff", "overlap0", "f")

vars<-rownames(mod$summary)

varsx<-vars[grep('^beta.\\[', vars)]

# All tracked nodes
round(mod$summary[,cols], 3)

# Beta nodes
round(mod$summary[varsx,cols], 3)
