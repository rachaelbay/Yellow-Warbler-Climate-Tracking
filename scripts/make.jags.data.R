make.jags.data<-function(write=FALSE, lag2=FALSE) {
	require(dplyr)

	# dim[1] = year
	# dim[2] = weather variable
	# dim[3] = route
	clim.dat <- readRDS("data/PopulationTrendsClimateArray.rds") 
	clim.dat <- adply(clim.dat, c(3,1), .id = c("rte", "Year"))
	
	mn.clim <- aggregate(cbind(AverageMonthlyTemp, AverageMonthlyPrecip) ~ rte, clim.dat, mean)
	names(mn.clim)[2:3] <- c("mn.temp", "mn.precip")
	clim.dat <- merge(clim.dat, mn.clim)
	clim.dat$temp.dev <- clim.dat$AverageMonthlyTemp - clim.dat$mn.temp
	clim.dat$precip.dev <- clim.dat$AverageMonthlyPrecip - clim.dat$mn.precip
	clim.dat$year <-  as.numeric(as.factor(clim.dat$Year)) + 1 ### CHANGE FOR LAGGED EFFECT!!!!
		
	# match climate to bird data
	dat <- merge(clim.dat, data.frame(rte = bbs.dat$rte, year = bbs.dat$year, 
	                                       count = bbs.dat$count, obser = bbs.dat$obser,
	                                       firstyr = bbs.dat$firstyr, strat = bbs.dat$strat))

	if (lag2) {
		climdat.lag2<-clim.dat
		climdat.lag2$year <-  as.numeric(as.factor(clim.dat$Year)) + 2 ### CHANGE FOR LAGGED EFFECT!!!!
		climdat.lag2<-dplyr::rename(climdat.lag2, precip.dev2=precip.dev, AverageMonthlyPrecip2= AverageMonthlyPrecip)
		dat2 <- merge(dat, select(climdat.lag2, rte, year, precip.dev2, AverageMonthlyPrecip2))
	 	dat<-dat2
		dat$precip.dev.std2 <- as.numeric(scale(dat$precip.dev2))
	}
	
	if (lag2==FALSE) {dat$precip.dev.std2<-NA ; dat$AverageMonthlyPrecip2<-NA}
	
	dat$Year <- as.numeric(as.character(dat$Year)) + 1 ### CHANGE FOR LAGGED EFFECT!!!!
	dat <- subset(dat,Year>=1968) ### CHANGE FOR LAGGED EFFECT!!!!
	dat$year <- dat$year-2 #min(dat$year)+1 ### CHANGE FOR LAGGED EFFECT!!!! #'year is an index, this gets things to start at 1.
	dat$obser <- as.numeric(as.factor(dat$obser))### CHANGE FOR LAGGED EFFECT!!!!
	
	dat$temp.std <- (dat$AverageMonthlyTemp - mean(dat$AverageMonthlyTemp))/sd(dat$AverageMonthlyTemp)
	dat$precip.std <- (dat$AverageMonthlyPrecip - mean(dat$AverageMonthlyPrecip))/sd(dat$AverageMonthlyPrecip)
	
	dat$temp.dev.std <- (dat$temp.dev - mean(dat$temp.dev))/sd(dat$temp.dev)
	dat$precip.dev.std <- (dat$precip.dev - mean(dat$precip.dev))/sd(dat$precip.dev)

	dat$rte <-factor(dat$rte)

	#---------------------------------------------------------------------------------
	# bundle data, set inits, parameters to monitor, run model
	#---------------------------------------------------------------------------------
	nyear <- max(dat$year) ### CHANGE FOR LAGGED EFFECT!!!!
	nobs <- length(unique(dat$obs)) ### CHANGE FOR LAGGED EFFECT!!!!
	ncount <- length(dat$count) ### CHANGE FOR LAGGED EFFECT!!!!
	nrte <- length(unique(dat$rte))
	
	if(write) {write.csv(dat, file="data/processed/raw_df.csv")}

	# Extract strata wide mean climate effects based on raw climate data (i.e. not subset to sites with data)	
	strats<-tapply(bbs.dat$strat, bbs.dat$rte, mean)
	
	get.mean.precip.dev<-function(strata, lag=1) {
		hld<-clim.dat[clim.dat$rte %in% names(which(strats== strata)),]
		strat.yr.mns<-tapply(hld$precip.dev, hld$Year, mean)
		lagyrs<-seq(1968-lag, length.out=nyear)
		strat.yr.mns[as.character(lagyrs)]
	}
	
	strat.precip<-(sapply(1:5, get.mean.precip.dev, lag=1)- mean(dat$precip.dev))/sd(dat$precip.dev)
	strat.precip2<-(sapply(1:5, get.mean.precip.dev, lag=2)- mean(dat$precip.dev2))/sd(dat$precip.dev2)

	# 1 - model trend + clim covars
	jags.data <- list(count = dat$count, year = dat$year, obser = dat$obser,
					  rte=as.numeric(dat$rte), nrte=nrte, 
	                  nyears = nyear,strat = dat$strat, firstyr = dat$firstyr, ### CHANGE nyears for lagged effect
	                  ncounts= ncount, areaweight = bbs.dat$areaweight, ### CHANGE ncount for lagged effect
	                  nonzeroweight = bbs.dat$nonzeroweight, nobservers = nobs, ### CHANGE nobservers for lagged effect
	                  nstrata = bbs.dat$nstrata, fixedyear = bbs.dat$fixedyear, 
	                  temp = dat$temp.dev.std, mn_temp = as.numeric(scale(dat$mn.temp)),
	                  precip = dat$precip.dev.std, mn_precip = as.numeric(scale(dat$mn.precip)),
	                  precip2 = dat$precip.dev.std2, 
	                  abs.precip = as.numeric(scale(dat$AverageMonthlyPrecip)),
	                  abs.precip2= as.numeric(scale(dat$AverageMonthlyPrecip2)),
	                  strat.precip= strat.precip, strat.precip2= strat.precip2)

	if (lag2==FALSE) {jags.data$precip2<-NULL}

	jags.data
}            