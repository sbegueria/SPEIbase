library(ncdf4)

#ff <- list.files('..','SPEIpen',full.names=TRUE)
#system(paste('gzip -d',ff))

# open the netCDF files to read
ff <- list.files('./outputNcdf','*.nc',full.names=TRUE)
d <- list()
for (i in 1:length(ff)) d[[i]] <- nc_open(ff[i])

# a matrix of coordinates for which we want to extract the data
co <- matrix(NA,ncol=2,nrow=720*360)
colnames(co) <- c('lon','lat')
co[,1] <- rep(d[[1]]$dim$lon$vals,360)
co[,2] <- sort(rep(d[[1]]$dim$lat$vals,720),decreasing=TRUE)
co <- co[83260:259200,]

toText <- function(co, d) {
	require(ncdf4)
  lons <- d[[1]]$dim$lon$vals
	lats <- d[[1]]$dim$lat$vals
	tims <- d[[1]]$dim$time$vals
	lo <- co[1]
	la <- co[2]
  if (file.exists(paste('./outputTxt/spei_',lo,'_',la,'.csv.gz',sep='')))
    return()
	x <- which(lons==lo)
	y <- which(lats==la)
	z <- ncvar_get(d[[1]],'spei',c(x,y,1),c(1,1,1))
	t <- length(tims)
	if (is.na(z)) {
		rm(lons,lats,lo,la,x,y,z)
		return()
	}
	a <- matrix(NA,ncol=48,nrow=t)
	scales <- formatC(1:48, width=2, format='d', flag='0')
	colnames(a) <- paste('SPEI',scales,sep='')
	for (s in 1:48) {
		a[,s] <- ncvar_get(d[[s]],'spei',c(x,y,1),c(1,1,-1))
	}
	b <- 'Description: Standardized Precipitation Evapotranspiration Index (SPEI) at timescales between 1 and 48 months. Each column is a time series of monthly SPEI values starting on January 1901. Units: standard (Z) scores, i.e. normally distributed values with mean 0 and unit standard deviation. Missing value: NA. Not a number: nan. Version: 2.5 - July 2017. Comment: Using Penman equation for ET. Creators: Santiago Beguería - santiago.begueria_add_the_at_symbol_csic.es and Sergio Vicente-Serrano - svicen_add_the_at_symbol_ipe.csic.es.'
	c <- paste('./outputTxt/spei_',lo,'_',la,'.csv',sep='')
	write.table(b,c,row.names=FALSE,col.names=FALSE)
	write.table(round(a,5),c,append=TRUE,row.names=FALSE,col.names=TRUE,sep=',',quote=FALSE)
	system(paste('gzip',c))
	#return(print(paste(lo,la,sep=', ')))
	rm(lons,lats,lo,la,x,y,z,a,b,c)
	gc()
	return()
}
# toText(co[20000,],d)

apply(co, MARGIN=1, FUN=toText, d)

# Parallel: doesn't work due to multiple access to same files
#library(snowfall)
#sfInit(parallel=TRUE,cpus=10)
#sfApply(co,margin=1,fun=toText,d)
#sfStop()

#system(paste('gzip',ff))

