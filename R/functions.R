#' Compute SPEI from very large arrayed datasets in netCDF format.
#' @param sca       Integer. The time scale of the SPEI.
#' @param inPre     Character vector. Path to the input precipitation netCDF file.
#' @param outFile   Character vector. Path of the output netCDF file.
#' @param inEtp     Character vector. Path to the input evapotranspiration netCDF file.
#' @param title     Character vector. Title of the dataset.
#' @param comment   Character vector. A comment.
#' @param inMask    Character vector. Path to a netCDF mask file.
#' @param block     Integer. Number of latitude blocks to be processed at the
#' same time.
#'
#' @return Computes the SPEI time series and stores it in outFile following
#' the same data structure of inPre.
#'
#' @section References
#'
#' @examples
#' # see script 'analysis.R'
#'
#' @export
spei.nc <- function(sca, inPre, outFile, inEtp=NA, title=NA, comment=NA,
                    inMask=NA, block=18, tlapse=NA) {
  require(SPEI)
  require(ncdf4)
  require(snowfall)
  
	# Read data
	pre.nc <- nc_open(inPre) # nc_open does not work with file connections
	if (!is.na(inEtp)) { # there is ETP input data, so compute the SPEI
		etp.nc <- nc_open(inEtp)
		isSPEI <- TRUE
	} else {
		isSPEI <- FALSE
	}
	
	if (!is.na(inMask)) { # there is a mask
		mask.nc <- nc_open(inMask)
		isMask <- TRUE
	} else {
		isMask <- FALSE
	}

	# Create variable and output file
	if (!is.na(tlapse[1])) {
		spi.dim.tme <- ncdim_def(pre.nc$dim$time$name, pre.nc$dim$time$units,
				pre.nc$dim$time$vals[tlapse[1]:tlapse[2]])
	} else {
		spi.dim.tme <- ncdim_def(pre.nc$dim$time$name, pre.nc$dim$time$units,
				pre.nc$dim$time$vals)
	}
	if (isSPEI) {
		nam <- 'spei'
		longnam <- 'Standardized Precipitation-Evapotranspiration Index'
	} else {
		nam <- 'spi'
		longnam <- 'Standardized Precipitation Index'
	}
	out.nc.var <- ncvar_def(
		name=nam,
		units='z-values',
		dim=list(pre.nc$dim$lon, pre.nc$dim$lat, spi.dim.tme),
		missval=1.e30,
		longname=longnam,
		prec='float',
		compression=9)
	out.nc <- nc_create(outFile, out.nc.var)
	ncatt_put(out.nc,0,'Title',title)
	ncatt_put(out.nc,0,'Version','2.5')
	ncatt_put(out.nc,0,'Id',outFile)
	ncatt_put(out.nc,0,'Summary',paste('Global dataset of the Standardized
	  Precipitation-Evapotranspiration Index (SPEI) at the ', sca,
	  '-month', ifelse(sca==1,'','s'), ' time scale. ', comment, sep=''))
	ncatt_put(out.nc,0,'Keywords','drought, climatology, SPEI, Standardized
	  Precipitation-Evapotranspiration Index')
	ncatt_put(out.nc,0,'Institution','Consejo Superior de Investigaciones
	          Científicas, CSIC')
	ncatt_put(out.nc,0,'Url','http://sac.csic.es/spei')
	ncatt_put(out.nc,0,'Creators','Santiago Beguería <santiago.begueria@csic.es>
	          and Sergio Vicente-Serrano <svicen@ipe.csic.es>')
	ncatt_put(out.nc,0,'Software','Created in R using the SPEI package
	          (https://cran.r-project.org/web/packages/SPEI/;
	          https://github.com/sbegueria/SPEI)')
	b <- match.call()
	ncatt_put(out.nc,0,'Call',
		paste(b[[1]],'(',names(b)[[2]],'=',b[[2]],', ',names(b)[[3]],'=',b[[3]],
		', ',names(b)[[4]],'=',b[[4]],')',sep=''))
	ncatt_put(out.nc,0,'Date',date())
	ncatt_put(out.nc,0,'Reference','Beguería S., Vicente-Serrano S., Reig F., Latorre B. (2014) Standardized precipitation evapotranspiration index (SPEI) revisited: Parameter fitting, evapotranspiration models, tools, datasets and drought monitoring. International Journal of Climatology 34, 3001-3023.')
	ncatt_put(out.nc,0,'Reference2','Vicente-Serrano S.M., Beguería S., López-Moreno J.I. (2010) A Multiscalar Drought Index Sensitive to Global Warming: The Standardized Precipitation Evapotranspiration Index. Journal of Climate 23, 1696–1718.')
	ncatt_put(out.nc,0,'Reference3','Beguería S., Vicente-Serrano S., Angulo-Martínez M. (2010) A multi-scalar global drought data set: the SPEIbase. Bulletin of the American Meteorological Society 91(10), 1351–1356.')
	ncatt_put(out.nc,'lon','long_name','longitude')
	ncatt_put(out.nc,'lon','limits',paste(out.nc$dim$lon$vals[1],
		out.nc$dim$lon$vals[out.nc$dim$lon$len],sep=', '))
	ncatt_put(out.nc,'lon','convention','Coordinates are referred to cell centres')
	ncatt_put(out.nc,'lat','long_name','latitude')
	ncatt_put(out.nc,'lat','limits',paste(out.nc$dim$lat$vals[1],
		out.nc$dim$lat$vals[out.nc$dim$lat$len],sep=', '))
	ncatt_put(out.nc,'lat','convention','Coordinates are referred to cell centres')
	ncatt_put(out.nc,'time','limits',paste(out.nc$dim$time$vals[1],
		out.nc$dim$time$vals[out.nc$dim$lat$len],sep=', '))

	# Dimensions
	dimlon <- out.nc.var$dim[[1]]$len
	dimlat <- out.nc.var$dim[[2]]$len
	dimtme <- out.nc.var$dim[[3]]$len
	
	# A vector with the number of days per month
	if (isSPEI) {
	  ndays <- Hmisc::monthDays(as.Date(out.nc.var$dim[[3]]$vals,
	                                    origin='1900-1-1'))
	}
	
	n <- block # latitudes block (number of lines read each time)
 	for (j in seq(dimlat,1,-n)-n+1) { # j=265 (lat=42.25º), lon=361 (lon=0.25)
		start <- c(1,j,ifelse(is.na(tlapse[1]),1,tlapse[1]))
		count <- c(dimlon,n,dimtme)
		if (isSPEI) {
		  # Read precipitation and ET data and calculate balance
		  x <- ncvar_get(etp.nc, varid=etp.nc$var[[1]]$name, start, count)
      x <- ncvar_get(pre.nc, varid=pre.nc$var[[1]]$name, start, count) - x
		} else {
		  # Read precipitation data
			x <- ncvar_get(pre.nc, varid=pre.nc$var[[1]]$name, start, count)
		}
		if (isMask) { # use a mask file
			msk <- ncvar_get(mask.nc, 'mask', c(1,j), c(dimlon,n))
		}
		# Convert from matrix to list
		x.list <- vector('list', dimlon*n)
		for (i in 1:dimlon) {
			for (h in 1:n-1) {
				if (!isMask) {
					x.list[[h*dimlon+i]] <- x[i,{h+1},]
				} else {
					if (msk[i,{h+1}]==1) {
						x.list[[h*dimlon+i]] <- x[i,{h+1},]
					} else {
						x.list[[h*dimlon+i]] <- rep(NA,dimtme)
					}
				}
			}
		}

		# Compute series of anomalies
		speii <- function(d,s) {spei(d, s, na.rm=TRUE)$fitted}
		x.list <- sfLapply(x.list, speii, sca)
		# Convert back to matrix
		for (i in 1:dimlon) {
			for (h in 1:n-1) {
				part <- x.list[[h*dimlon+i]]
				part[is.nan(part/part)] <- NA
				x[i,{h+1},] <- part
			}
		}
		
		# Store into output netCDF file
		ncvar_put(out.nc, nam, x, c(1,j,1), count)
		nc_sync(out.nc)
	} # next j

	# close output netCDF file
	nc_close(out.nc)
}







