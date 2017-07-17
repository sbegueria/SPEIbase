# -----------------------------
# FUNCTIONS
# -----------------------------

spei.nc <- function(outFile, sca, inPre, inEtp=NA, inMask=NA, method=NA, nCPU=1, append=F,
	block=18, tlapse=NA) {

	require(SPEI)
	require(ncdf4)
#	require(snowfall)

	# Leer datos
	pre.nc <- nc_open(inPre)
	if (!is.na(inEtp)) { # there is ETP input data, so compute the SPEI
		etp.nc <- nc_open(inEtp)
		isSPEI <- T
	} else {
		isSPEI <- F
	}
	if (!is.na(inMask)) { # there is a mask
		mask.nc <- nc_open(inMask)
		isMask <- T
	} else {
		isMask <- F
	}

	# Crea la variable y archivo de salida
	if (append) {
		spi.nc <- nc_open(outFile, write=T)
		nam <- 'spei'
	} else {
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
		spi.nc.var <- ncvar_def(
			name=nam,
			units='z-values',
			dim=list(pre.nc$dim$lon, pre.nc$dim$lat, spi.dim.tme),
			missval=1.e30,
			longname=longnam,
			prec='float')
#			prec='float',
#			compression=9)
		#spi.nc <- nc_create(paste('./spei', sca, method, '.nc', sep=''), spi.nc.var)
		spi.nc <- nc_create(outFile, spi.nc.var)
		ncatt_put(spi.nc, nam, attname='description',
			attval=paste('Using ',method,' method for parameter estimation',sep=''))
		ncatt_put(spi.nc, nam, attname='version', attval='1.0')
		ncatt_put(spi.nc, nam, attname='creator',
			attval='Santiago Begueria <santiago.begueria@csic.es>')
		ncatt_put(spi.nc, 'lon', attname='limits',			attval=paste(spi.nc$dim$lon$vals[1],
			spi.nc$dim$lon$vals[spi.nc$dim$lon$len],sep=', '))
		ncatt_put(spi.nc, 'lat', attname='limits',			attval=paste(spi.nc$dim$lat$vals[1],
			spi.nc$dim$lat$vals[spi.nc$dim$lat$len],sep=', '))
		ncatt_put(spi.nc, 'time', attname='limits',			attval=paste(spi.nc$dim$time$vals[1],
			spi.nc$dim$time$vals[spi.nc$dim$time$len],sep=', '))
	}
	#print.ncdf4(spi.nc)

	# Dimensiones
	dimlon <- spi.nc$dim$lon$len
	dimlat <- spi.nc$dim$lat$len
	dimtme <- spi.nc$dim$time$len

	# Inicia ejecucion en paralelo
#	sfInit(parallel=ifelse(nCPU==1,F,T), cpus=nCPU)

	n <- block # latitudes block (number of lines read each time)
 	for (j in seq(dimlat,1,-n)-n+1) { # j=265 (lat=42.25º), lon=361 (lon=0.25)

		# Read precipitation and ET time series and calculate balance
		start <- c(1,j,ifelse(is.na(tlapse[1]),1,tlapse[1]))
		count <- c(dimlon,n,dimtme)
		x <- ncvar_get(pre.nc, varid=pre.nc$var[[1]]$name, start, count) -
			ncvar_get(etp.nc, varid=etp.nc$var[[1]]$name, start, count)
		#x <- ifelse(x==-9.99000e+01,NA,x)
		if (isMask) { # use a mask file
			msk <- ncvar_get(mask.nc, 'mask', c(1,j), c(dimlon,n))
		}
		# convert from matrix to a list
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

		# Compute SPI series for cell (x,y)
		#for (hh in 1:length(x.list)) kk <- spi(x.list[[hh]], sca, 12, method)
		#kk <- spi(x.list[[hh]], sca, 12, method)
		#kk <- spi(x, sca, 12)
		speii <- function(d,s) {require(SPEI); spei(d,s,na.rm=TRUE)}
		x.list <- sfLapply(x.list,speii,sca)

		# Store the SPEI series in the netCDF file
		# convert back to matrix
		for (i in 1:dimlon) {
			for (h in 1:n-1) {
				part <- x.list[[h*dimlon+i]]
				part[is.nan(part/part)] <- NA
				x[i,{h+1},] <- part
			}
		}
		# store into ncdf file
		ncvar_put(spi.nc, nam, x, c(1,j,1), count)
		#nc_sync(spi.nc)

	} # next j

	# close output ncdf
	nc_close(spi.nc)
	
	# stop cluster
#	sfStop()
}

# -----------------------------
# ANALYSIS
# -----------------------------

inPre <- './inputData/cru_ts_3_10.1901.2009.pre.dat.nc'
inTho <- './inputData/etpThornthwaite.nc'
inHar <- './inputData/etpHargreaves.nc'
inPen <- './inputData/etpPenman.nc'
#inMask <- './mask.nc'

library(snowfall)
sfInit(parallel=TRUE, cpus=10)
for (i in 1:48) {
	spei.nc(paste('SPEItho',i,'.nc',sep=''), i, inPre, inTho, nCPU=10, block=12)
	spei.nc(paste('SPEIhar',i,'.nc',sep=''), i, inPre, inHar, nCPU=10, block=12)
	spei.nc(paste('SPEIpen',i,'.nc',sep=''), i, inPre, inPen, nCPU=10, block=12)
}
sfStop()
