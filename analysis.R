# -----------------------------
# FUNCTIONS
# -----------------------------

spei.nc <- function(sca, inPre, inEtp, outFile, title, comment=NA,
                    inMask=NA, block=18, tlapse=NA) {
  
  require(SPEI)
  require(ncdf4)
  require(snowfall)
  
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
  
  # Create the variable and the output file
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
  ncatt_put(out.nc,0,'Version','2.2')
  ncatt_put(out.nc,0,'Id',outFile)
  ncatt_put(out.nc,0,'Summary',paste('Global dataset of the Standardized Precipitation-Evapotranspiration Index (SPEI) at the ',
                                     sca,'-month',ifelse(sca==1,'','s'),' time scale. ',comment,sep=''))
  ncatt_put(out.nc,0,'Keywords','drought, climatology, SPEI, Standardized Precipitation-Evapotranspiration Index')
  ncatt_put(out.nc,0,'Institution','Consejo Superior de Investigaciones Científicas, CSIC')
  ncatt_put(out.nc,0,'Url','http://sac.csic.es/spei')
  ncatt_put(out.nc,0,'Creator','Santiago Beguería <santiago.begueria@csic.es>')
  ncatt_put(out.nc,0,'Software','Created in R using the SPEI package (http://cran.r-project.org/web/packages/SPEI)')
  b <- match.call()
  ncatt_put(out.nc,0,'Call',
            paste(b[[1]],'(',names(b)[[2]],'=',b[[2]],', ',names(b)[[3]],'=',b[[3]],
                  ', ',names(b)[[4]],'=',b[[4]],')',sep=''))
  ncatt_put(out.nc,0,'Date',date())
  ncatt_put(out.nc,0,'Reference','Vicente-Serrano S.M., Beguería S., López-Moreno J.I. (2010) A Multiscalar Drought Index Sensitive to Global Warming: The Standardized Precipitation Evapotranspiration Index. Journal of Climate 23, 1696–1718.')
  ncatt_put(out.nc,0,'Reference2','Beguería S., Vicente-Serrano S., Angulo-Martínez M. (2010) A multi-scalar global drought data set: the SPEIbase. Bulletin of the American Meteorological Society 91(10), 1351–1356.')
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
  dimlon <- out.nc$dim$lon$len
  dimlat <- out.nc$dim$lat$len
  dimtme <- out.nc$dim$time$len
  
  n <- block # latitudes block (number of lines read each time)
  for (j in seq(dimlat,1,-n)-n+1) { # j=265 (lat=42.25º), lon=361 (lon=0.25)
    # Read precipitation and ET time series and calculate balance
    start <- c(1,j,ifelse(is.na(tlapse[1]),1,tlapse[1]))
    count <- c(dimlon,n,dimtme)
    if (isSPEI) {
      x <- ncvar_get(pre.nc, varid=pre.nc$var[[1]]$name, start, count) -
        ncvar_get(etp.nc, varid=etp.nc$var[[1]]$name, start, count)
    } else {
      x <- ncvar_get(pre.nc, varid=pre.nc$var[[1]]$name, start, count)
    }
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
    speii <- function(d,s) {require(SPEI); spei(d,s,na.rm=TRUE)$fitted}
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
    ncvar_put(out.nc, nam, x, c(1,j,1), count)
    nc_sync(out.nc)
  } # next j
  
  # close output ncdf
  nc_close(out.nc)
}



# ---------
# ANALYSIS
# ---------

inPre <- './inputData/cru_ts3.20.1901.2011.pre.dat.nc'
inPen <- './inputData/etpPenman/etpPenmanCRUTS3.2.nc'
#inMask <- './mask.nc'

library(snowfall)
sfInit(parallel=TRUE, cpus=24)
for (i in c(1:48)) {
  tit <- paste('Global ',i,'-month',ifelse(i==1,'','s'),' SPEI, z-values, 0.5 degree',sep='')
  spei.nc(i,
          './inputData/cru_ts3.20.1901.2011.pre.dat.nc',
          './inputData/etpPenman/etpPenmanCRUTS3.2.nc',
          paste('./spei',i,'.nc',sep=''),
          tit, 'Using CRU 3.2 precipitation and potential evapotranspiration data')
}
sfStop()
