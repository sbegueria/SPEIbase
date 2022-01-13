# Compute SPEI at various time scales from very large netCDF files,
# using parallel computing. The resulting netCDF files are stored
# to disk on directory './outputNcdf/'.

source('./R/functions.R')
library(SPEI)
library(ncdf4)
library(snowfall)
library(Hmisc)

# Init a parallel computing cluster; modify the parameter `cpus` to the
# desired number of cores; otherwise, use all available cores
sfInit(parallel=TRUE, cpus=detectCores())
sfExport(list='spei', namespace='SPEI')

# Compute SPEI at all time scales between 1 and 48 and store to disk
for (i in c(1:48)) {
#for (i in c(1,3,6,12)) {
    spei.nc(
      sca=i,
		  inPre='./inputData/cru_ts4.03.1901.2018.pre.dat.nc',
		  inEtp='./inputData/cru_ts4.03.1901.2018.pet.dat.nc',
		  outFile=paste('./outputNcdf/spei',
		              formatC(i, width=2, format='d', flag='0'),'.nc',sep=''),
		  title=paste('Global ',i,'-month',
		            ifelse(i==1,'','s'),' SPEI, z-values, 0.5 degree',sep=''),
		  comment='Using CRU TS 4.03 precipitation and potential evapotranspiration data',
		  block=36,
		  inMask=NA,
		  tlapse=NA
	  )
  gc()
}

# Stop the parallel cluster
sfStop()







