# Compute SPEI at various time scales from very large netCDF files,
# using parallel computing. The resulting netCDF files are stored
# to disk on directory './outputNcdf/'.

if (!require('pacman')) install.packages('pacman')
pacman::p_load(SPEI, ncdf4, snowfall, parallel, Hmisc)
source('./R/functions.R')

# Init a parallel computing cluster; modify the parameter `cpus` to the
# desired number of cores; otherwise, use all available cores
#sfInit(parallel=TRUE, cpus=detectCores())
sfInit(parallel=TRUE, cpus=4)
sfExport(list='spei', namespace='SPEI')

# Compute SPEI at all time scales between 1 and 48 and store to disk
#for (i in c(1:48)) {
for (i in c(1,3,6,12)) {
    spei.nc(
      sca=i,
		  inPre='./inputData/cru_ts4.05.1901.2020.pre.dat.nc',
		  inEtp='./inputData/cru_ts4.05.1901.2020.pet.dat.nc',
		  outFile=paste('./outputNcdf/spei',
		              formatC(i, width=2, format='d', flag='0'),'.nc',sep=''),
		  title=paste('Global ',i,'-month',
		            ifelse(i==1,'','s'),' SPEI, z-values, 0.5 degree',sep=''),
		  comment='Using CRU TS 4.05 precipitation and potential evapotranspiration data',
		  block=36,
		  inMask=NA,
		  tlapse=NA
	  )
  gc()
}

# Stop the parallel cluster
sfStop()
