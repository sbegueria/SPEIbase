# Compute SPEI at various time scales from very large netCDF files,
# using parallel computing. The resulting netCDF files are stored
# to disk on directory './outputNcdf/'.

if (!require('pacman')) install.packages('pacman')
pacman::p_load(ncdf4, snowfall, parallel, Hmisc, devtools)

# Using SPEI package version 1.8.0. Please, note that this will replace any
# other version of the SPEI package that you might have installed!
devtools::install_github('sbegueria/SPEI@v1.8.0')

# A function to efficiently compute the SPEI over a large netCDF file (using
# multiple cores).
source('./functions.R')

# Initialize a parallel computing cluster; modify the parameter `cpus` to the
# desired number of cores; otherwise, it will use all available cores.
sfInit(parallel=TRUE, cpus=detectCores())
sfExport(list='spei', namespace='SPEI')

# Compute SPEI at all time scales between 1 and 48 and store to disk.
for (i in c(1:48)) {#  for (i in c(6)) {
    spei.nc(
      sca=i,
		  inPre='./inputData/cru_ts4.07.1901.2022.pre.dat.nc',
		  inEtp='./inputData/cru_ts4.07.1901.2022.pet.dat.nc',
		  outFile=paste('./outputNcdf/spei',
		              formatC(i, width=2, format='d', flag='0'),'.nc',sep=''),
		  title=paste('Global ',i,'-month',
		            ifelse(i==1,'','s'),' SPEI, z-values, 0.5 degree',sep=''),
		  comment='Using CRU TS 4.07 precipitation and potential evapotranspiration data',
	    	  version='2.9.0',
		  block=36,
		  inMask=NA,
		  tlapse=NA
	  )
  gc()
}

# Stop the parallel cluster
sfStop()

# Print session information
sessionInfo()
