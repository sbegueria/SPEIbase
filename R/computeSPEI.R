# Compute SPEI at various time scales from very large netCDF files,
# using parallel computing.

source('./R/functions.R')
library(SPEI)
library(ncdf4)
library(snowfall)

sfInit(parallel=TRUE, cpus=10)
sfExport(list='spei', namespace='SPEI')

# Compute all time scales between 1 and 48
for (i in c(1:48)) {
#for (i in c(1,3,6,12)) {
    spei.nc(
		sca=i,
		inPre='./inputData/cru_ts3.24.01.1901.2015.pre.dat.nc',
		inEtp='./inputData/cru_ts3.24.01.1901.2015.pet.dat.nc',
		outFile=paste('./outputNcdf/spei',formatC(i, width=2, format='d', flag='0'),'.nc',sep=''),
		title=paste('Global ',i,'-month',ifelse(i==1,'','s'),' SPEI, z-values, 0.5 degree',sep=''),
		comment='Using CRU TS 4.00 precipitation and potential evapotranspiration data',
		block=24,
		inMask=NA,
		tlapse=NA
	)
  gc()
}

sfStop()







