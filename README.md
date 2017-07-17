# SPEIbase

R code used in generating the SPEI global database,
http://spei.csic.es/database.html.

In order to run the script (R/computeSPEI.R), you need to replace the fake
files in the 'inputData' directory with the real ones containing the data.
These can be downloaded from the website of the Climatic Research Unit (CRU),
University of East Anglia.

The output files, named 'SPEI01.nc', 'SPEI02.nc', etc, will be stored in the
'outputNcdf' directory.

The output files are in netCDF v4 format, which allows for data compression.
It is possible to convert the files to the netCDF v3 format used in many
GIS softwares, by using the nccopy program by unidata.
For instance, to convert a netCDF-4 format file foo4c.nc to file foo3.nc in
netCDF-3 you can use:

nccopy -k classic foo4c.nc foo3.nc.

Another way would be to install the netCDF operators (NCO) toolset from unidata,
and then use the ncks utility:

ncks -3 foo4c.nc foo3.nc.

In fact, ncks allows for much greater functionality.
For example, if one wants to extract the first 100 times from an SPEIbase file:

ncks -d time,0,100 spei_12.nc output_file.nc

would generate a (smaller) netCDF file with only those timesteps.
In a similar fashion, it is possible to use ncks to select a specific
geographical region, and other useful options.
