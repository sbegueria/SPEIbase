# SPEIbase [![DOI](https://zenodo.org/badge/97464490.svg)](https://zenodo.org/badge/latestdoi/97464490)

R code used in generating the SPEI global database.


## Details

The SPEI global database, http://spei.csic.es/database.html, offers long-time,
robust information about drought conditions at the global scale, with a 0.5 
degrees spatial resolution and a monthly time resolution.
It has a multi-scale character, providing SPEI time-scales between 1 and 48 
months. 
Currently it covers the period between January 1901 and December 2023.
The SPEI is the Standardized Precipitation-Evapotranspiration Index, defined
in the following research papers:

* **Vicente-Serrano S.M., Santiago Beguería, Juan I. López-Moreno, (2010) A 
Multi-scalar drought index sensitive to global warming: The Standardized 
Precipitation Evapotranspiration Index - SPEI. Journal of Climate 23: 
1696-1718.**
* **Beguería, S., Vicente-Serrano, S.M., Fergus Reig, Borja Latorre. 
Standardized Precipitation Evapotranspiration Index (SPEI) revisited (2014): 
parameter fitting, evapotranspiration models, kernel weighting, tools, 
datasets and drought monitoring. International Journal of Climatology, 
34: 3001-3023**

The data are accessible from the SPEI web site, http://spei.csic.es, but we
have shared the R code on this repository to allow anyone to reproduce the
data set.

If you use this code or parts of it, it'd be nice if you cite the repo as:
* **Beguería S. (2017) SPEIbase: R code used in generating the SPEI global 
database, doi:10.5281/zenodo.834462**.

If you use the SPEI data set in your research, please cite the following papers:
* **Beguería, S., Vicente-Serrano, S.M. y Angulo, M., (2010): A multi-scalar 
global drought data set: the SPEIbase: A new gridded product for the analysis 
of drought variability and impacts. Bulletin of the American Meteorological 
Society. 91, 1351-1354**
* **Vicente-Serrano, S.M., Beguería, S., López-Moreno, J.I., Angulo, M., El 
Kenawy, A. (2010): A new global 0.5° gridded dataset (1901-2006) of a 
multiscalar drought index: comparison with current drought index datasets 
based on the Palmer Drought Severity Index. Journal of Hydrometeorology. 
11: 1033-1043**


## Using the code

The script `R/computeSPEI.R` computes the global SPEI data set at different
time scales. One netCDF file covering the whole globe and time period is
generated for each time scale, e.g. `spei01.nc` for a time scale of 1 month,
etc. Output files are stored on `/outputNcdf`.
These are the global files that can be downloaded from http://spei.csic.es/database.html.

Before running the script it is necessary to replace the fake files in the
`/inputData` directory with the real ones containing the data, which are:
`cru_ts4.08.1901.2023.pet.dat.nc` and
`cru_ts4.08.1901.2023.pre.dat.nc`.
These files can be downloaded from the website of the Climatic Research Unit
(CRU) of the University of East Anglia.
Check http://www.cru.uea.ac.uk/data/ for more information on this data set.
You might need to decompress the files before computing the SPEI.

The script `R/outputTxt` generates, from the .nc files computed previously,
additional files containing the SPEI time series at the scales between 1 and
48 months, for every single cell in the data set.
These files are stored in plain text, comma-separated, files, and stored
under `/outputTxt`.
The files names (e.g. `spei_-0.25_5.25.csv.gz`) contain the central coordinates
of the grid cell.
These are the files that can be downloaded as 'single location' data at
http://spei.csic.es/database.html.



## Using the dataset

The output files are in netCDF v4 format, which allows for data compression.
The files can be read in R using the ncdf4 package.

However, some widely used GIS packages do not provide suport yet for netCDF v4.
It is possible to convert the files to the netCDF v3 format used in many
GIS softwares, by using the nccopy program by unidata.
For instance, to convert a netCDF-4 format file foo4c.nc to file foo3.nc in
netCDF-3 you can use:

```
nccopy -k classic foo4c.nc foo3.nc.
```

Another way would be to install the netCDF operators (NCO) toolset from unidata,
and then use the ncks utility:

```
ncks -3 foo4c.nc foo3.nc.
```

In fact, ncks allows for much greater functionality.
For example, if one wants to extract the first 100 times from an SPEIbase file:

```
ncks -d time,0,100 spei_12.nc output_file.nc
```

would generate a (smaller) netCDF file with only those timesteps.
In a similar fashion, it is possible to use ncks to select a specific
geographical region, and other useful options.


## Version history

* SPEIbase v2.10.0: 1) Based on the CRU TS 4.08 dataset, spanning the period
between January 1901 to December 2023. 2) Minor changes to the CRS definition
metadata (crs 4979 to 4326). 3) Using SPEI package version 1.8.1.
* SPEIbase v2.9.0: 1) Based on the CRU TS 4.07 dataset, spanning the period
between January 1901 to December 2022.
* SPEIbase v2.8.0: 1) Based on the CRU TS 4.06 dataset, spanning the period
between January 1901 to December 2021. 2) Using SPEI package version 1.8.0.
* SPEIbase v2.7: 1) Based on the CRU TS 4.05 dataset, spanning the period
between January 1901 to December 2020. 2) Using TLMoments::PWM instead of
lmomco::pwm.ub for calculating distribution parameters.
* SPEIbase v2.6: 1) Based on the CRU TS 4.03 dataset, spanning the period
between January 1901 to December 2018.
* SPEIbase v2.5: 1) Based on the CRU TS 3.24.01 dataset, extending the temporal
range of the SPEIbase up to December 2015. 2) Corrected an important bug on
versions 2.2 to 2.4 of the dataset that prevented correctly reading the ETo data
in mm/month.
* SPEIbase v2.4: 1) Based on the CRU TS 3.23 dataset, extending the temporal 
range of the SPEIbase up to December 2014.
* SPEIbase v2.3: 1) Based on the CRU TS 3.22 dataset, extending the temporal 
range of the SPEIbase up to December 2013.
* SPEIbase v2.2: 1) The CRU TS 3.2 dataset has been used, extending the data 
range of the SPEIbase up to December 2011. 2) Potential evapotranspiration data 
from the CRU TS 3.2 dataset has been used instead of computing our own.
* SPEIbase v2.1: 1) The revised dataset CRU TS 3.10.01 for precipitation is 
used. 2) Data on surface pressure and wind from 20th Century Reanalysis v.2 
until December 2009 has been used for computing Penman PET. 3) An error while 
reading some data sources that yielded no data at longitudes 179.25 and 179.75 
has been corrected.
* SPEIbase v2.0: 1) Data has been extended to the period 1901-2010. 2) The 
FAO-56 Penman-Monteith's method has been used for computing PET instead of 
Thornthwaite in SPEIbase v1.0. 3) Unbiased probability weighted moments (ub-pwm) 
method has been used for fitting the log-Logistic distribution, instead of the 
sub-optimal plotting-position pwm method used in version 1.0. 4) The whole 
World is put in one single netCDF file.
* SPEIbase v1.0: The global gridded SPEI dataset is available at time scales 
between 1 and 48 months, with spatial resolution of 0.5º lat/lon and temporal 
coverage between January 1901 and December 2006.


## Got problems?

Feel free to [write an issue](https://github.com/sbegueria/SPEIbase/issues)
if you have any questions or problems.


## Copyright and license

The code on this repository is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, using version 3 or any upgrade
(https://www.gnu.org/licenses/gpl-3.0.en.html).
