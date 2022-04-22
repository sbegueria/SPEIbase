#' Compute SPEI from very large arrayed data sets in netCDF format.
#' @param sca       Integer. The time scale of the SPEI.
#' @param inPre     Character vector. Path to the input precipitation netCDF file.
#' @param outFile   Character vector. Path of the output netCDF file.
#' @param inEtp     Character vector. Path to the input evapotranspiration netCDF file.
#' @param title     Character vector. Title of the dataset.
#' @param comment   Character vector. A comment.
#' @param inMask    Character vector. Path to a netCDF mask file.
#' @param block     Integer. Number of latitude blocks to be processed at the
#' same time. Must be an integer dividend of 360.
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
  #require(SPEI)
  require(ncdf4)
  require(snowfall)
  require(Hmisc)
  
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
	ncatt_put(out.nc,0,'Version','2.6')
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
	ncatt_put(out.nc,0,'Reference2','Vicente-Serrano S.M., Beguería S., López-Moreno J.I. (2010) A Multiscalar Drought Index Sensitive to Global Warming: The Standardized Precipitation Evapotranspiration Index. Journal of Climate 23, 1696–1718.')
	ncatt_put(out.nc,0,'Reference3','Beguería S., Vicente-Serrano S., Angulo-Martínez M. (2010) A multi-scalar global drought data set: the SPEIbase. Bulletin of the American Meteorological Society 91(10), 1351–1356.')
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
		  x <- ncvar_get(etp.nc, varid=etp.nc$var[[1]]$name, start, count) * ndays
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



spei <- function(data, scale, kernel=list(type='rectangular',shift=0),
                 distribution='log-Logistic', fit='ub-pwm', na.rm=FALSE, 
                 ref.start=NULL, ref.end=NULL, x=FALSE, params=NULL, ...) {
  
  scale <- as.numeric(scale)
  na.rm <- as.logical(na.rm)
  x <- as.logical(x)
  #if (!exists("data",inherits=F) | !exists("scale",inherits=F)) {
  #	stop('Both data and scale must be provided')
  #}
  if (anyNA(data) && na.rm==FALSE) {
    stop('Error: Data must not contain NAs')
  }
  if (!(distribution %in% c('log-Logistic', 'Gamma', 'PearsonIII'))) {
    stop('Distrib must be one of "log-Logistic", "Gamma" or "PearsonIII"')
  }
  if (!(fit %in% c('max-lik', 'ub-pwm', 'pp-pwm'))) {
    stop('Method must be one of "ub-pwm" (default), "pp-pwm" or "max-lik"')
  }
  if ( (!is.null(ref.start) && length(ref.start)!=2) | (!is.null(ref.end) && length(ref.end)!=2) ) {
    stop('Start and end of the reference period must be a numeric vector of length two.')
  }
  
  if (!is.ts(data)) {
    data <- ts(as.matrix(data), frequency = 12)
  } else {
    data <- ts(as.matrix(data), frequency=frequency(data), start=start(data))
  }
  m <- ncol(data)
  fr <- frequency(data)
  
  
  coef = switch(distribution,
                "Gamma" = array(NA,c(2,m,fr),list(par=c('alpha','beta'),colnames(data),NULL)),
                "log-Logistic" = array(NA,c(3,m,fr),list(par=c('xi','alpha','kappa'),colnames(data),NULL)),
                "PearsonIII" = coef <- array(NA,c(3,m,fr),list(par=c('mu','sigma','gamma'),colnames(data),NULL))
  )
  
  dim_one = ifelse(distribution == "Gamma", 2, 3)
  
  if (!is.null(params)) {
    if (dim(params)[1]!=dim_one | dim(params)[2]!=m | dim(params)[3]!=12) {
      stop(paste('parameters array should have dimensions (', dim_one, ', ', m, ', 12)',sep=' '))
    }
  }
  
  # Loop through series (columns in data)
  if (!is.null(ref.start) && !is.null(ref.end)) {
    data.fit <- window(data,ref.start,ref.end)	
  } else {
    data.fit <- data
  }
  std <- data*NA
  for (s in 1:m) {
    # Cumulative series (acu)
    acu <- data.fit[,s]
    acu.pred <- data[,s]
    if (scale>1) {
      wgt <- kern(scale,kernel$type,kernel$shift)
      acu[scale:length(acu)] <- rowSums(embed(acu,scale)*wgt,na.rm=na.rm)
      acu[1:(scale-1)] <- NA
      acu.pred[scale:length(acu.pred)] <- rowSums(embed(acu.pred,scale)*wgt,na.rm=na.rm)
      acu.pred[1:(scale-1)] <- NA
    }
    
    # Loop through the months
    for (c in (1:fr)) {
      # Filter month m, excluding NAs
      f <- which(cycle(acu)==c)
      f <- f[!is.na(acu[f])]
      ff <- which(cycle(acu.pred)==c)
      ff <- ff[!is.na(acu.pred[ff])]
      
      # Monthly series, sorted
      month <- sort.default(acu[f], method="quick")
      
      if (length(month)==0) {
        std[f] <- NA
        next()
      }
      
      if (is.null(params)) {
        month_sd = sd(month,na.rm=TRUE)
        if (is.na(month_sd) || (month_sd == 0)) {
          std[f] <- NA
          next
        }
        
        if(distribution != "log-Logistic"){
          pze <- sum(month==0)/length(month)
          month = month[month > 0]
        }
        
        # Stop early and assign NAs if month's data is length < 4
        if(length(month) < 4){
          std[ff,s] = NA
          coef[,s,c] <- NA
          next
        }
        
        # Calculate probability weighted moments based on fit with lmomco or TLMoments
        pwm = switch(fit,
                     "pp-pwm" = pwm.pp(month,-0.35,0, nmom=3),
                     #pwm.ub(month, nmom=3)
                     TLMoments::PWM(month, order=0:2)
        )
        
        # Check L-moments validity
        lmom <- pwm2lmom(pwm)
        if ( !are.lmom.valid(lmom) || anyNA(lmom[[1]]) || any(is.nan(lmom[[1]])) ){
          next
        }
        
        # lmom fortran functions need specific inputs L1, L2, T3
        # this is handled by lmomco internally with lmorph
        fortran_vec = c(lmom$lambdas[1:2], lmom$ratios[3])
        
        # Calculate parameters based on distribution with lmom then lmomco
        f_params = switch(distribution,
                          "log-Logistic" = tryCatch(lmom::pelglo(fortran_vec), error = function(e){ parglo(lmom)$para }),
                          "Gamma" = tryCatch(lmom::pelgam(fortran_vec), error = function(e){ pargam(lmom)$para }),
                          "PearsonIII" = tryCatch(lmom::pelpe3(fortran_vec), error = function(e){ parpe3(lmom)$para })
        )
        
        # Adjust if user chose log-Logistic and max-lik
        if(distribution == 'log-Logistic' && fit=='max-lik'){
          f_params = parglo.maxlik(month, f_params)$para
        }
      } else {
        
        f_params = as.vector(params[,s,c])
        
      }
      
      # Calculate cdf based on distribution with lmom
      cdf_res = switch(distribution,
                       "log-Logistic" = lmom::cdfglo(acu.pred[ff], f_params),
                       "Gamma" = lmom::cdfgam(acu.pred[ff], f_params),
                       "PearsonIII" = lmom::cdfpe3(acu.pred[ff], f_params)				  				
      )
      
      std[ff,s] = qnorm(cdf_res)
      coef[,s,c] <- f_params
      
      # Adjust if user chose Gamma or PearsonIII
      if(distribution != 'log-Logistic'){ 
        std[ff,s] = qnorm(pze + (1-pze)*pnorm(std[ff,s]))
      }
      
    } # next c (month)
  } # next s (series)
  colnames(std) <- colnames(data)
  
  z <- list(call=match.call(expand.dots=FALSE),
            fitted=std,coefficients=coef,scale=scale,kernel=list(type=kernel$type,
                                                                 shift=kernel$shift,values=kern(scale,kernel$type,kernel$shift)),
            distribution=distribution,fit=fit,na.action=na.rm)
  if (x) z$data <- data
  if (!is.null(ref.start)) z$ref.period <- rbind(ref.start,ref.end)
  
  class(z) <- 'spei'
  return(z)
}


spei_OLD <- function (data, scale, kernel = list(type = "rectangular", shift = 0),
          distribution = "log-Logistic", fit = "ub-pwm", na.rm = FALSE,
          ref.start = NULL, ref.end = NULL, x = FALSE, params = NULL,
          ...)
{
  scale <- as.numeric(scale)
  na.rm <- as.logical(na.rm)
  x <- as.logical(x)
  if (sum(is.na(data)) > 0 & na.rm == FALSE) {
    stop("Error: Data must not contain NAs")
  }
  if (distribution != "log-Logistic" & distribution != "Gamma" &
      distribution != "PearsonIII") {
    stop("Distrib must be one of \"log-Logistic\", \"Gamma\" or \"PearsonIII\"")
  }
  if (fit != "max-lik" & fit != "ub-pwm" & fit != "pp-pwm") {
    stop("Method must be one of \"ub-pwm\" (default), \"pp-pwm\" or \"max-lik\"")
  }
  if ({
    !is.null(ref.start) & length(ref.start) != 2
  } | {
    !is.null(ref.end) & length(ref.end) != 2
  }) {
    stop("Start and end of the reference period must be a numeric vector of length two.")
  }
  if (!is.ts(data)) {
    data <- ts(as.matrix(data), frequency = 12)
  }
  else {
    data <- ts(as.matrix(data), frequency = frequency(data),
               start = start(data))
  }
  m <- ncol(data)
  fr <- frequency(data)
  if (distribution == "Gamma") {
    coef <- array(NA, c(2, m, fr), list(par = c("alpha",
                                                "beta"), colnames(data), NULL))
    if (!is.null(params)) {
      if (dim(params)[1] != 2 | dim(params)[2] != m | dim(params)[3] !=
          12) {
        stop(paste("parameters array should have dimensions (2,",
                   m, "12)", sep = " "))
      }
    }
  }
  if (distribution == "log-Logistic") {
    coef <- array(NA, c(3, m, fr), list(par = c("xi", "alpha",
                                                "kappa"), colnames(data), NULL))
    if (!is.null(params)) {
      if (dim(params)[1] != 3 | dim(params)[2] != m | dim(params)[3] !=
          12) {
        stop(paste("parameters array should have dimensions (3,",
                   m, "12)", sep = " "))
      }
    }
  }
  if (distribution == "PearsonIII") {
    coef <- array(NA, c(3, m, fr), list(par = c("mu", "sigma",
                                                "gamma"), colnames(data), NULL))
    if (!is.null(params)) {
      if (dim(params)[1] != 3 | dim(params)[2] != m | dim(params)[3] !=
          12) {
        stop(paste("parameters array should have dimensions (3,",
                   m, "12)", sep = " "))
      }
    }
  }
  if (!is.null(ref.start) & !is.null(ref.end)) {
    data.fit <- window(data, ref.start, ref.end)
  }
  else {
    data.fit <- data
  }
  std <- data * NA
  for (s in 1:m) {
    acu <- data.fit[, s]
    acu.pred <- data[, s]
    if (scale > 1) {
      wgt <- kern(scale, kernel$type, kernel$shift)
      acu[scale:length(acu)] <- rowSums(embed(acu, scale) *
                                          wgt, na.rm = na.rm)
      acu[1:{
        scale - 1
      }] <- NA
      acu.pred[scale:length(acu.pred)] <- rowSums(embed(acu.pred,
                                                        scale) * wgt, na.rm = na.rm)
      acu.pred[1:{
        scale - 1
      }] <- NA
    }
    for (c in (1:fr)) {
      f <- which(cycle(acu) == c)
      f <- f[!is.na(acu[f])]
      ff <- which(cycle(acu.pred) == c)
      ff <- ff[!is.na(acu.pred[ff])]
      month <- sort(acu[f])
      if (length(month) == 0) {
        std[f] <- NA
        (next)()
      }
      if (is.null(params)) {
        if (is.na(sd(month, na.rm = TRUE)) | (sd(month,
                                                 na.rm = TRUE) == 0)) {
          std[f] <- NA
          (next)()
        }
        if (distribution == "log-Logistic") {
          if (fit == "pp-pwm") {
            pwm <- pwm.pp(month, -0.35, 0)
          }
          else {
            pwm <- (month)
          }
          lmom <- pwm2lmom(pwm)
          if (!are.lmom.valid(lmom) | is.na(sum(lmom[[1]])) |
              is.nan(sum(lmom[[1]]))) {
            (next)()
          }
          llpar <- parglo(lmom)
          if (fit == "max-lik") {
            llpar <- parglo.maxlik(month, llpar$para)
          }
          std[ff, s] <- qnorm(pglo(acu.pred[ff], llpar))
          coef[, s, c] <- llpar$para
        } else if (distribution == "Gamma" | distribution ==
                 "PearsonIII") {
          zeros <- sum(month == 0)
          pze <- sum(month == 0)/length(month)
          if (fit == "pp-pwm") {
            pwm <- pwm.pp(month[month > 0], -0.35, 0)
          }
          else {
            pwm <- pwm.ub(month[month > 0])
          }
          lmom <- pwm2lmom(pwm)
          if (!are.lmom.valid(lmom) | is.na(sum(lmom[[1]])) |
              is.nan(sum(lmom[[1]]))) {
            (next)()
          }
          if (distribution == "Gamma") {
            gampar <- pargam(lmom)
            std[ff, s] <- qnorm(cdfgam(acu.pred[ff],
                                       gampar))
            std[ff, s] <- qnorm(pze + (1 - pze) * pnorm(std[ff,
                                                            s]))
            coef[, s, c] <- gampar$para
          }
          else if (distribution == "PearsonIII") {
            p3par <- parpe3(lmom)
            std[ff, s] <- qnorm(cdfpe3(acu.pred[ff],
                                       p3par))
            std[ff, s] <- qnorm(pze + (1 - pze) * pnorm(std[ff,
                                                            s]))
            coef[, s, c] <- p3par$para
          }
        }
      } else {
        if (dim(params)[1] != 3 & dim(params)[2] != m &
            dim(params)[3] != 12) {
          stop(paste("params should be an array with dimensions (3,",
                     m, ",12)", sep = " "))
        }
        coef[, s, c] <- params[, s, c]
        if (distribution == "log-Logistic") {
          std[ff, s] <- qnorm(pglo(acu.pred[ff], list(type = "glo",
                                                      para = params[, s, c], source = "user")))
        } else {
          if (distribution == "Gamma") {
            std[ff, s] <- qnorm(cdfgam(acu.pred[ff],
                                       list(type = "gam", para = params[, s, c],
                                            source = "user")))
            std[ff, s] <- qnorm(pze + (1 - pze) * pnorm(std[ff,
                                                            s]))
          }
          else if (distribution == "PearsonIII") {
            std[ff, s] <- qnorm(cdfpe3(acu.pred[ff],
                                       list(type = "pe3", para = params[, s, c],
                                            source = "user")))
            std[ff, s] <- qnorm(pze + (1 - pze) * pnorm(std[ff,
                                                            s]))
          }
        }
      }
    }
  }
  colnames(std) <- colnames(data)
  z <- list(call = match.call(expand.dots = FALSE), fitted = std,
            coefficients = coef, scale = scale, kernel = list(type = kernel$type,
                                                              shift = kernel$shift, values = kern(scale, kernel$type,
                                                                                                  kernel$shift)), distribution = distribution,
            fit = fit, na.action = na.rm)
  if (x)
    z$data <- data
  if (!is.null(ref.start))
    z$ref.period <- rbind(ref.start, ref.end)
  class(z) <- "spei"
  return(z)
}




