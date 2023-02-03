library(rugarch)
##########################################################################
##
##	rate()
## 	Calculate coverage rate of one-step rolling forecasts form  
##		an uGARCHforecast class of object
##  Eg: fore = ugarchforecast(fit)
## 		use method: rate(fore)
##
##########################################################################
rate = function(object){
    UseMethod("rate", object)
}
setMethod("rate", 
    signature(object = "uGARCHforecast"),
    function(object){
        dist = object@model$modeldesc$distribution
        distpar = object@model$pars[c("lambda", "skew", "shape"),1];
        cri.u = qdist(dist, c(.95,.975), lambda = distpar["lambda"], skew = distpar["skew"], shape = distpar["shape"])
		cri.l = qdist(dist, c(.05,.025), lambda = distpar["lambda"], skew = distpar["skew"], shape = distpar["shape"])
        T0 = object@model$modeldata$T
        n.roll = object@forecast$n.roll
        ind.new = (T0+1):(T0+1+n.roll)
        y.new = object@model$modeldata$data[ind.new]
        y.hat = fitted(object)[1,]; 
        sig = sigma(object)[1,];
        out.90 = c(mean(y.new < (y.hat + cri.l[1]*sig)), mean(y.new > (y.hat + cri.u[1]*sig)))
        out.90 = c(1-sum(out.90), out.90)
        out.95 = c(mean(y.new < (y.hat + cri.l[2]*sig)), mean(y.new > (y.hat + cri.u[2]*sig)))
        out.95 = c(1-sum(out.95), out.95)
        out = rbind(out.95, out.90)
        dimnames(out)[[1]] = c("95% PI", "90% PI")
        dimnames(out)[[2]] = c("coverage", "below PI", "beyond PI")
        cat(paste("\n      One-Step Rolling Forecast   \n"))
        cat(paste("----------------------------------\n", sep = ""))
        print(round(out,4))
        invisible(out)
        }
)

##########################################################################
##
##	rate0()
## 	Calculate coverage rate of one-step rolling forecasts form  
##		an arfimaforecast class of object
##  Eg: fore = arfimaforecast(fit)
## 		use method: rate(fore)
##
##########################################################################
rate0 = function(object){
    UseMethod("rate0", object)
}
setMethod("rate0", 
    signature(object = "ARFIMAforecast"),
    function(object){
        dist = object@model$modeldesc$distribution
        distpar = object@model$pars[c("lambda", "skew", "shape"),1];
        cri.u = qdist(dist, c(.95,.975), lambda = distpar["lambda"], skew = distpar["skew"], shape = distpar["shape"])
		cri.l = qdist(dist, c(.05,.025), lambda = distpar["lambda"], skew = distpar["skew"], shape = distpar["shape"])
        T0 = object@model$modeldata$T
        n.roll = object@forecast$n.roll
        ind.new = (T0+1):(T0+1+n.roll)
        y.new = object@model$modeldata$data[ind.new]
        y.hat = fitted(object)[1,]; 
        sig = object@model$pars["sigma",1];
        out.90 = c(mean(y.new < (y.hat + cri.l[1]*sig)), mean(y.new > (y.hat + cri.u[1]*sig)))
        out.90 = c(1-sum(out.90), out.90)
        out.95 = c(mean(y.new < (y.hat + cri.l[2]*sig)), mean(y.new > (y.hat + cri.u[2]*sig)))
        out.95 = c(1-sum(out.95), out.95)
        out = rbind(out.95, out.90)
        dimnames(out)[[1]] = c("95% PI", "90% PI")
        dimnames(out)[[2]] = c("coverage", "below PI", "beyond PI")
        cat(paste("\n      One-Step Rolling Forecast   \n"))
        cat(paste("----------------------------------\n", sep = ""))
        print(round(out,4))
        invisible(out)
        }
)

#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
##   This file is based on R package rugarch by Alexios Ghalanos.
##
##   Define the function showShort for print shorter summary report
##   for ARFIMAfit object
##   
#################################################################################

showShort0 = function(object){
    UseMethod("showShort0", object)
}
setMethod("showShort0",
		signature(object = "ARFIMAfit"),
		function(object){
			model = object@model
			modelinc = model$modelinc
			cat(paste("\n*----------------------------------*", sep = ""))
			cat(paste("\n*          ARFIMA Model Fit        *", sep = ""))
			cat(paste("\n*----------------------------------*", sep = ""))
			cat("\nMean Model\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Distribution\t:", model$modeldesc$distribution,"\n")
			if(object@fit$convergence == 0){
				cat("\nOptimal Parameters")
				cat(paste("\n------------------------------------\n",sep=""))
				print(round(object@fit$matcoef,6), digits = 5)
				cat("\nRobust Standard Errors:\n")
				print(round(object@fit$robust.matcoef,6), digits = 5)
				if(!is.null(object@fit$hessian.message)){
					cat(paste("\n", object@fit$hessian.message))
				}
				cat("\nLogLikelihood :", object@fit$LLH, "\n")
				stdresid = object@fit$residuals/coef(object)["sigma"]
				itestm = t(infocriteria(object))
				cat("\nInformation Criteria")
				cat(paste("\n------------------------------------\n",sep=""))
				print(itestm,digits=5)
				cat("\nWeighted Ljung-Box Test on Standardized Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp1 = .weightedBoxTest(stdresid, p = 1, df = sum(modelinc[2:3]))
				print(tmp1, digits = 4)
				cat("\nH0 : No serial correlation\n")
				cat("\nWeighted Ljung-Box Test on Standardized Squared Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp2 = .weightedBoxTest(stdresid, p = 2, df = sum(modelinc[8:9]))
				print(tmp2, digits = 4)
				cat("\n\nARCH LM Tests")
				cat(paste("\n------------------------------------\n",sep=""))
				L2 = .archlmtest(stdresid, lags = 2)
				L5 = .archlmtest(stdresid, lags = 5)
				L10 = .archlmtest(stdresid, lags = 10)
				alm = matrix(0,ncol = 3,nrow = 3)
				alm[1,1:3] = c(L2$statistic, L2$parameter, L2$p.value)
				alm[2,1:3] = c(L5$statistic, L5$parameter, L5$p.value)
				alm[3,1:3] = c(L10$statistic, L10$parameter, L10$p.value)
				colnames(alm) = c("Statistic", "DoF", "P-Value")
				rownames(alm) = c("ARCH Lag[2]", "ARCH Lag[5]", "ARCH Lag[10]")
				print(alm,digits = 4)
			} else{
				cat("\nConvergence Problem:")
				cat("\nSolver Message:", object@fit$message,"\n\n")
				
			}
			return(invisible(object))
})

#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
##   This file is based on R package rugarch by Alexios Ghalanos.
##
##   Define the function showShort for print shorter summary report
##   for uGARCHfit object
##   
#################################################################################
#----------------------------------------------------------------------------------
# univariate show method / seperate for fit,sim and forecast
#----------------------------------------------------------------------------------
# fit show
showShort = function(object){
	UseMethod("showShort", object)
}
setMethod("showShort",
#setMethod("show",
		signature(object = "uGARCHfit"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			modelinc = object@model$modelinc
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*          GARCH Model Fit        *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\nGARCH Model\t: ", vmodel,"(", modelinc[8], ",", modelinc[9], ")\n", sep=""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model\t: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			cat("Mean Model\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Distribution\t:", model$modeldesc$distribution,"\n")
			if(object@fit$convergence == 0){
				cat("\nOptimal Parameters")
				cat(paste("\n------------------------------------\n",sep=""))
				print(round(object@fit$matcoef,6), digits = 5)
				if(!is.null(object@fit$hessian.message)){
					cat(paste("\n", object@fit$hessian.message))
				}
				cat("\nLogLikelihood :", object@fit$LLH, "\n")
				stdresid = object@fit$residuals/object@fit$sigma
				itestm = t(infocriteria(object))
				cat("\nInformation Criteria")
				cat(paste("\n------------------------------------\n",sep=""))
				print(itestm,digits=5)
				cat("\nWeighted Ljung-Box Test on Standardized Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp1 = .weightedBoxTest(stdresid, p = 1, df = sum(modelinc[2:3]))
				print(tmp1, digits = 4)
				cat(paste("d.o.f=", sum(modelinc[2:3]), sep=""))
				cat("\nH0 : No serial correlation\n")
				cat("\nWeighted Ljung-Box Test on Standardized Squared Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp2 = .weightedBoxTest(stdresid, p = 2, df = sum(modelinc[8:9]))
				print(tmp2, digits = 4)
				cat(paste("d.o.f=", sum(modelinc[8:9]), sep=""))
				cat("\n\nWeighted ARCH LM Tests")
				cat(paste("\n------------------------------------\n",sep=""))
				gdf = sum(modelinc[8:9])
				L2 = .weightedarchlmtest(residuals(object), sigma(object), lags  = gdf+1, fitdf=gdf)
				L5 = .weightedarchlmtest(residuals(object), sigma(object), lags  = gdf+3, fitdf=gdf)
				L10 = .weightedarchlmtest(residuals(object), sigma(object), lags = gdf+5, fitdf=gdf)
				alm = matrix(0, ncol = 4, nrow = 3)
				alm[1,1:4] = as.numeric(c(L2$statistic, L2$parameter, L2$p.value))
				alm[2,1:4] = as.numeric(c(L5$statistic, L5$parameter, L5$p.value))
				alm[3,1:4] = as.numeric(c(L10$statistic, L10$parameter, L10$p.value))
				colnames(alm) = c("Statistic", "Shape", "Scale", "P-Value")
				rownames(alm) = c(paste("ARCH Lag[",gdf+1,"]",sep=""), paste("ARCH Lag[",gdf+3,"]",sep=""), paste("ARCH Lag[",gdf+5,"]",sep=""))
				print(alm,digits = 4)
				cat("\nElapsed time :", object@fit$timer,"\n\n")
			} else{
				cat("\nConvergence Problem:")
				cat("\nSolver Message:", object@fit$message,"\n\n")

			}
			invisible(object)
		})
        
#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
##   This file is part of the R package rugarch.
##
##   The R package rugarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rugarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

# Q-Statistics on Standardized Residuals
.box.test = function(stdresid, p=1, df = 0)
{
	if(any(!is.finite(stdresid))) stdresid[!is.finite(stdresid)]=0
	# p=1 normal case, p=2 squared std. residuals
	# Q-Statistics on Standardized Residuals
	#H0 : No serial correlation ==> Accept H0 when prob. is High [Q < Chisq(lag)]
	box10 = Box.test(stdresid^p, lag = 1, type = "Ljung-Box", fitdf = 0)
	box15 = Box.test(stdresid^p, lag = df+1, type = "Ljung-Box", fitdf = df)
	box20 = Box.test(stdresid^p, lag = df+5, type = "Ljung-Box", fitdf = df)
	LBSR<-matrix(NA,ncol=2,nrow=3)
	LBSR[1:3,1] = c(box10$statistic[[1]],box15$statistic[[1]],box20$statistic[[1]])
	LBSR[1:3,2] = c(box10$p.value[[1]],box15$p.value[[1]],box20$p.value[[1]])
	rownames(LBSR) = c(paste("Lag[1]",sep=""), paste("Lag[p+q+1][",df+1,"]",sep=""), paste("Lag[p+q+5][",df+5,"]",sep=""))
	colnames(LBSR) = c("statistic","p-value")
	return(LBSR)
}


.weightedBoxTest = function(stdresid, p=1, df = 0)
{
	if(any(!is.finite(stdresid))) stdresid[!is.finite(stdresid)]=0
	# p=1 normal case, p=2 squared std. residuals
	# Q-Statistics on Standardized Residuals
	#H0 : No serial correlation ==> Accept H0 when prob. is High [Q < Chisq(lag)]
	box10 = Weighted.Box.test(stdresid, lag = 1, type = "Ljung-Box", fitdf = 0, if(p==2) sqrd.res = TRUE else sqrd.res = FALSE)
	box15 = Weighted.Box.test(stdresid, lag = max(2, 2*df+df-1), type = "Ljung-Box", fitdf = df, if(p==2) sqrd.res = TRUE else sqrd.res = FALSE)
	box20 = Weighted.Box.test(stdresid, lag = max(5, 4*df+df-1), type = "Ljung-Box", fitdf = df, if(p==2) sqrd.res = TRUE else sqrd.res = FALSE)
	LBSR<-matrix(NA,ncol=2,nrow=3)
	LBSR[1:3,1] = c(box10$statistic[[1]],box15$statistic[[1]],box20$statistic[[1]])
	LBSR[1:3,2] = c(box10$p.value[[1]],box15$p.value[[1]],box20$p.value[[1]])
	rownames(LBSR) = c(paste("Lag[1]",sep=""), paste("Lag[2*(p+q)+(p+q)-1][",max(2, 2*df+df-1),"]",sep=""), paste("Lag[4*(p+q)+(p+q)-1][",max(5, 4*df+df-1),"]",sep=""))
	colnames(LBSR) = c("statistic","p-value")
	return(LBSR)
}


.archlmtest = function (x, lags, demean = FALSE)
{
	if(any(!is.finite(x))) x[!is.finite(x)] = 0
	x = as.vector(x)
	if(demean) x = scale(x, center = TRUE, scale = FALSE)
	lags = lags + 1
	mat = embed(x^2, lags)
	arch.lm = summary(lm(mat[, 1] ~ mat[, -1]))
	STATISTIC = arch.lm$r.squared * length(resid(arch.lm))
	names(STATISTIC) = "Chi-squared"
	PARAMETER = lags - 1
	names(PARAMETER) = "df"
	PVAL = 1 - pchisq(STATISTIC, df = PARAMETER)
	METHOD = "ARCH LM-test"
	result = list(statistic = STATISTIC, parameter = PARAMETER,
			p.value = PVAL, method = METHOD)
	class(result) = "htest"
	return(result)
}


.weightedarchlmtest = function (x, sigma, lags, fitdf = 2, demean = FALSE)
{
	if(any(!is.finite(x))) x[!is.finite(x)] = 0
	x = as.vector(x)
	if(demean) x = scale(x, center = TRUE, scale = FALSE)
	result = Weighted.LM.test(x, sigma^2, lag = lags, type = c("correlation", "partial")[1], fitdf = fitdf, weighted=TRUE) 
	return(result)
}

######
## Weighted.Box.test -- The Weighted Box-Pierce, Ljung-Box, Monti, McLeod-Li type test for fitted ARMA and detection of nonlinear Models
##
## Quick intro for the input parameters
##
##      x       -- the residuals or initial unfitted time series
##     lag      -- The lag we wish to test at, this is m or M in the literature, default=1
##     type     -- The type of test, default="Box-Pierce"
##    fitdf     -- The number of ARMA parameter that have been fit to the series x, default=0
##
## The following are boolean flags to perform transformations to detect nonlinear processes, all default=FALSE
##
##    sqrt.res  -- Should the residuals be squared to detect for nonlinear models?  If so, fitdf is ignore, the type "Ljung-Box" is really the "McLeod-Li" type test
## log.sqrd.res -- Take the log of the squared residuals, generally simulations show this is less powerful than the sqrd.res or abs.res
##    abs.res   -- Take the absolute value of the residuals, similar to the above two boolean flags.  Simulations indicate abs.res is more powerful at detecting long memory processes
##
## For backward compatability
##
##   weighted   -- For backward compatability, If set to FALSE, you will perform the traditional Box-Pierce, Ljung-Box or Monti or McLeod Li type if using a transformation, default=TRUE.
######

Weighted.Box.test <- function (x, lag = 1, type = c("Box-Pierce", "Ljung-Box", "Monti"), fitdf = 0, sqrd.res = FALSE, log.sqrd.res = FALSE, abs.res = FALSE, weighted=TRUE)
{
 ### Error Checking
 ###
    if (NCOL(x) > 1) 
       stop("x is not a vector or univariate time series");
    if (lag < 1)
       stop("Lag must be positive");
    if (fitdf < 0)
       stop("Fitdf cannot be negative");
    if (fitdf >= lag)
       stop("Lag must exceed fitted degrees of freedom");
    if( (sqrd.res && log.sqrd.res) || (sqrd.res && abs.res) || (log.sqrd.res && abs.res) )
       stop("Only one option of: sqrd.res, log.sqrd.res or abs.res can be selected");

 ### Find the type
    DNAME <- deparse(substitute(x))
    type <- match.arg(type)

 ### Transform if checking for nonlinear
    if(abs.res) {
       x <- abs(x);
    }
    if(sqrd.res || log.sqrd.res) {
       x <- x^2;
    }
    if(log.sqrd.res) {
       x <- log(x);
    }

 if(weighted) 
 {
    if( type == "Monti") {
       METHOD <- "Weighted Monti test (Gamma Approximation)"
       cor <- acf(x, lag.max = lag, type="partial", plot=FALSE, na.action=na.pass)
       obs <- cor$acf[1:lag];
    }
    else {
       cor <- acf(x, lag.max = lag, type="correlation", plot=FALSE, na.action=na.pass)
       obs <- cor$acf[2:(lag + 1)];
    }

    if(type == "Ljung-Box") {
       METHOD <- "Weighted Ljung-Box test (Gamma Approximation)"
    }

    n <- sum(!is.na(x))

    weights <- (lag - 1:lag+1)/(lag);

    if (type == "Box-Pierce") {
       METHOD <- "Weighted Box-Pierce test (Gamma Approximation)"
       STATISTIC <- n * sum(weights*obs^2)
    }
    else {
       STATISTIC <- n * (n + 2) * sum(weights*(1/seq.int(n - 1, n - lag)*obs^2))
    }

    if(sqrd.res) {
       fitdf <- 0;
       names(STATISTIC) <- "Weighted X-squared on Squared Residuals for detecting nonlinear processes"
    }
    else if(log.sqrd.res) { 
       fitdf <- 0;
       names(STATISTIC) <- "Weighted X-squared on Log-Squared Residuals for detecting nonlinear processes";
    }
    else if(abs.res) { 
       fitdf <- 0;
       names(STATISTIC) <- "Weighted X-squared on Absolute valued Residuals for detecting nonlinear processes";
    }
    else {
       names(STATISTIC) <- "Weighted X-squared on Residuals for fitted ARMA process"
    }

    shape <- (3/4)*(lag+1)^2*lag/(2*lag^2 + 3*lag + 1 - 6*lag*fitdf);
    scale <- (2/3)*(2*lag^2 + 3*lag + 1 - 6*lag*fitdf)/lag/(lag+1);

    PARAMETER <- c(shape, scale);
    names(PARAMETER) <- c("Shape", "Scale")

    PVAL <- 1 - pgamma(STATISTIC, shape=shape, scale=scale)
    names(PVAL) <- "Approximate p-value"
 }
 else #Not weighted
 {
    if( type == "Monti") {
       METHOD <- "Monti test"
       cor <- acf(x, lag.max = lag, type="partial", plot=FALSE, na.action=na.pass)
       obs <- cor$acf[1:lag];
    }
    else {
       cor <- acf(x, lag.max = lag, type="correlation", plot=FALSE, na.action=na.pass)
       obs <- cor$acf[2:(lag + 1)];
    }

    if(type == "Ljung-Box") {
       METHOD <- "Ljung-Box test"
    }

    n <- sum(!is.na(x))


    if (type == "Box-Pierce") {
       METHOD <- "Box-Pierce test"
       STATISTIC <- n * sum(obs^2)
    }
    else {
       STATISTIC <- n * (n + 2) * sum((1/seq.int(n - 1, n - lag)*obs^2))
    }

    if(sqrd.res) {
       fitdf <- 0;
       names(STATISTIC) <- "X-squared on Squared Residuals for detecting nonlinear processes"
    }
    else if(log.sqrd.res) { 
       fitdf <- 0;
       names(STATISTIC) <- "X-squared on Log-Squared Residuals for detecting nonlinear processes";
    }
    else if(abs.res) { 
       fitdf <- 0;
       names(STATISTIC) <- "X-squared on Absolute valued Residuals for detecting nonlinear processes";
    }
    else {
       names(STATISTIC) <- "X-squared on Residuals for fitted ARMA process"
    }

    mydf <- lag - fitdf;

    PARAMETER <- c(mydf);
    names(PARAMETER) <- c("df")

    PVAL <- 1 - pchisq(STATISTIC, df=mydf)
    names(PVAL) <- "p-value"

 }
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
       p.value = PVAL, method = METHOD, data.name = DNAME), 
       class = "htest")
}

######
## Weighted.LM.test -- The Weighted Li-Mak type test
##
## Quick intro for the input parameters
##
##      x       -- the residuals or initial unfitted time series
##     h.t      -- the sample conditional variances
##     lag      -- The lag we wish to test at, this is m or M in the literature, default=1
##     type     -- The type of test, here just based on acfs or partial-acfs, default="correlation"
##    fitdf     -- The number of ARCH parameter that have been fit to the series x, default=1
##                 **Note** If fitdf=0, h.t is empty, you want to be performing one of the test in Weighted.Box.test() 
##
## For backward compatability
##
##   weighted   -- For backward compatability, If set to FALSE, you will perform the Li-Mak test, default=TRUE
######

Weighted.LM.test <- function (x, h.t, lag = 1, type = c("correlation", "partial"), fitdf = 1, weighted=TRUE) 
{
 ### Error Checking
 ###
    if (NCOL(x) > 1) 
       stop("x is not a vector or univariate time series");
    if (fitdf >= lag)
       stop("Lag must exceed fitted degrees of freedom");
    if (fitdf < 1)
       stop("Fitted degrees of freedom must be positive");
    if( !(length(x)==length(h.t)) )
       stop("Length of x and h.t must match");

    DNAME <- deparse(substitute(x))
    type <- match.arg(type)

    x <- x^2/h.t;

    if( type == "partial") {
       cor <- acf(x, lag.max = lag, type="partial", plot=FALSE, na.action=na.pass)
       obs <- cor$acf[1:lag];
    }
    else {
       cor <- acf(x, lag.max = lag, type="correlation", plot=FALSE, na.action=na.pass)
       obs <- cor$acf[2:(lag + 1)];
    }


    if(type == "correlation" && weighted) {
       METHOD <- "Weighted Li-Mak test on autocorrelations (Gamma Approximation)"
    }
    else if(type == "partial" && weighted) {
       METHOD <- "Weighted Li-Mak test on partial autocorrelations (Gamma Approximation)"
    }
    else if(type == "correlation" && !weighted) {
       METHOD <- "Li-Mak test on autocorrelations (Chi-Squared Approximation)"
    }
    else {
       METHOD <- "Li-Mak test on partial autocorrelations (Chi-Squared Approximation)"
    }

    n <- sum(!is.na(x))
    if(weighted) {
       weights <- (lag - (fitdf+1):lag + (fitdf+1) )/lag;
       obs <- obs[(fitdf+1):lag];
       STATISTIC <- n * sum(weights*obs^2);
       names(STATISTIC) <- "Weighted X-squared on Squared Residuals for fitted ARCH process";
       shape <- (3/4)*(lag + fitdf + 1)^2*(lag - fitdf)/(2*lag^2 + 3*lag + 2*lag*fitdf + 2*fitdf^2 + 3*fitdf + 1);
       scale <- (2/3)*(2*lag^2 + 3*lag + 2*lag*fitdf + 2*fitdf^2 + 3*fitdf + 1)/(lag*(lag + fitdf + 1));
       PARAMETER <- c(shape, scale);
       names(PARAMETER) <- c("Shape", "Scale")
    }
    else {
       weights <- rep(1,(lag-fitdf) );
       obs <- obs[(fitdf+1):lag];
       STATISTIC <- n * sum(weights*obs^2);
       names(STATISTIC) <- "X-squared on Squared Residuals for fitted ARCH process"
       shape <- (lag-fitdf)/2;          # Chi-squared df in Gamma form.
       scale <- 2;
       PARAMETER <- c((lag-fitdf));
       names(PARAMETER) <- c("Degrees of Freedom");
    }

    PVAL <- 1 - pgamma(STATISTIC, shape=shape, scale=scale)
    names(PVAL) <- "Approximate p-value"

    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
       p.value = PVAL, method = METHOD, data.name = DNAME), 
       class = "htest")
}
