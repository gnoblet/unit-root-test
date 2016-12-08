###            TERM PAPER          ###
### Econometrics 3B : Times Series ###
###         Guillaume Noblet       ###

rm(list=ls())


# 1 Packages, functions, and others  ------------------------------------------------------------



# 1.1 Packages ------------------------------------------------------------


library("forecast")      # Used for Arima and auto.arima
library("tseries")       
library("lmtest") 
library("dynlm")
library("vars")
library("ArgumentCheck") # Necessary to switch to "checkmate" early 2017
library("seasonal")      # May we use it? 


# 1.2 Functions -----------------------------------------------------------




###            AUTO.ADF             ###
### Automatic function for ADF-test ###
### To be used only with ts objects ###

require(urca)
require(ArgumentCheck)

# Arguments to be set:
# data: vector to be tested for unit root;
# lags: lags selection to be achieved ("AIC": Akaike, "BIC": Bayes, 
#   "Fixed": default lags);
# maxlag: number of lags to be included.
# alpha: level of significance

auto.adf <- function(data, lags="AIC", maxlag=36, alpha=0.05){
  
  # Checking arguments
  #* Establish a new 'ArgCheck' object
  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if data is not a ts object
  if (is.ts(data)==FALSE)
    ArgumentCheck::addError(
      msg = paste0("The object to be tested for a unit root is not a ts object.", 
                   " auto.adf is designed to run only objects of class 'ts'.") ,
      argcheck = Check
    )
  #* Add an error if maxlag > 36
  if (maxlag > 36) 
    ArgumentCheck::addError(
      msg = "'maxlag' must be lower than 36, which is accurate for monthly data",
      argcheck = Check
    )
  #* Add an error if alpha is not 0.01, 0.05 or 0.1
  if (alpha != 0.01 & alpha != 0.05 & alpha != 0.1)
    ArgumentCheck::addError(
      msg = "'alpha' is the signficance level and  must be either 0.01, 0.05, or 0.1",
      argcheck = Check
    )
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)
  
  # Pre-allocation of space
  name <- deparse(substitute(data))      # get the name of the data set 
  type <- c("trend", "drift", "none")    # defining all types
  # Setting critical value corresponding to significance level
  # cv is then the corresponding column in the distribution table
  if (alpha == 0.01){
    cv <- 2 
  }
  else if (alpha == 0.1){
    cv <- 4
  }
  else if (alpha == 0.05){
    cv <- 3
  }
  
  # Loop over the different types of tests
  results <- list(trend=NA, drift=NA, none=NA)  
  mod <- list(trend=NA, drift=NA, none=NA)
  goto <- 0
  for (i in 1:3){
    model <- ur.df(data, type=type[i], lags=maxlag, selectlags=lags)
    results[[i]] <- cbind(t(model@teststat), model@cval)
    mod[[i]] <- summary(model)
  }
  
  # Decision tree for unit root
  #* Test for a unit root with a trend and drift model
  if (results$trend[1,1] < results$trend[1,cv]){
    ur <- paste0("tau3 is rejected at a ", alpha, " significance level. The ts object ", name, " has no unit root.", 
                 " The chosen model is with trend and drift.")
  }
  else if (results$trend[3,1] > results$trend[3,cv]){
    if (mod$trend@testreg$coefficients[2,4] < 0 & mod$trend@testreg$coefficients[2,4] < qt(0.05, mod$trend@testreg$df[2])){
      ur <- paste0("tau 3 is rejected at a ", alpha, "significance level using a t-test. The ts object ", name, " has no unit root.")
    }
    else {
      goto <- 2
    }
  }
  else {
    goto <- 2
  }
  #* Test for a unit root with a drift only model
  if (goto == 2) {
    if (results$drift[1,1] < results$drift[1,cv])
      ur <- paste0("tau2  is rejected at a ", alpha, " significance level. The ts object ", name, " has no unit root.", 
                   " The chosen model is with drift only")
    else if (results$drift[2,1] > results$drift[2,cv]){
      if (mod$drift@testreg$coefficients[2,4] < 0 & mod$drift@testreg$coefficients[2,4] <  qt(0.05, mod$drift@testreg$df[2])){
        ur <- paste0("tau2 is rejected at a ", alpha, " significance level performing a t-test. The ts object ", name, " has no unit root.")
      }
      else {
        goto <- 3
      }
    }
    else {
      goto <- 3
    }
  }
  #* Test for a unit root with a drift only model
  if (goto == 3) {
    if (results$none[1,1] < results$none[1, cv]) {
      ur <- paste0("tau1  is rejected at a ", alpha, " significance level. The ts object ", name, " has no unit root.", 
                   " The chosen model is with neither trend nor drift.")
    }
    else {
      ur <- paste0("The ts object ", name, " has a unit root; take its first difference and re-run the test.") 
    }
  }
  return(cat(ur))
}




###             AUTO.PP             ###
### Automatic function for pp.test  ###
### To be used only with ts objects ###

require(urca)
require(ArgumentCheck)

# Arguments to be set:
# data: vector to be tested for unit root;
# alpha: level of significance (default: 0.05)

# Phillips-Perron test:
# Null hypothesis: a time series is integrated of order 1 (unit root)
# The null hypothesis is rejected at a alpha significance level if the p-value 
# is smaller than the critical value

auto.pp <- function(data, alpha = 0.05){
  
  # Checking arguments
  #* Establish a new 'ArgCheck' object
  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if data is not a ts object
  if (is.ts(data)==FALSE)
    ArgumentCheck::addError(
      msg = paste0("The object to be tested for a unit root is not a ts object.", 
                   " auto.pp.test is designed to run only objects of class 'ts'.") ,
      argcheck = Check
    )
  #* Add an error if alpha is not 0.01, 0.05 or 0.1
  if (alpha != 0.01 & alpha != 0.05 & alpha != 0.1)
    ArgumentCheck::addError(
      msg = "'alpha' is the signficance level and  must be either 0.01, 0.05, or 0.1",
      argcheck = Check
    )
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)
  
  # Pre-allocation of space
  name <- deparse(substitute(data))      # get the name of data
  type = c("trend", "constant")
  # Setting critical value corresponding to significance level
  # cv is then the corresponding column in the distribution table
  if (alpha == 0.01){
    cv <- 1 
  }
  else if (alpha == 0.1){
    cv <- 3
  }
  else if (alpha == 0.05){
    cv <- 2
  }
  
  
  # Loop over the different types of tests
  results <- list(trend=NA, drift=NA) 
  for (i in 1:2){
    model <- ur.pp(data, model = type[i], type = "Z-tau")
    results[[i]] <- c(model@teststat, model@cval[cv])
  }

  # Test
  # Decision algorithm for unit root
  # Model with trend
  if (results$trend[1] < results$trend[2]){
    ur <- paste0("Phillips-Perron test: ", name, " has no unit root at a ", alpha,
                 " significance level for a with trend and drift model.")
  }
  else if (results$drift[1] < results$drift[2]){
    ur <- paste0("Phillips-Perron test: ", name, " has no unit root at a ", alpha,
                 " significance level for a with drit only model.")
  }
  else {
    ur <- paste0("Phillips-Perron test: the ts object ", name, " has a unit root; take its first difference and re-run the test.")
  }
  
  return(cat(ur))
}




###           KPSS.LINEAR.TEST            ###
###   Automatic function for kpss.test    ###
###   To be used only with ts objects     ###

require(urca)
require(ArgumentCheck)

# Arguments to be set:
# data: vector to be tested for unit root;
# alpha: level of significance (default: 0.05)

# KPSS test:
# Null hypothesis: the time series is weakly stationary
# Alternative: the time series is non-stationary
# The null is rejected if the p-value is above a critical size 

auto.kpss <- function(data, alpha = 0.05){
  
  # Checking arguments
  #* Establish a new 'ArgCheck' object
  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if data is not a ts object
  if (is.ts(data)==FALSE)
    ArgumentCheck::addError(
      msg = paste0("The object to be tested for a unit root is not a ts object.", 
                   " auto.pp.test is designed to run only objects of class 'ts'.") ,
      argcheck = Check
    )
  #* Add an error if alpha is not 0.01, 0.05 or 0.1
  if (alpha != 0.01 & alpha != 0.05 & alpha != 0.1)
    ArgumentCheck::addError(
      msg = "'alpha' is the signficance level and  must be either 0.01, 0.05, or 0.1",
      argcheck = Check
    )
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)
  
  # Pre-allocation of space
  name <- deparse(substitute(data))  # get the name of data
  type = c("tau", "mu")              # tau: trend and drif; mu: drift only
  # Setting critical value corresponding to significance level
  # cv is then the corresponding column in the distribution table
  if (alpha == 0.01){
    cv <- 1 
  }
  else if (alpha == 0.1){
    cv <- 3
  }
  else if (alpha == 0.05){
    cv <- 2
  }
  
  # Loop over the different types of tests
  results <- list(trend=NA, drift=NA) 
  for (i in 1:2){
    model <- ur.kpss(data, type = type[i])
    results[[i]] <- c(model@teststat, model@cval[cv])
  }
  
  # Test
  # Decision algorithm for unit root
  # Model with trend
  if (results$trend[1] < results$trend[2]){
    ur <- paste0("KPSS test: ", name, " has no unit root at a ", alpha,
                 " significance level for a with trend and drift model.")
  }
  else if (results$drift[1] < results$drift[2]){
    ur <- paste0("KPSS test: ", name, " has no unit root at a ", alpha,
                 " significance level for a with drit only model.")
  }
  else {
    ur <- paste0("KPSS test: the ts object ", name, " has a unit root; take its first difference and re-run the test.")
  }
  
  return(cat(ur))
}




###                    IRF.PLOT.SVAREST             ###
###  Plot impulse response functions and their CIs  ###
###       To be used only with svarest objects      ###

require(vars)

# This function plot ALL the shocks regarding the SVAR model. It cannot be used if one
# is just interested in some shocks.

# Arguments to be used:
# model: svarest object;
# n.ahead: time ahead after the shocks;

# Default values are:
# – bootstrap = TRUE, bootstraping method to get the confidence intervals;
# – cumul = TRUE, cumulative impulse response function;
# – ci = 0.95, confidence interval is 0.95;
# – zero = FALSE; if TRUE, force abline(0) to be displayed;
# – linte.type = 3, dotted line type (see "lty" in par{graphics} for more info)

irf.plot.svarest <- function(model, n.ahead, bootstrap = TRUE, cumul = TRUE, ci = 0.95, zero = FALSE, line.type = 3){
  # Computing the impulse response function
  IRF <- irf(model, n.ahead = n.ahead, boot = bootstrap, cumulative = cumul, ci = ci)
  # Preparing the plots
  names <- colnames(model$var$y)
  M <- length(names)
  # Plotting the IRFs
  par(mfrow=c(M,M),oma = c(2,2,1,1), mar =c(2,2,2,2), mgp = c(2,1,0), xpd = FALSE)
  for (i in 1:M){
    for (j in 1:M){
      if (zero == FALSE) {
        ymax <- max(IRF$irf[[i]][,j], IRF$Lower[[i]][,j], IRF$Upper[[i]][,j])
        ymin <- min(IRF$irf[[i]][,j], IRF$Lower[[i]][,j], IRF$Upper[[i]][,j])
        plot(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$irf[[i]][,j], type = "n", ylim = c(ymin, ymax), xlab = "",
             ylab= "", main = paste0("Shock from ", names[i], " to ", names[j]))
        polygon(c(x = seq(along.with = IRF$irf[[i]][,j]), rev(seq(along.with = IRF$irf[[i]][,j]))), c(IRF$Lower[[i]][,j], rev(IRF$Upper[[i]][,j])),
                col = "gray90", border = NA)
        lines(x = c(0,seq(along.with = IRF$irf[[i]][,j])), y = rep(0, 1 + length(IRF$irf[[i]][,j])), col = "black")
        lines(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$irf[[i]][,j], col = "blue")
        lines(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$Lower[[i]][,j], col = "red", lty = line.type)
        lines(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$Upper[[i]][,j], col = "red", lty = line.type)
      }
      else {
        ymax <- max(IRF$irf[[i]][,j], IRF$Lower[[i]][,j], IRF$Upper[[i]][,j], 0)
        ymin <- min(IRF$irf[[i]][,j], IRF$Lower[[i]][,j], IRF$Upper[[i]][,j], 0)
        plot(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$irf[[i]][,j], type = "n", ylim = c(ymin, ymax), xlab = "",
             ylab= "", main = paste0("Shock from ", names[i], " to ", names[j]))
        polygon(c(seq(along.with = IRF$irf[[i]][,j]), rev(seq(along.with = IRF$irf[[i]][,j]))), c(IRF$Lower[[i]][,j], rev(IRF$Upper[[i]][,j])),
                col = "gray90", border = NA)
        lines(x = c(0,seq(along.with = IRF$irf[[i]][,j])), y = rep(0, 1 + length(IRF$irf[[i]][,j])), col = "black")
        lines(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$irf[[i]][,j], col = "blue")
        lines(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$Lower[[i]][,j], col = "red", lty = line.type)
        lines(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$Upper[[i]][,j], col = "red", lty = line.type)
      }
    }
  }
  par(mfrow=c(1,1))
}




###                   IRF.PLOT.VAREST              ###
### Plot impulse response functions and their CIs  ###
###       To be used only with varest objects      ###

require(vars)

# This function plot ALL the shocks regarding the SVAR model. It cannot be used if one
# is just interested in some shocks.

# Arguments to be used:
# model: varest object;
# n.ahead: time ahead after the shocks;

# Default values are:
# – bootstrap = TRUE, bootstraping method to get the confidence intervals;
# – cumul = TRUE, cumulative impulse response function;
# – ci = 0.95, confidence interval is 0.95;
# – zero = FALSE; if TRUE, force abline(0) to be displayed;
# – line.type = 3, dotted line type (see "lty" in par{graphics} for more info)

irf.plot.varest <- function(model, n.ahead, bootstrap = TRUE, cumul = TRUE, ci = 0.95, zero = FALSE, line.type = 3){
  # Computing the impulse response function
  IRF <- irf(model, n.ahead = n.ahead, boot = bootstrap, cumulative = cumul, ci = ci)
  # Preparing the plots
  names <- colnames(model$y)
  M <- length(names)
  # Plotting the IRFs
  par(mfrow=c(M,M),oma = c(2,2,1,1), mar =c(2,2,2,2), mgp = c(2,1,0), xpd = FALSE)
  for (i in 1:M){
    for (j in 1:M){
      if (zero == FALSE) {
        ymax <- max(IRF$irf[[i]][,j], IRF$Lower[[i]][,j], IRF$Upper[[i]][,j])
        ymin <- min(IRF$irf[[i]][,j], IRF$Lower[[i]][,j], IRF$Upper[[i]][,j])
        plot(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$irf[[i]][,j], type = "n", ylim = c(ymin, ymax), xlab = "",
             ylab= "", main = paste0("Shock from ", names[i], " to ", names[j]))
        polygon(c(x = seq(along.with = IRF$irf[[i]][,j]), rev(seq(along.with = IRF$irf[[i]][,j]))), c(IRF$Lower[[i]][,j], rev(IRF$Upper[[i]][,j])),
                col = "gray90", border = NA)
        lines(x = c(0,seq(along.with = IRF$irf[[i]][,j])), y = rep(0, 1 + length(IRF$irf[[i]][,j])), col = "black")
        lines(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$irf[[i]][,j], col = "blue")
        lines(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$Lower[[i]][,j], col = "red", lty = line.type)
        lines(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$Upper[[i]][,j], col = "red", lty = line.type)
      }
      else {
        ymax <- max(IRF$irf[[i]][,j], IRF$Lower[[i]][,j], IRF$Upper[[i]][,j], 0)
        ymin <- min(IRF$irf[[i]][,j], IRF$Lower[[i]][,j], IRF$Upper[[i]][,j], 0)
        plot(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$irf[[i]][,j], type = "n", ylim = c(ymin, ymax), xlab = "",
             ylab= "", main = paste0("Shock from ", names[i], " to ", names[j]))
        polygon(c(seq(along.with = IRF$irf[[i]][,j]), rev(seq(along.with = IRF$irf[[i]][,j]))), c(IRF$Lower[[i]][,j], rev(IRF$Upper[[i]][,j])),
                col = "gray90", border = NA)
        lines(x = c(0,seq(along.with = IRF$irf[[i]][,j])), y = rep(0, 1 + length(IRF$irf[[i]][,j])), col = "black")
        lines(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$irf[[i]][,j], col = "blue")
        lines(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$Lower[[i]][,j], col = "red", lty = line.type)
        lines(x = seq(along.with = IRF$irf[[i]][,j]), y = IRF$Upper[[i]][,j], col = "red", lty = line.type)
      }
    }
  }
  par(mfrow=c(1,1))
}






###                   GROWTH.RATE                   ###
###    Calculate the growth rates of ts objects     ###
###   If set up, calculate an aggregation and then  ###
###                 the growth rates                ###   

require(stats)
require(ArgumentCheck)

# Default values:
# n: 1, n should be the "speaking" number of lags
# percent.change = FALSE
# aggregate = FALSE
# freq = frequency of the ts object
# bymean = FALSE

growthrate <- function(ts, n = 1, percent.change = FALSE, agg = FALSE, freq = attributes(ts)$tsp[3], bymean = FALSE){
  # Checking arguments
  #* Establish a new 'ArgCheck' object
  Check <- ArgumentCheck::newArgCheck()
  
  if (is.logical(agg) == FALSE)
    ArgumentCheck::addError(
      msg = "agg must be a boolean",
      argcheck = Check
    )

  if (is.logical(percent.change) == FALSE) 
    ArgumentCheck::addError(
      msg = "percent.change muse be a boolean",
      argcheck = Check
    )
  
  if (is.logical(bymean) == FALSE)
    ArgumentCheck::addError(
      msg = "bymean must be a boolean",
      argcheck = Check
    )  
  
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)
  
  # Computing the growth rate with respect to n lags
  if (agg == FALSE) {
    if (percent.change == FALSE){
      gr <- ts / lag(ts, k = - n ) - 1
    }
    else {
      gr <- (ts / lag(ts, k = - n ) - 1) * 100
    }    
  }
  else if (agg == TRUE) {
    if (percent.change == FALSE){
      if (bymean == FALSE) {
        gr <- aggregate(ts, nfrequency = freq)
        gr <- gr / lag(gr, k = - n ) - 1
      }
      else {
        gr <- aggregate(ts, nfrequency = freq, mean)
        gr <- gr / lag(gr, k = - n ) - 1
      }
    }
    else {
      if (bymean == FALSE) {
        gr <- aggregate(ts, nfrequency = freq)
        gr <- (gr / lag(gr, k = - n ) - 1) * 100
      }
      else {
        gr <- aggregate(ts, nfrequency = freq, mean)
        gr <- (gr / lag(gr, k = - n ) - 1) * 100
     }
    }
  }
  return(gr)
}





# 2 Data -----------------------------------------------------------------




# 2.1 Preparation of the ts objects ---------------------------------------




# 2.1.1 Crude Oil Domestic First Purchase Price ---------------------------



## This data set is also from the EIA
## Data set with different variables regarding oil price

## Read the .csv file and import it as a data.frame
cop <- read.csv("EIA_cop.csv")

## Keep only date, value and the Column_Order
cop <- subset(cop, select = c(YYYYMM, Value, Column_Order))

## Keep only crude oil domestic first purchase price
cop <- subset(cop, Column_Order == 1)

## Remove all data whose date ends with 13 (does not include 2013), which are yearly data before 1973
cop <- subset(cop, !grepl("^.+13$", YYYYMM))

## Window from 1973-01 to 2015-04
cop <- subset(cop, YYYYMM < 201601 & YYYYMM > 197312)

## Check for unknown values
table(cop$Value == "Not Available")
# 12 values are unknown, there are the ones from year 1973
# Here we can't assume that the price was 0
# Hence, set the start year for all variables to 1974-01

## Set new starting date
cop <- subset(cop, YYYYMM > 197312)
start.year <- 1974
start.month <- 1

## Construct the ts object
cop <- ts(as.numeric(as.character(cop$Value)), 
          start = c(start.year, start.month), frequency = 12)



# 2.1.2 Real gross domestic product ---------------------------------------



## The real GDP data are quarterly seasonally-adjusted FRED data using the 2009 
## chained-method in order to control for inflation
## Unit: billion of dollars
## From 1947-01 to 2016-01

## Read the .csv file and import it as a data.frame
gdp <- read.csv2("FRED_gdp_US.csv")

## Create the ts-object
gdp <- ts(gdp$VALUE, start = c(1947,1), frequency = 4)

## Window gdp to the same starting and ending date as the other variables
gdp <- window(gdp, start = c(start.year,start.month))
gdp <- window(gdp, end = c(2015,4))



# 2.1.3 Energy  consumption ---------------------------------------------------



## Import the energy consumption's data set from the EIA (Energy information administration)
## The data are monthly data from January 1973 to January 2016 for all variables
## The data.frame object contains the following columns:
##  – the data sets' primary key: MSN;
##  – the time of the observation (mostly monthly data): YYYYMM;
##  – the value of the observation: Value;
##  – the order of the column which relates to colomn 5: Column_Order;
##  – the unit of the observations: Unit in Quadrillion Btu
## More info about column order:
##  – 1  Coal consumption
##  – 2  Natural Gas Consumption (Excluding Supplemental Gaseous Fuels)
##  – 3  Petroleum Consumption (Excluding Biofuels)
##  – 4  Total Fossil Fuels Consumption
##  – 5  Nuclear Electric Power Consumption
##  – 6  Hydroelectric Power Consumption
##  – 7  Geothermal Energy Consumption
##  – 8  Solar/PV Energy Consumption
##  – 9  Wind Energy Consumption
##  – 10 Biomass Energy Consumption
##  – 11 Total Renewable Energy Consumption
##  – 12 Total Primary Energy Consumption

## Read the .csv file and import it as a data.frame
cons <- read.csv("EIA_pecbs.csv")

## Keep only date, value and the Column_Order
cons <- subset(cons, select = c(YYYYMM, Value, Column_Order))

## Remove all data whose date ends with 13 (does not include 2013),
## which are yearly data before 1973
cons <- subset(cons, !grepl("^.+13$", YYYYMM))

## Remove January 2016 and start from January 1974
cons <- subset(cons, YYYYMM < 201601 & YYYYMM > 197312)

## Check for missing values
# Note: in the database, unknown values are denoted by "Not Available"
table(cons$Value == "Not Available")
# There are 252 missing values (mostly for renewable energies with 
# 0 consumption around)

## Set NAs equal to 0
## The assumption made here is that some unknown data (notably for the
## renewable energy sources) were due to a lack of information implied
## by a lack of consumption
cons$Value[cons$Value == "Not Available"] <- 0

## Coal consumption ##
cons.coal <- subset(cons, Column_Order == 1)
# extract beginning date in order to construct the ts object
#start.date <- cons.coal$YYYYMM[1]
#start.year <- as.integer(substring(start.date, 1, 4))
#start.month <- as.integer(substring(start.date, 5, 6))
# construct the ts object
cons.coal <- ts(as.numeric(as.character(cons.coal$Value)), 
                start = c(start.year, start.month), 
                frequency = 12)

## Petroleum consumption excluding biofuels
cons.pet <- subset(cons, Column_Order == 3)
# construct the ts object with the same starting date as coal
cons.pet <- ts(as.numeric(as.character(cons.pet$Value)),
               start = c(start.year, start.month),
               frequency = 12)


## Total Fossil Fuels Consumption (nuclear, petroleum, coal)
cons.fossil <- subset(cons, Column_Order == 4)
# construct the ts object with the same starting date as coal
cons.fossil <- ts(as.numeric(as.character(cons.fossil$Value)), 
                  start = c(start.year, start.month), 
                  frequency = 12)



## Total Renewable Energy Consumption
cons.renew <- subset(cons, Column_Order == 11)
# construct the ts object with the same starting date as coal
cons.renew <- ts(as.numeric(as.character(cons.renew$Value)),
                 start = c(start.year, start.month), 
                 frequency = 12)

## Total Primary Energy Consumption
cons.prim <- subset(cons, Column_Order == 12)
# construct the ts object with the same starting date as coal
cons.prim <- ts(as.numeric(as.character(cons.prim$Value)), 
                start = c(start.year, start.month),
                frequency = 12)

rm(cons)



# 2.2 Seasonality, growth rates, and aggregation --------------------------



# 2.2.1 Deseason of the energy consumption ------------------



## It is clear and logical that seasons have a strong impact on energy consumptions
## We remove the seasonality by assuming a linear seasonality and then aggregate as quarterly data
## The additive model is useful when the seasonal variation is relatively constant over time.
## Here, for energy consumption, we assume that we have a seasonal variation (weather and temperature).
## And that this season does not change over time, even if the variation might grow or shrink.
## We use monthly dummies, which are captured by the use of "season()" in dynlm().

## We also aggregate the data as quarterly data, since "gdp" is a quarterly data set


## Deseason coal consumption
ds.cons.coal.result <- dynlm(cons.coal ~ 1 + season(cons.coal)) 
cons.coal <- coefficients(ds.cons.coal.result)[1]                  + 
                sum((coefficients(ds.cons.coal.result))[2:12]/12)  +
                residuals(ds.cons.coal.result)
rm(ds.cons.coal.result)

## Deseason total fossil fuel consumption
ds.cons.fossil.result <- dynlm(cons.fossil ~ 1 + season(cons.fossil)) 
cons.fossil <- coefficients(ds.cons.fossil.result)[1]                  + 
                  sum((coefficients(ds.cons.fossil.result))[2:12]/12)  +
                  residuals(ds.cons.fossil.result)
rm(ds.cons.fossil.result)

## Deseason total petroleum consumption
ds.cons.pet.result <- dynlm(cons.pet ~ 1 + season(cons.pet)) 
cons.pet <- coefficients(ds.cons.pet.result)[1]                  + 
                     sum((coefficients(ds.cons.pet.result))[2:12]/12)  +
                      residuals(ds.cons.pet.result)
rm(ds.cons.pet.result)

## Deseason renewable energy consumption
ds.cons.renew.result <- dynlm(cons.renew ~ 1 + season(cons.renew)) 
cons.renew <- coefficients(ds.cons.renew.result)[1]                  + 
                 sum((coefficients(ds.cons.renew.result))[2:12]/12)  +
                 residuals(ds.cons.renew.result)
rm(ds.cons.renew.result)

## Deseason primary energy consumption
ds.cons.prim.result <- dynlm(cons.prim ~ 1 + season(cons.prim)) 
cons.prim <- coefficients(ds.cons.prim.result)[1]                  +
                sum((coefficients(ds.cons.prim.result))[2:12]/12)  +
                residuals(ds.cons.prim.result)
rm(ds.cons.prim.result)



# 2.2.2 Deseason the final purchase domestic crude oil price --------



## Deseason primary energy consumption
#ds.cop.result <- dynlm(cop ~ 1 + season(cop)) 
#ds.cop <- coefficients(ds.cop.result)[1]                  +
#  sum((coefficients(ds.cop.result))[2:12]/12)  +
#  residuals(ds.cop.result)
#rm(ds.cop.result)
#ds.cop <- aggregate(ds.cop, nfrequency = 4)

## With the above commented lines, we checked for seasonality.
## Plot both variables as:
## ts.plot(cop, ds.cop, col = c("black", "red"))
## Removing a "hypothetical" seasonality seems to actually include a new seasonality.
## Moreover since 1975, the Strategic Petroleum Reserve (SPR) is an emergency fuel 
## storage of oil maintained by the United States Department of Energy.
## It was implemented after the oil embargo in 1973 to mitigate future temporary
## supply disruptions. Thus, there is no effect of an increase in demand, since 
## the supply cannot be disrupted: this is our argument in favor of no seasonality
## of the crude oil price.
## It is also the fact that price of crude oil (not the final purchase domestic ones
## though) is mostly fixed by the OPEC countries (even if their supply power is 
## decreasing). Hence their seasonality is not a consumption seasonality that could 
## take place in the US.



# 2.2.3 Growth rates and aggregation --------------------------------------



## Since gdp is a quarterly data set, we need to aggregate the 
## other objects from monthly data to quarterly data.
## We also want to compute the growth rates. Use the function "growthrate" we wrote.



## Growth rate of crude oil price, mean of three months
g.cop <- growthrate(cop, agg = TRUE, freq = 4, bymean = TRUE)

## Growth rate of gross domestic product
g.gdp <- growthrate(gdp)

## Growth rate of energy consumptions
g.cons.coal   <- growthrate(cons.coal,   agg = TRUE, freq = 4)
g.cons.fossil <- growthrate(cons.fossil, agg = TRUE, freq = 4)
g.cons.pet    <- growthrate(cons.pet,    agg = TRUE, freq = 4)
g.cons.renew  <- growthrate(cons.renew,  agg = TRUE, freq = 4)
g.cons.prim   <- growthrate(cons.prim,   agg = TRUE, freq = 4)
g.cons.renew2001 <- window(g.cons.renew, start = c(2001,1))

# 2.3 Plot data -----------------------------------------------------------

# "#pdf()" and "#dev.off()" are used to save the plot when not running in batch mode

#pdf("US_energyconsumption.pdf", width = 9, height = 7.7)
par(mfrow = c(3,2), oma = c(2,2,1,1), mar =c(2,2,2,2), mgp = c(2,1,0), xpd = NA)
ts.plot(gdp, ylab = "Billion dollars", xlab = "", lwd = 0.5,
        main = "US gross domestic product (1974-2015)")
ts.plot(cop, ylab = "Dollars", xlab = "", lwd = 0.5,
        main = "US final purchase domestic crude oil price (1974-2015)")
ts.plot(cons.fossil, "ylab" = "Quadrillion Btu", xlab = "",
        main = "US consumption of fossil fuels (1974-2015)", lwd = 0.5)
ts.plot(cons.pet, "ylab" = "Quadrillion Btu", xlab = "",
        main = "US consumption of petroleum (1974-2015)", lwd = 0.5)
ts.plot(cons.renew, "ylab" = "Quadrillion Btu", xlab = "",
        main = "US consumption of renewable energy (1974-2015)", lwd = 0.5)
ts.plot(cons.prim, "ylab" = "Quadrillion Btu", xlab = "",
        main = "US consumption of primary energy (1974-2015)", lwd = 0.5)
par(mfrow = c(1,1))
#dev.off()

       



# 3 Unit root test and cointegration --------------------------------------




# 3.1 ADF-test ----------------------------------------------------------



## Perform the Augmented Dickey Fuller test
## Use the auto.adf function, loaded previously, that we have constructed
## All the tests are performed at a 0.05 significance level without any further precision
## The Lag selection is achieved according to the Akaike "AIC" (default in auto.adf())
## As a rule of thumb, we use a maximum lags of 3 times the frequency, that is 12

## Energy consumptions
auto.adf(g.cons.coal,      maxlag = 12)  # no unit root (trend and drift)
auto.adf(g.cons.pet,       maxlag = 12)  # no unit root (trend and drift)
auto.adf(g.cons.fossil,    maxlag = 12)  # no unit root (trend and drift)
auto.adf(g.cons.renew,     maxlag = 12)  # no unit root (trend and drift)
auto.adf(g.cons.prim,      maxlag = 12)  # no unit root (trend and drift)
auto.adf(g.cons.renew2001, maxlag = 12)  # no unit root (trend and drift)

## Gross domestic product
auto.adf(g.gdp, maxlag = 12)  # no unit root (trend and drift)

## Crude oil price
auto.adf(g.cop, maxlag = 12)  # no unit root (trend and drift)



# 3.2 PP-test -------------------------------------------------------------



## Perform the Phillips-Perron test
## Use the auto.pp function, loaded previously, that we have constructed
## All the tests are performed at a 0.05 significance level without any further precision
## Note: perform worse than the ADF in finite samples (Davidson and MacKinnon, 2004)

## Energy consumptions
auto.pp(g.cons.coal)      # no unit root (trend and drift)
auto.pp(g.cons.pet)       # no unit root (trend and drift)
auto.pp(g.cons.fossil)    # no unit root (trend and drift)
auto.pp(g.cons.renew)     # no unit root (trend and drift)
auto.pp(g.cons.prim)      # no unit root (trend and drift)
auto.pp(g.cons.renew2001) # no unit root (trend and drift)

## Gross domestic product
auto.pp(g.gdp)  # no unit root (trend and drift)

## Crude oil price
auto.pp(g.cop)  # no unit root (trend and drift)



# 3.3 KPSS-test -----------------------------------------------------------



## Perform the Kwiatkowski and al. test
## Use the auto.kpss function, loaded previously, that we have constructed
## All the tests are performed at a 0.05 significance level without any further precision

## Energy consumptions
auto.kpss(g.cons.coal)      # no unit root (trend and drift)
auto.kpss(g.cons.pet)       # no unit root (trend and drift)
auto.kpss(g.cons.fossil)    # no unit root (trend and drift)
auto.kpss(g.cons.renew)     # no unit root (trend and drift)
auto.kpss(g.cons.prim)      # no unit root (trend and drift)
auto.kpss(g.cons.renew2001) # no unit root (trend and drift)

## Gross domestic product
auto.kpss(g.gdp)  # no unit root (trend and drift)

## Crude oil price
auto.kpss(g.cop)  # no unit root (trend and drift)



# 3.4 Cointegration ------------------------------------------------------



## There is no need to test for cointegration since every variable is I(0).
## However, we would have performed the following tests.

## Among two variables, the Engel-Granger test:
## result <- dynlm(y ~ x)
## coef <- coeftest(result)
## coef
## We want now to test whether the residuals are stationary or not
## We cannot use the ADF-test, since it assumes that data are generated by
## a stochastic process.
## Thus we use the EG-ADF test
## EG <- ur.df(result$residuals, type = "drift", lags=12, selectlags="AIC")
## EG
## Reject the null if tau2 is lower than the following critical values in 
## in Woolridge table 18.4 are: 
#####################################################
###    Significance level   10%      5%      1%   ###
###    Critical value      -3.04    -3.34  -3.90  ###
#####################################################
## Note: the null hypothesis is: there is a unit root, which implies
## no cointegration.

## Among all the variables to be used in the VAR analysis, in order to obtain
## the number of cointegration relations
## cajo <- ca.jo(data, type = "trace", ecdet="trend", K=5)
## stats.cajo <- cbind(cajo@teststat, cajo@cval)
## stats.cajo




# 4 VAR models ------------------------------------------------------------




# 4.1 VAR primary energy --------------------------------------------------



## Merge the data in one ts object (same time scale)
var.prim <- ts.intersect(g.cons.prim, g.gdp)
## Selection criteria for the number of lags
var.select <- VARselect(var.prim, lag.max=12)$selection
## Get the minimum of the criteria and save it in variable p
criterium <- min(var.select)
## Note: all criteria suggest 3 lags

## Estimate the VAR with ols
var.prim.ols <- VAR(var.prim, p = criterium)

## Multivariate Portmanteau test
serial.test(var.prim.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.06414
## Null hypothesis of no autocorrelation is accepted.

## Multivariate test for normality in the residuals
normality.test(var.prim.ols, multivariate.only = F)$jb.uni
## p-values far smaller than 0.05
## Null of normality is rejected.

## Root test
roots(var.prim.ols)
## The roots are real and MOSTLY distinct and 
## in the positive part of the unit interval.
## The estimated model is stable.

## Granger-causality tests

estprim.hac <- grangertest(g.gdp ~ g.cons.prim, order = criterium, vcov = NeweyWest)
estprim.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## The p-value is 0.02503. The null is rejected at 0.05 significance level, thus 
## g.cons.prim granger-causes g.gdp

estpriminv.hac <- grangertest(g.cons.prim ~ g.gdp, order = criterium, vcov = NeweyWest)
estpriminv.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## The p-value is 3.059e-11. The null is rejected
## g.gdp granger-causes g.cons.prim



# 4.2 VAR fossil fuel -----------------------------------------------------



## Merge the data in one ts object (same time scale)
var.fossil <- ts.intersect(g.cons.fossil, g.gdp)
## Selection criteria for the number of lags
var.select <- VARselect(var.fossil, lag.max=12)$selection
## Get the minimum of the criteria and save it in variable p
criterium <- min(var.select)
## Note: all criteria suggest 3 lags

## Estimate the VAR with ols
var.fossil.ols <- VAR(var.fossil, p = criterium)

## Multivariate Portmanteau test
serial.test(var.fossil.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.118
## Null hypothesis of no autocorrelation is accepted at a 0.05 significance level

## Multivariate test for normality in the residuals
normality.test(var.fossil.ols, multivariate.only = F)$jb.uni
## p-values far smaller than 0.05
## Null of normality is rejected.

## Root test
roots(var.fossil.ols)
## The roots are real and MOSTLY distinct and 
## in the positive part of the unit interval.
## The estimated model is stable.

## Granger-causality tests

estfossil.hac <- grangertest(g.gdp ~ g.cons.fossil, order = criterium, vcov = NeweyWest)
estfossil.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## The p-value is 0.06189. The null is NOT rejected at 0.05 significance level, 
## but is rejected at a 0.1 significance level.
## Thus, for convenience (and the purpose of this term paper), 
## we will say that g.cons.fossil granger-causes g.gdp

estfossilinv.hac <- grangertest(g.cons.fossil ~ g.gdp, order = criterium, vcov = NeweyWest)
estfossilinv.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## The p-value is 3.085e-11. The null is rejected.



# 4.3 VAR coal ------------------------------------------------------------



## Merge the data in one ts object (same time scale)
var.coal <- ts.intersect(g.cons.coal, g.gdp)
## Selection criteria for the number of lags
var.select <- VARselect(var.coal, lag.max=12)$selection
## Get the minimum of the criteria and save it in variable p
criterium <- min(var.select)
## Note: AIC suggests 12 lags and Bayes 1

## Estimate the VAR with ols
var.coal.ols <- VAR(var.coal, p = criterium)

## Multivariate Portmanteau test
serial.test(var.coal.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.001586
## Null hypothesis of no autocorrelation is rejected


## We re-run the var model with one more lag and test again for autocorrelation


## Estimate the VAR with ols
var.coal.ols <- VAR(var.coal, p = criterium + 1)

## Multivariate Portmanteau test
serial.test(var.coal.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.001586
## Null hypothesis of no autocorrelation is rejected


## We re-run the var model with one more lag and test again for autocorrelation


## Estimate the VAR with ols
var.coal.ols <- VAR(var.coal, p = criterium + 2)

## Multivariate Portmanteau test
serial.test(var.coal.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.001586
## Null hypothesis of no autocorrelation is rejected at a 0.05 significance level
## but accepted at a 0.01 significance level
## criterium is 1
## So we needed to add 2 lags, that is to use 3 lags

## Multivariate test for normality in the residuals
normality.test(var.coal.ols, multivariate.only = F)$jb.uni
## p-values far smaller than 0.05
## Null of normality is rejected.

## Root test
roots(var.coal.ols)
## The roots are real and MOSTLY distinct and 
## in the positive part of the unit interval.
## The estimated model is stable.

## Granger-causality test
estcoal.hac <- grangertest(g.gdp ~ g.cons.coal, order = criterium + 2, vcov = NeweyWest)
estcoal.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## The p-value is 0.937. The null is NOT rejected at all, thus 
## g.cons.coal DO NOT granger-causes g.gdp
## This might be due to the fact that the use of coal decreases during the past 15 years.
## However, even while testing for a window from 1974Q1 to 1989Q4, the p-value decreases 
## to 0.8909, without rejecting the hypothesis.
## Check for our own information if g.gdp granger-causes g.cons.coal
estcoalinv.hac <- grangertest(g.cons.coal ~ g.gdp, order = criterium + 2, vcov = NeweyWest)
estcoalinv.hac
## g.gdp granger-causes g.cons.coal at a 0.01 significance level
## We can assume that as the GDP increased, as the price of extracting coal increased, and
## as politicians started to take care about pollution, it has decreased the use of coal.

## Hence, for further consideration, we would just use the shock from g.cons.coal to g.gdp



# 4.4 VAR pet -------------------------------------------------------------



## Merge the data in one ts object (same time scale)
var.pet <- ts.intersect(g.cons.pet, g.gdp)
## Selection criteria for the number of lags
var.select <- VARselect(var.pet, lag.max=12)$selection
## Get the minimum of the criteria and save it in variable p
criterium <- min(var.select)
## Note: AIC suggest 12, Schwarz 3, and HQ 4

## Estimate the VAR with ols
var.pet.ols <- VAR(var.pet, p = criterium)

## Multivariate Portmanteau test
serial.test(var.pet.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.0002837
## Null hypothesis of no autocorrelation is rejected


## We re-run the var model with one more lag and test again for autocorrelation


## Estimate the VAR with ols
var.pet.ols <- VAR(var.pet, p = criterium + 1)

## Multivariate Portmanteau test
serial.test(var.pet.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.00519
## Null hypothesis of no autocorrelation is rejected


## We re-run the var model with one more lag and test again for autocorrelation


## Estimate the VAR with ols
var.pet.ols <- VAR(var.pet, p = criterium + 2)

## Multivariate Portmanteau test
serial.test(var.pet.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.001253
## Null hypothesis of no autocorrelation is rejected


## We re-run the var model with one more lag and test again for autocorrelation


## Note: we re-do this procedure until "criterium + 6" which is 9 lags
## Estimate the VAR with ols
var.pet.ols <- VAR(var.pet, p = criterium + 6)

## Multivariate Portmanteau test
serial.test(var.pet.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.02038
## Null hypothesis of no autocorrelation is accepted at a 0.01 significance
## level, but rejected at a 0.05 significance level

## Multivariate test for normality in the residuals
normality.test(var.pet.ols, multivariate.only = F)$jb.uni
## Null of normality is rejected for g.gdp
## Nult of normality is accepted for g.cons.pet

## Root test
roots(var.pet.ols)
## The roots are real and MOSTLY distinct and 
## in the positive part of the unit interval.
## The estimated model is stable.
## Note: some roots are close to the unity (0.96, e.g.)

## Granger-causality test
estpet.hac <- grangertest(g.gdp ~ g.cons.pet, order = criterium + 6, vcov = NeweyWest)
estpet.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## The p-value is 2.706e-09. The null is rejected

## Granger-causality test
estpetinv.hac <- grangertest(g.cons.pet ~ g.gdp, order = criterium + 6, vcov = NeweyWest)
estpetinv.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## The p-value is 7.748e-3. The null is rejected



# 4.5 VAR renew -----------------------------------------------------------



## Merge the data in one ts object (same time scale)
var.renew <- ts.intersect(g.cons.renew, g.gdp)
## Selection criteria for the number of lags
var.select <- VARselect(var.renew, lag.max=12)$selection
## Get the minimum of the criteria and save it in variable p
criterium <- min(var.select)
## Note: all criteria suggest 1 lag

## Estimate the VAR with ols
var.renew.ols <- VAR(var.renew, p = criterium)

## Multivariate Portmanteau test
serial.test(var.renew.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.01389
## Null hypothesis of no autocorrelation is accepted at a 0.01 significance level
## but not at  a 0.05 significance level

## Multivariate test for normality in the residuals
normality.test(var.renew.ols, multivariate.only = F)$jb.uni
## Null of normality is rejected for g.gdp
## Null of normality is accepted for g.cons.renew

## Root test
roots(var.renew.ols)
## The roots are real and MOSTLY distinct and 
## in the positive part of the unit interval.
## The estimated model is stable.


## Granger-causality test
estrenew.hac <- grangertest(g.gdp ~ g.cons.renew, order = criterium, vcov = NeweyWest)
estrenew.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## The p-value is 0.063. The null is NOT rejected at a 0.1 significance level



## Would using a new window for recent years, where more renewable energy sources are used,
## yield a significant granger-causality test
## From the plot of cons.renew,  2001 is the increasing starting year.
## Window data from 2001Q1 to 2015Q4
## That is 4*15 = 60 observations, which is enough to estimate a VAR model (higher number than 50)
g.gdp2001 <- window(g.gdp, start = c(2001,1))
var.renew2001 <- ts.intersect(g.gdp, g.cons.renew2001)
## Selection criteria for the number of lags
var.select <- VARselect(var.renew2001, lag.max=12)$selection
## Get the minimum of the criteria and save it in variable p
criterium <- min(var.select)
## Note: all criteria but AIC suggest 1 lag, and AIC suggests 10 lags

## Estimate the VAR with ols
var.renew2001.ols <- VAR(var.renew2001, p = criterium)

## Multivariate Portmanteau test
serial.test(var.renew2001.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.1363
## Null hypothesis of no autocorrelation is rejected at a 0.05 significance level, 
## but not at a 0.01 significance level
## The null is rejected until criterium + 2, i.e., 3 lags.

## Thus the VAR model we estimate would be with "p = criterium + 2", which gives:

## Estimate the VAR with ols
var.renew2001.ols <- VAR(var.renew2001, p = criterium + 2)

## Multivariate Portmanteau test
serial.test(var.renew2001.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.1718
## Null hypothesis of no autocorrelation is accepted


## Multivariate test for normality in the residuals
normality.test(var.renew2001.ols, multivariate.only = F)$jb.uni
## Null of normality is rejected for g.gdp
## Nult of normality is accepted for g.cons.renew2000

## Root test
roots(var.renew2001.ols)
## The roots are real and MOSTLY distinct and 
## in the positive part of the unit interval.
## The estimated model is stable.


## Granger-causality tests

estrenew2001.hac <- grangertest(g.gdp2001 ~ g.cons.renew2001, order = criterium + 2, vcov = NeweyWest)
estrenew2001.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## The p-value is 0.01274. The null is rejected at a 0.05 significance level
## Thus, g.cons.renew granger-causes g.gdp for recent years, not for the entire period.

estrenew2001inv.hac <- grangertest(g.cons.renew2001 ~ g.gdp2001, order = criterium + 2, vcov = NeweyWest)
estrenew2001inv.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## The p-value is 6.447e-3. The null is rejected.



# 4.6 VAR renew and fossil  -----------------------------------------------



## Granger-causality test
## Tests were done for g.cons.renew and g.cons.fossil
## Recall:
## – g.cons.fossil granger-causes g.gdp at a 0.1 significance level
## – g.cons.renew DOES NOT granger-cause g.gdp
## – BUT g.cons.renew2001 granger-cause g.gdp at a 0.05 significance level


## For the entire period from 1974Q1 to 2015Q4


## Merge the data in one ts object (same time scale)
var.renewfossil <- ts.intersect(g.gdp, g.cons.renew, g.cons.fossil)
## Selection criteria for the number of lags
var.select <- VARselect(var.renewfossil, lag.max=12)$selection
## Get the minimum of the criteria and save it in variable p
criterium <- min(var.select)
## Note: all criteria suggest 3 lags

## Estimate the VAR with ols
var.renewfossil.ols <- VAR(var.renewfossil, p = criterium)

## Multivariate Portmanteau test
serial.test(var.renewfossil.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.08769
## Null hypothesis of no autocorrelation is accepted at a 0.05 significance level

## Multivariate test for normality in the residuals
normality.test(var.renewfossil.ols, multivariate.only = F)$jb.uni
## Null of normality is rejected for g.gdp
## Null of normality is accepted for g.cons.renew
## Null of normality is rejected for g.cons.fossil

## Root test
roots(var.renewfossil.ols)
## The two roots are real and distinct and 
## in the positive part of the unit interval.
## The estimated model is stable.


## For the recent period from 2001Q1 to 2015Q4

## Merge the data in one ts object (same time scale)
var.renew2001fossil <- ts.intersect(g.gdp, g.cons.renew2001, g.cons.fossil)
## Selection criteria for the number of lags
var.select <- VARselect(var.renew2001fossil, lag.max=12)$selection
## Get the minimum of the criteria and save it in variable p
criterium <- min(var.select)
## Note: all criteria suggest 1 lag, but AIC which suggest 12

## Estimate the VAR with ols
var.renew2001fossil.ols <- VAR(var.renew2001fossil, p = criterium)

## Multivariate Portmanteau test
serial.test(var.renew2001fossil.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.1203
## Null hypothesis of no autocorrelation is accepted 

## Multivariate test for normality in the residuals
normality.test(var.renew2001fossil.ols, multivariate.only = F)$jb.uni
## Null of normality is rejected for g.gdp
## Null of normality is accepted for g.cons.renew
## Null of normality is accepted for g.cons.fossil

## Root test
roots(var.renew2001fossil.ols)
## The roots are real and distinct and 
## in the positive part of the unit interval.
## The estimated model is stable.


# Granger causality test of g.cons.fossil for the new period:
g.cons.fossil2001 <- window(g.cons.fossil, start= c(2001,1))

estfossil2001.hac <- grangertest(g.gdp2001 ~ g.cons.fossil2001, order = criterium +2, vcov = NeweyWest)
estfossil2001.hac
## p-value = .3002. The null hypothesis of DO NOT granger-cause is accepted.

## The problem here is:
## for the time period 1974Q1-2015Q4, fossil granger-causes gdp, but not renew
## for the time period 2001Q1–2015Q4, renew granger-causes gdp, but not fossil



# 4.7 VAR renew and pet ---------------------------------------------------



## We test this VAR model for the recent years

## Merge the data in one ts object (same time scale)
var.renew2001pet <- ts.intersect(g.gdp, g.cons.renew2001, g.cons.pet)
## Selection criteria for the number of lags
var.select <- VARselect(var.renew, lag.max=12)$selection
## Get the minimum of the criteria and save it in variable p
criterium <- min(var.select)
## Note: HQ and SC suggest 1 lag, and AIC and FPE suggest respectively 12 and 9 lags

## Estimate the VAR with ols
var.renew2001pet.ols <- VAR(var.renew2001pet, p = criterium +2)

## Multivariate Portmanteau test
serial.test(var.renew2001pet.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.1714
## Null hypothesis of no autocorrelation is accepted after adding 2 lags to the criterium
## Hence, number of lags to be used later on is 3.

## Multivariate test for normality in the residuals
normality.test(var.renew2001pet.ols, multivariate.only = F)$jb.uni
## Null of normality is accepted for g.gdp
## Null of normality is accepted for g.cons.renew
## Null of normality is accepted at a 0.01 significane level, not at a 0.05 SL

## Root test
roots(var.renew2001pet.ols)
## The roots are real and MOSTLY distinct and 
## in the positive part of the unit interval.
## The estimated model is stable.


## Granger-causality tests

## Granger-causality test for g.cons.renew2001
## Test already done with 3 lags
## Recall: g.cons.renew2001 granger-causes g.gdp at a 0.05 significance level

## Granger-causality test for g.cons.pet2001
g.cons.pet2001 <- window(g.cons.pet, start = c(2001,1))
estrenew2001pet.hac <- grangertest(g.gdp2001 ~ g.cons.pet2001, order = criterium, vcov = NeweyWest)
estrenew2001pet.hac
estrenew2001pet.hac <- grangertest(g.gdp2001 ~ g.cons.pet2001, order = criterium + 2, vcov = NeweyWest)
estrenew2001pet.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## 1 lag = crit: p-value is 0.081111. The null is rejected at a 0.1 significance level
## 3 lags: p-value is 0.1301. The null is NOT rejected
## Thus if we control for autocorrelation, that is lags = 3:
## g.cons.renew2001 granger-causes g.gdp, otherwise no (see previous tests); and
## we will not use this model for the SVAR/IRF



# 4.8 VAR primary energy and oil price ------------------------------------



## Merge the data in one ts object (same time scale)
var.primcop <- ts.intersect(g.gdp, g.cons.prim, g.cop)
## Selection criteria for the number of lags
var.select <- VARselect(var.primcop, lag.max=12)$selection
## Get the minimum of the criteria and save it in variable p
criterium <- min(var.select)
## Note: all criteria suggest 3 lags

## Estimate the VAR with ols
var.primcop.ols <- VAR(var.primcop, p = criterium)

## Multivariate Portmanteau test
serial.test(var.primcop.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.1634
## Null hypothesis of no autocorrelation is accepted at a 0.05 significance level

## Multivariate test for normality in the residuals
normality.test(var.primcop.ols, multivariate.only = F)$jb.uni
## Null of normality is rejected for all variables

## Root test
roots(var.primcop.ols)
## The roots are real and MOSTLY distinct and 
## in the positive part of the unit interval.
## The estimated model is stable.


## Granger-causality tests

## Granger-causality for primary ernergy
## Done earlier with the same number of lags
## Recall: g.cons.prim granger-causes g.gdp



## Granger-causality tests for crude oil price

estcop.hac <- grangertest(g.gdp ~ g.cop, order = criterium , vcov = NeweyWest)
estcop.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## The p-value is 0.02503. The null is accepted at a 0.05 significance level,
## but rejected at a 0.1 significance level.
## Thus, for convenience (and the purpose of this term paper), 
## we will say that g.cop granger-causes g.gdp

estcopinv.hac <- grangertest(g.cop ~ g.gdp, order = criterium , vcov = NeweyWest)
estcopinv.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## The p-value is 0.5079. The null is accepted.
## If one chooses this model, one should not look at the effect of g.gdp on g.cop
## One could think of g.cop as an exogenous parameter not caused by an increase in 
## the US gdp.



# 4.9 VAR renew and oil price ---------------------------------------------



## Merge the data in one ts object (same time scale)
var.renew2001cop <- ts.intersect(g.gdp, g.cons.renew2001, g.cop)
## Selection criteria for the number of lags
var.select <- VARselect(var.renew2001cop, lag.max=12)$selection
## Get the minimum of the criteria and save it in variable p
criterium <- min(var.select)
## Note: all criteria but SC (1 lag) suggest 12 lags

## Estimate the VAR with ols
var.renew2001cop.ols <- VAR(var.renew2001cop, p = criterium)

## Multivariate Portmanteau test
serial.test(var.renew2001cop.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.3634
## Null hypothesis of no autocorrelation is accepted 

## Multivariate test for normality in the residuals
normality.test(var.renew2001cop.ols, multivariate.only = F)$jb.uni
## Null of normality is rejected for g.gdp
## Null of normality is accepted for g.cop
## Null of normality is accepted for g.cons.renew2001 at a 0.01 significance level

## Root test
roots(var.renew2001cop.ols)
## The roots are real and distinct and 
## in the positive part of the unit interval.
## The estimated model is stable.

## Granger-causality tests

## Granger-causality test for g.renew2001
## We previously tested the granger causality for 3 lags
## Here, criterium is 1.
estrenew2001.hac <- grangertest(g.gdp2001 ~ g.cons.renew2001, order = criterium, vcov = NeweyWest)
estrenew2001.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## The p-value is 0.9758 The null is rejected.
## g.cons.renew2001 DOES NOT granger-cause g.gdp with only one lag

## Granger-causality test for g.cop2001
g.cop2001 <- window(g.cop, start = c(2001,1))
estcop2001.hac <- grangertest(g.gdp2001 ~ g.cop2001, order = criterium, vcov = NeweyWest)
estcop2001.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## The p-value is 0.6384. The null is accepted 
## g.cop2001 DOES NOT granger-cause g.gdp 

## Hence we do not test for the inverse granger-causality, since we would not use this model anyway



# 4.10 VAR pet and oil price -----------------------------------------------



## Merge the data in one ts object (same time scale)
var.petcop <- ts.intersect(g.gdp, g.cons.pet, g.cop)
## Selection criteria for the number of lags
var.select <- VARselect(var.petcop, lag.max=12)$selection
## Get the minimum of the criteria and save it in variable p
criterium <- min(var.select)
## Note: AIC and FPE suggest 5 lags, the others 3

## Estimate the VAR with ols
var.petcop.ols <- VAR(var.petcop, p = criterium + 1)

## Multivariate Portmanteau test
serial.test(var.petcop.ols, lags.pt = 12, type = "PT.asymptotic")
## p-value = 0.1634
## Null hypothesis of no autocorrelation is accepted at a 0.05 significance level

## Multivariate test for normality in the residuals
normality.test(var.petcop.ols, multivariate.only = F)$jb.uni
## Null of normality is rejected for all variables

## Root test
roots(var.petcop.ols)
## The roots are real and MOSTLY distinct and 
## in the positive part of the unit interval.
## The estimated model is stable.

## Granger-causality tests

## Granger-causality for petroleum
## Done earlier with the same number of lags
## Recall: g.cons.pet granger-causes g.gdp
## Recall: g.gdp granger-causes g.cons.pet

## Granger-causality for crude oil price
estcop.hac <- grangertest(g.gdp ~ g.cop, order = criterium  + 1 , vcov = NeweyWest)
estcop.hac
estpet.hac <- grangertest(g.gdp ~ g.cons.pet, order = criterium + 1 , vcov = NeweyWest)
estpet.hac
## Null hypothesis for the F-test: DO NOT granger cause 
## g.cop DOES NOT granger-cause g.gdp
## g.pet granger-causes g.gdp at a 0.05 significance level
estpetcop.hac <- grangertest(g.cons.pet ~ g.cop, order = criterium  +1, vcov = NeweyWest)
estpetcop.hac
## g.cop granger-causes g.gdp at a 0.05 significance level
## It would be here possible to use g.cop as an instrument variable for g.pet
## We will use this VAR model by further arguing that the p-value of 0.1642, which is
## far from a neat significance level, still allow us to derive the SVAR and IRF.
## Also, for the matter of this term paper, we state that this is enough of a p-value
## not that bad, since we want to look at the IRFs.



# 4.11 Summary on VAR models ----------------------------------------------



## GDP ~ Primary energy
## g.cons.prim granger-causes g.gdp (0.05 SL)
## g.gdp granger-causes g.gdp (strongly)
## Use this model

## GDP ~ Fossil Fuel
## g.cons.fossil granger-causes g.gdp (0.1 SL)
## g.gdp granger-causes g.gdp (strongly)
## Use this model

## GDP ~ Coal
## g.cons.coal DOES NOT granger-causes g.gdp
## g.gdp granger-causes g.cons.coal
## DO NOT use this model

## GDP ~ Petroleum 
## g.cons.pet granger-causes g.gdp (strongly)
## g.gdp granger-causes g.cons.pet (strongly)
## Use this model

## GDP ~ Renewable energies
## g.cons.renew2001 granger-causes g.gdp2001 (0.05 SL)
## g.gdp2001 granger-causes g.cons.renew2001 (0.1 SL)
## Use this model
## Note: possibility of IRF comparison between g.cons.renew and g.cons.renew2001

## GDP ~ Renewable energies + Fossil fuel
## Both DO NOT granger-cause g.gdp for the same time period
## Both are granger-caused for the same time period
## DO NOT use this model

## GDP ~ Renewable energies + Petroleum
## DO NOT use this model (see section 4.7)

## GDP ~ Primary energy + Crude oil price
## g.cons.prim granger-causes g.gdp (0.05 SL)
## g.gdp granger-causes g.gdp (strongly)
## g.cop granger-causes g.gdp (0.1 SL)
## g.gdp DOES NOT granger-causes g.cop
## Use this model (carefully)

## GDP ~ Renewable energies + Crude oil price (period from 2001Q1)
## Both g.cons.renew2001 and g.cop2001 do not cause g.gdp2001, with the number
## of lags selected
## DO NOT use this model

## GDP ~ Petroleum + Crude oil price
## Use this model (see section 4.10)




# 5 SVAR and IRF ----------------------------------------------------------




## Shut down the current device
##dev.off()
## Note: cannot be used in batch mode
## Note: same for the following functions such as "#pdf()" or "#dev.off()"



# 5.1 GDP ~ Primary energy ------------------------------------------------



## We do not set up the restriction(s) and use a SVAR model since well-ordering
## the variables in the VAR model allows to use the Cholesky decomposition directly
## implemented with function "irf {vars}".

## Assume that g.cons.prim does not respond contemporaneously to a shock in
## g.gdp (that is also the right assumption since we want the effect on g.gdp)

## Use "irf.plot" to plot the IRFs and their confidence intervals
## Plot and save (not in batch mode) for 12 quarters ahead, that is for 3 years
#pdf("irf.prim.pdf", width = 10, height = 7) 
irf.plot.varest(var.prim.ols, n.ahead = 12)
#dev.off()



# 5.2 GDP ~ Fossil fuel ---------------------------------------------------



## We do not set up the restriction(s) and use a SVAR model since well-ordering
## the variables in the VAR model allows to use the Cholesky decomposition directly
## implemented with function "irf {vars}".

## Assume that g.cons.fossil does not respond contemporaneously to a shock in
## g.gdp (that is also the right assumption since we want the effect on g.gdp)

## Plot for 12 quarters ahead, that is for 3 years
#pdf("irf.fossil.pdf", width = 10, height = 7) 
irf.plot.varest(var.fossil.ols, n.ahead = 12)
#dev.off()



# 5.3 GDP ~ Petroleum -----------------------------------------------------



## We do not set up the restriction(s) and use a SVAR model since well-ordering
## the variables in the VAR model allows to use the Cholesky decomposition directly
## implemented with function "irf {vars}".

## Assume that g.cons.pet does not respond contemporaneously to a shock in
## g.gdp (that is also the right assumption since we want the effect on g.gdp)

## Plot for 12 quarters ahead, that is for 3 years
#pdf("irf.pet.pdf", width = 10, height = 7) 
irf.plot.varest(var.pet.ols, n.ahead = 12)
#dev.off()



# 5.4 GDP ~ Renewable energies --------------------------------------------



## We do not set up the restriction(s) and use a SVAR model since well-ordering
## the variables in the VAR model allows to use the Cholesky decomposition directly
## implemented with function "irf {vars}".

## Assume that g.cons.renew does not respond contemporaneously to a shock in
## g.gdp (that is also the right assumption since we want the effect on g.gdp)

## Plot for 12 quarters ahead, that is for 3 years
#pdf("irf.renew.pdf", width = 10, height = 7) 
irf.plot.varest(var.renew.ols, n.ahead = 12)
#dev.off()



# 5.5 GDP ~ Renewable energies (from 2001Q1) ------------------------------



## We do not set up the restriction(s) and use a SVAR model since well-ordering
## the variables in the VAR model allows to use the Cholesky decomposition directly
## implemented with function "irf {vars}".

## Assume that g.cons.renew2001 does not respond contemporaneously to a shock in
## g.gdp (that is also the right assumption since we want the effect on g.gdp)

## Plot for 12 quarters ahead, that is for 3 years
#pdf("irf.renew2001.pdf", width = 10, height = 7) 
irf.plot.varest(var.renew2001.ols, n.ahead = 12, zero = FALSE)
#dev.off()



# 5.6 GDP ~ Primary energy + Crude oil price ------------------------------



## We do not set up the restriction(s) and use a SVAR model since well-ordering
## the variables in the VAR model allows to use the Cholesky decomposition directly
## implemented with function "irf {vars}".

## Set up the restriction(s) for the SVAR model
## Since the restriction are different than with the Cholesky decomposition,
## we set up a SVAR model with our own restriction and then use irf.plot.svarest
Arestric <- matrix(NA, 3,3)
## Here we need 3*2/2 = 3 restrictions
## Since g.gdp DOES NOT granger-causes g.cop, this is our first restriction
Arestric[3,1] <- 0
## Assume that g.cons.prim does not respond contemporaneously to a shock in
## g.gdp (that is also the right assumption since we want the effect on g.gdp)
Arestric[2,1] <- 0
## Assume that g.cop does not respond contemporaneously to a shock in g.cons.prim. That
## is probably true since g.cons.prim is an aggregate of all kind of energy sources.
## Thus such a shock would be diluted, in the beginning at least.
Arestric[3,2] <-0

svar.primcop <- SVAR(var.primcop.ols, Amat = Arestric, Bmat = NULL, max.iter = 10000)
## Note:  The A-model is just identified. No test possible.
svar.primcop

## Plot for 12 quarters ahead, that is for 3 years
#pdf("irf.primcop.pdf", width = 15, height = 10) 
irf.plot.svarest(svar.primcop, n.ahead = 12)
#dev.off()



# 5.7 GDP ~ Petroleum + Crude oil price -----------------------------------



## Set up the restriction(s) for the SVAR model
Arestric <- matrix(NA, 3,3)
## Here we need 3*2/2 = 3 restrictions
## Assume that g.cons.prim and g.cop do not respond contemporaneously to a shock in
## g.gdp (that is also the right assumption since we want the effect on g.gdp)
Arestric[2,1] <- 0
Arestric[3,1] <- 0
## Assume that g.cop does not respond contemporaneously to a shock in g.cons.pet (see
## section 2.2.2 for an explanation)
Arestric[3,2] <-0

svar.petcop <- SVAR(var.petcop.ols, Amat = Arestric, Bmat = NULL, max.iter = 10000)
## Note:  The A-model is just identified. No test possible.
svar.petcop

## Plot for 12 quarters ahead, that is for 3 years
#pdf("irf.petcop.pdf", width = 10, height = 6.7) 
irf.plot.svarest(svar.petcop, n.ahead = 12)
#dev.off() 