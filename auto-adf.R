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
      msg = "'alpha' is the significance level and  must be either 0.01, 0.05, or 0.1",
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