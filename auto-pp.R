

###             AUTO.PP             ###
### Automatic function for pp.test  ###
### to be used only with ts objects ###

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