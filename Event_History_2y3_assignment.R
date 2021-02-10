#### Loading libraries ####
library(tidyverse)
library(survival)


#### Importing data ####
breast <- read.table("breast.txt")

glimpse(breast)


#### Q1. Estimate a proportional hazard model with a Weibull distribution as baseline. ####

# Use patient with group=Good as reference
breast %>% 
  mutate(x1 = as.factor(ifelse(group == "Medium", 1, 0)),
         x2 = as.factor(ifelse(group == "Poor", 1, 0))) -> breast

# Generating a logL for Weibull: 2 parameters + 2 covariates (binary)

# Loglikelihood function with right-censoring
logL.W <- function(ini_par, event, observation, x1, x2){
  ## Weibull parameters
  a <- ini_par[1]
  b <- ini_par[2]
  ## regression coefficients
  beta <- ini_par[3:4]
  ## design matrix
  X <- cbind(x1, x2)
  ## regression part
  Xbeta <- X%*%beta
  expXbeta <- exp(Xbeta)
  
  ## In right-censoring we have the following LogL:
  # \sum_{i=1}^{n} (\delta_i * log(h(y_i)) - H(y_i))
  # here we have the event in the exact different position: so, we change delta to (1-delta)
  # Being H(t) = (t/b)^a
  
  # So, first we have the observations for those censored:
  lht <- (event * (log(a)-log(observation)+a*log(observation)-a*log(b)+Xbeta))
  # Now, the hazard for those that died:
  H0y <- (((observation/b)^a) * expXbeta)
  # Everything together
  logL <-   -sum( lht - H0y )
  return(logL)
}

# Setting initial parameters
ini_par <- c(1,10,1,1)

# Optimization
opt.weibull <- optim(ini_par, logL.W, event=breast$event,
                      observation=breast$exit,
                      x1=breast$x1, x2=breast$x2,
                      control=list(maxit=10^3))


