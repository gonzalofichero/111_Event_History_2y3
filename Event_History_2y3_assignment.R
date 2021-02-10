#### Loading libraries ####
library(tidyverse)
library(survival)
library(stargazer)



#### Exercise 1 ####


##### Importing data #####
breast <- read.table("breast.txt")

glimpse(breast)


##### Q1. Estimate a proportional hazard model with a Weibull distribution as baseline. #####

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
                      control=list(maxit=10^3), hessian=T)


# Optimized parameters
a.hat <- opt.weibull$par[1]
b.hat <- opt.weibull$par[2]
beta1.hat <- opt.weibull$par[3]
beta2.hat <- opt.weibull$par[4]



##### Q2. Provide a table with estimated coefficients and 95% confidence interval. Are the three group statistically diferent? #####

V <- solve(opt.weibull$hessian)


# CI for a and b (assuming 1-alpha=95%)
se.a <- sqrt(diag(V))[1]
se.b <- sqrt(diag(V))[2]

# CI for \beta_1 and \beta_2 (assuming 1-alpha=95%)
se.beta1 <- sqrt(diag(V))[3]
se.beta2 <- sqrt(diag(V))[4]

# Confidence intervals for all parameters
c(a.hat - 1.96*se.a , a.hat + 1.96*se.a) 
c(b.hat - 1.96*se.b , b.hat + 1.96*se.b)
c(beta1.hat - 1.96*se.beta1 , beta1.hat + 1.96*se.beta1)
c(beta2.hat - 1.96*se.beta2 , beta2.hat + 1.96*se.beta2)



##### Q3. Evaluate and plot the log-hazard for the three prognostic groups over the following times: ##### 

t <- seq(min(breast$exit), max(breast$exit), length=100)

# Creating the loghazard for the specific found parameters
loghazard <- function(t,x1,x2) {
  ## regression coefficients
  beta <- c(beta1.hat, beta2.hat)
  ## design matrix
  X <- cbind(x1, x2)
  ## regression part
  Xbeta <- X%*%beta
  
  lh <- (log(a.hat)-log(t)+a.hat*log(t)-a.hat*log(b.hat)+Xbeta)
  #lh <- (log(a.hat)+b.hat*t+Xbeta)
  return(lh)
}
# Check
loghazard(2, 0, 0)

# Creating data frame of simulation
data.frame(cbind(t, rep(0,100), rep(0,100), rep("Good", 100))) -> good.w
names(good.w) <- c("t", "x1", "x2", "group")

data.frame(cbind(t, rep(1,100), rep(0,100), rep("Medium", 100))) -> medium.w
names(medium.w) <- c("t", "x1", "x2", "group")

data.frame(cbind(t, rep(0,100), rep(1,100), rep("Poor", 100))) -> poor.w
names(poor.w) <- c("t", "x1", "x2", "group")


# All together
rbind(good.w, medium.w, poor.w) -> simul.weibull

# All to numeric in new dataframe
indx <- sapply(simul.weibull, is.factor)
simul.weibull.2 <- cbind(data.frame(lapply(simul.weibull[indx,1:3], function(x) as.numeric(as.character(x)))),simul.weibull$group)


# Running the function of the loghazard
for (i in 1:300){
  simul.weibull.2$logh[i] <- loghazard(simul.weibull.2$t[i], simul.weibull.2$x1[i], simul.weibull.2$x2[i])
}
names(simul.weibull.2) <- c("t", "x1", "x2", "group", "logh")

# Plotting
simul.weibull.2 %>% 
  ggplot(aes(x=t, y=logh, color=group)) + geom_line() +
  xlim(1.5,2) + ylim(100,150)



#### Exercise 2 ####

##### Importing data #####
students <- read.table("students.txt", header = T)

glimpse(students)

# Setting Event to numeric
students$event <- as.numeric(students$event)

##### Q1. Tabulate a frequency distribution of the duration variable or show this information in a barplot. Are tied observations rare in this data set? #####
students %>% 
  ggplot(aes(x=dur)) + geom_histogram(bins = max(students$dur)) +
  facet_grid(~subj)


##### Q2. Estimate a discrete PH model #####

##### Q2.A. Create a new variable quarter, which measures intervals of 3 months' length #####
students$quarter <- cut(students$dur, breaks = seq(0,24,by=3),
                        include.lowest = TRUE, labels = F)


##### Q2.B. Set up the data set appropriately so that you can estimate a discrete-time PH model. Provide the dimensions of the new augmented dataset. #####

students$id <- 1:length(students$dur)

# Taking from Giancarlo's code. But there has to be a better way of doing it...
i <- 1
st2 <- cbind(rep(students$id[i], students$quarter[i]),
             rep(students$quarter[i], students$quarter[i]),
             rep(students$sex[i], students$quarter[i]),
             rep(students$subj[i], students$quarter[i]),
             1:students$quarter[i],
             c(rep(0, students$quarter[i]-1), students$event[i]))

for(i in 2:nrow(students)){
  st2 <- rbind(st2,
                 cbind(rep(students$id[i], students$quarter[i]),
                       rep(students$quarter[i], students$quarter[i]),
                       rep(students$sex[i], students$quarter[i]),
                       rep(students$subj[i], students$quarter[i]),
                       1:students$quarter[i],
                       c(rep(0, students$quarter[i]-1), students$event[i])
                 ))
}


st3 <- as.data.frame(st2)
names(st3) <- c("id","quarter","sex","subj","quarter","event")
head(st3)

# Check
glimpse(st3)


##### Q3. Estimate the model where we have a discrete unspecified baseline and where sex and subject interact #####

# Estimate both the logistic and the cloglog version of the model and create a table for comparing the estimated coefficients.

# Factorizing linear predictors

st3$subj <- as.factor(st3$subj)
st3$sex <- as.factor(st3$sex)
st3$quarter <- as.factor(st3$quarter)

###### Logistic Regression ######
logistic <- glm(event ~ quarter + sex*subj,
                family = binomial,
                data = st3)

###### CLogLog Regression ######
cloglog <- glm(event ~ quarter + sex*subj,
               family = binomial(link=cloglog),
               data = st3)

summary(logistic)
summary(cloglog)


######  Creating table of results ###### 
stargazer(logistic, cloglog)

