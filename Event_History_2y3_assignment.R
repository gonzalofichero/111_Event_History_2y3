#### Loading libraries ####
library(tidyverse)
library(survival)
library(stargazer)
library(survminer)
library(wesanderson)



#### Exercise 1 ####


##### Importing data #####
breast <- read.table("breast.txt")

glimpse(breast)


##### Q1. Estimate a proportional hazard model with a Weibull distribution as baseline. #####

# Use patient with group=Good as reference
breast %>% 
  mutate(x1 = ifelse(group == "Medium", 1, 0),
         x2 = ifelse(group == "Poor", 1, 0)) -> breast

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
  ggplot(aes(x=t, y=logh, color=group)) + geom_line()



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
names(st3) <- c("id","quarter","sex","subj","quarter2","event")
head(st3)

# Check
glimpse(st3)


##### Q2.C. Estimate the model where we have a discrete unspecified baseline and where sex and subject interact #####

# Estimate both the logistic and the cloglog version of the model and create a table for comparing the estimated coefficients.

# Factorizing linear predictors

st3$subj <- as.factor(st3$subj)
st3$sex <- as.factor(st3$sex)
st3$quarter2 <- as.factor(st3$quarter2)

###### Logistic Regression ######
logistic <- glm(event ~ quarter2 + sex*subj,
                family = binomial,
                data = st3)

###### CLogLog Regression ######
cloglog <- glm(event ~ quarter2 + sex*subj,
               family = binomial(link=cloglog),
               data = st3)

summary(logistic)
summary(cloglog)


######  Creating table of results ###### 
stargazer(logistic, cloglog, type="html")


##### Q3. The data were generated by a baseline hazard Weibull(2,15). #####
# Plot the estimated discrete baseline from the previous model (estimated using both approaches) and the true one in the same figure.

# Plot - Weibull!
x <- 1:8
h1 <- c(logistic$coef[1],
        logistic$coef[1] + logistic$coef[2:8])
h1 <- exp(h1) / (1 + exp(h1))

h2 <- c(cloglog$coef[1],
        cloglog$coef[1] + cloglog$coef[2:8])
h2 <- 1 - exp(-exp(h2))




plot(x,h1,col=1, pch=1,
     #ylim = c(0,0.55),
     xlab = "time", ylab="discrete hazard",main="hazard comparison")
points(h2,col=2,pch=2)



#### Exercise 3 ####

##### Importing data #####
swed_dead <- as.data.frame(read.table("SWIdeaths.txt", header = T))
swed_expo <- as.data.frame(read.table("SWIexposures.txt", header = T))


glimpse(swed_dead)
glimpse(swed_expo)

##### Data Wrangling #####
swed_dead.2 <- swed_dead %>% 
  pivot_longer(cols = starts_with("X"), names_to = "Year", names_prefix = "X", values_to = "Deaths", values_drop_na = F) 

swed_expo.2 <- swed_expo %>% 
  pivot_longer(cols = starts_with("X"), names_to = "Year", names_prefix = "X", values_to = "Exposure", values_drop_na = F) 

# Check NAs
sum(is.na(swed_dead.2))
sum(is.na(swed_expo.2))

# All together now!
swed.data <- left_join(swed_dead.2, swed_expo.2, by = c("Age","Year"))

glimpse(swed.data)


##### Q1. Fit a Makeham model on data in 1950 from age 30 to age 90 #####

# LogLikelihood function for Makeham
lnLmakeham <- function(para, Age, Death, Exposure){
  # parameters
  a <- para[1]
  b <- para[2]
  c <- para[3]
  mu <- (c + a * exp(b * Age))
  llk <- sum(Death * log(mu) - mu * Exposure)
  return(llk)
}


# Initial point
para <- c(0.00002, 0.1, 0.0001)
set.seed(42)

# Let's just keep the 30 to 90 years old people, only 1950
swed.data.3090 <- swed.data %>% filter(Age >= 30, Age <= 90, Year == "1950") %>% select(-Year)


# Optimization
make.mle <- optim(para,
                  lnLmakeham,
                  Age=swed.data.3090$Age, Death=swed.data.3090$Deaths, 
                  Exposure=swed.data.3090$Exposure, 
                  control=list(fnscale=-1),
                  #hessian = T,
                  method="Nelder-Mead")

make.mle

##### Q2. Create a table with estimated parameters from Q1 and associated 95% confidence intervals #####

# Optimized a is
a_hat <- make.mle$par[1]

# Optimized b is
b_hat <- make.mle$par[2]

# Optimized b is
c_hat <- make.mle$par[3]


# Build variance-covariance from Hessian
# Taking the Hessian with nlme library
library(nlme)

H <- fdHess(pars=make.mle$par, 
            fun=lnLmakeham,
            Age=swed.data.3090$Age, 
            Death=swed.data.3090$Deaths, 
            Exposure=swed.data.3090$Exposure)$Hessian


V <- solve(-H) # inverse of negative Hessian


#  square root of diagonal elements are the s.e. of a and b and c
sqrt(diag(V))
# parameter estimates correlated?
cov2cor(V)

# CI for a and b (assuming 1-alpha=95%)
se.a <- sqrt(diag(V))[1]
se.b <- sqrt(diag(V))[2]
se.c <- sqrt(diag(V))[3]

# Confidence intervals for both a and b (estimated) for Weibull function
c(a_hat - 1.96*se.a , a_hat + 1.96*se.a) 
c(b_hat - 1.96*se.b , b_hat + 1.96*se.b) 
c(c_hat - 1.96*se.c , c_hat + 1.96*se.c)



##### Q3. Plot estimated rates in 1950 along with the fitted rates from !1 and the 95% confidence interval #####

# Taking the empirical rate
real.swed_rate <- swed.data.3090$Deaths/ swed.data.3090$Exposure

# Generating the hazard rate with optimized parameters
mu.hat <- c_hat + a_hat * exp(b_hat*swed.data.3090$Age) 

# Calculating CIs for the optimized rate: delta-method
# Gradients:
g.h.makeham <- function(t){
  
  # Optimized parameters
  a_hat <-  make.mle$par[1]
  b_hat <-  make.mle$par[2]
  c_hat <-  make.mle$par[3]
  
  # Var-Cov matrix from optimization
  V <- solve(-H)
  
  # Gradient for parameter a for h(t) for Makeham:
  g.a <- exp(b_hat * t)
  # Gradient for parameter b for h(t) for Makeham:
  g.b <- a_hat * exp(b_hat * t) * t
  # Gradient for parameter c for h(t) for Makeham:
  g.c <- 1
  
  # Putting it together
  gradient <- c(g.a, g.b, g.c)
  
  # Matrix multiplication:
  se.q <- sqrt( t(gradient) %*% V %*% gradient)
  
  # Getting back what we need
  return(se.q)
  
}

# Create empty vector to fill
se.q.gradient <- rep(NA,length(swed.data.3090$Age))

# Loop through the possible values of t and apply the previous function to calculate
# standard errors for h(t) Makeham based on optimization
for (t in 1:length(swed.data.3090$Age)) {
  
  se.q.gradient[t] <- g.h.makeham(t+29)
  
}

# Creating low and high CI for h(t) hat
CI.make.low <- mu.hat - 1.96 * se.q.gradient
CI.make.high <- mu.hat + 1.96 * se.q.gradient

# Putting everything together
plot(swed.data.3090$Age, real.swed_rate, ylab="mu", main="Fitted Rates") 
lines(swed.data.3090$Age, mu.hat, col=3, lwd=2) 
lines(swed.data.3090$Age, CI.make.low, col=2, lwd=2)
lines(swed.data.3090$Age, CI.make.low, col=2, lwd=2)
legend("topleft", legend=c("Observed", "Fitted", "CI"), 
       col=c(1,3,2), lwd=c(1,1,1), 
       lty=c(1,1,1), pch=c(1,1,1))



#### Exercise 4 ####

##### Importing data #####
blind <- read.table("BlindData.csv", header = T)

glimpse(blind)


##### Q1. Read in the data and produce the Kaplan-Maier estimators of the survival curves for the four groups in the data. ##### 
# Which group has the lowest risk of experiencing blindness?
# Would you accept the proportionality assumptions?

######  Creating Survival object with censoring ######  
blind.surv <- Surv(time=blind$Time, event=blind$Status)

###### Kaplan-Meier fitting for 4 groups (Diabetes + Treatment) ###### 
KM <- survfit(blind.surv ~ Diabetes + Treatment , data=blind)

###### Plotting K-P curves ######
ggsurvplot(KM, data = blind, pval = F,
           legend.title = "Groups:",
           legend.labs = c("D=0, T=0", "D=0, T=1", "D=1, T=0", "D=1, T=1"),
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           palette = wes_palette(n=4, name="Darjeeling1"),
           ggtheme = theme_bw()
           )



##### Q2. Build a data frame that is appropriate for an analysis assuming a Piece-Wise Constant Hazard model. ##### 

# Splitting the data
tauj <- c(6,15,25,40,50,80)
blind2 <- survSplit(blind, cut=tauj, end="Time", start="start",
                   event="Status", episode="interval")

# We need to create the new variable: length spent by each individual in each interval
blind2$y.new <- blind2$Time - blind2$start
blind2$interval <- as.factor(blind2$interval)



##### Q3. Estimate a model for the hazard without any further covariates. ##### 

# Can you reject the hypothesis, that the hazard is constant overall?

# Fitting the Poisson model
fit1.blind <- glm(Status ~ interval, offset = log(y.new), family=poisson, data=blind2)
summary(fit1.blind)

anova(fit1.blind,test="Chisq")


##### Q4. Now build and estimate the following models: ##### 
# mT: Treatment
# mD: Diabetes
# mTD: Treatment + Diabetes
# mTDint: Treatment * Diabetes

###### Running regression models ######
mT <- glm(Status ~ interval + Treatment, offset = log(y.new), family=poisson, data=blind2)
mD <- glm(Status ~ interval + Diabetes, offset = log(y.new), family=poisson, data=blind2)
mTD <- glm(Status ~ interval + Treatment + Diabetes, offset = log(y.new), family=poisson, data=blind2)
mTDint <- glm(Status ~ interval + Treatment * Diabetes, offset = log(y.new), family=poisson, data=blind2)

###### Check regression summary ######
summary(mT)
summary(mD)
summary(mTD)
summary(mTDint)

stargazer(fit1.blind, mT, mD, mTD, mTDint)

###### ANOVA for 4 models ######
anova(fit1.blind, mT, test="Chisq")
anova(fit1.blind, mD, test="Chisq")
anova(mT, mTD, test="Chisq")
anova(mD, mTD, test="Chisq")
anova(mTD, mTDint, test="Chisq")

###### AIC 4 models ######
AIC(fit1.blind, mT, mD, mTD, mTDint)


###### All regression together in 1 table ######
stargazer(fit1.blind, mT, mD, mTD, mTDint)


###### Plotting loghazard interaction model ######
# Define tauj and the hazards for each cohort
tauj2 <- c(0,tauj)
M <- 6
# Baseline: only interval + intercept
haz_cohort1 <- mTDint$coefficients[1:M]
# Adding Treatment, No Diabetes:
haz_cohort2 <- mTDint$coefficients[1:M] + mTDint$coefficients[M+1]
# Adding Diabetes, No Treatment:
haz_cohort3 <- mTDint$coefficients[1:M] + mTDint$coefficients[M+2]
# Adding Diabetes, Treatment & the interaction:
haz_cohort4 <- mTDint$coefficients[1:M] + mTDint$coefficients[M+1] + mTDint$coefficients[M+2] + mTDint$coefficients[M+3]


# Plot the log-hazard over age for each of the three cohorts
plot(1, 1, t="n", ylim = range(haz_cohort1, haz_cohort2, haz_cohort3, haz_cohort4), xlim = range(tauj2), ylab = "log[h(t)]", xlab = "t")
for (i in 1:length(haz_cohort1)) {
  segments(x0=tauj2[i], y0=haz_cohort1[i], x1=tauj2[i+1], y1=haz_cohort1[i], lwd = 2, col = 1, lty = 1)
  segments(x0=tauj2[i], y0=haz_cohort2[i], x1=tauj2[i+1], y1=haz_cohort2[i], lwd = 2, col = 2, lty = 2)
  segments(x0=tauj2[i], y0=haz_cohort3[i], x1=tauj2[i+1], y1=haz_cohort3[i], lwd = 2, col = 3, lty = 2)
  segments(x0=tauj2[i], y0=haz_cohort4[i], x1=tauj2[i+1], y1=haz_cohort4[i], lwd = 2, col = 4, lty = 2)
}
abline(v=6, lty=3) # lines to mark cut points 
abline(v=15, lty=3)
abline(v=25, lty=3)
abline(v=40, lty=3)
abline(v=50, lty=3)
legend("bottomright", legend = c("Baseline", "Baseline + Treatment", 
                              "Baseline + Diabetes", "Baseline + T*D"), col = 1:4, lty = c(1,2,2,2), lwd = 2) 



##### Q5. Estimate a Cox PH model using the selected set of covariates from the previous task ##### 

###### Adding 0s as starting point for Cox ###### 
blind2$entry <- 0

###### Survival object ###### 
surv.dat <- Surv(time = blind2$entry,
                 time2 = blind2$Time,
                 event = blind2$Status)

###### Cox modeling ###### 
cox.blind <- coxph(surv.dat ~ Treatment * Diabetes, data=blind2)


###### Comparing PC model with interaction vs Cox model ######
summary(mTDint)
summary(cox.blind)

AIC(mTDint, cox.blind)


