#### Importing packages ####
library(tidyverse)


#### Contraceptive Use ####

# Loading package
fs <- read.table("fsMod.txt",header=TRUE)

glimpse(fs)


# Change reference level of desire var
fs$desire <- relevel(fs$desire, "yes")


##### Modeling #####

###### Solving problem as a t-test ######
ferN <- glm(cbind(contra,non.contra)~1, family=binomial(link=logit), data=fs)
summary(ferN)
# Ok, they are not the same. There is a statistically significant difference

# The estimate of the intercept is actually the logit of the sample: 
# sum(fs$contra)/sum(fs$total)
# log(p/(1-p))


# A 95% confidence interval for the overall probability of using contraception equal to:
exp(confint(ferN)) / (1+exp(confint(ferN)))


###### Solving problem as 1 var regression ######
fer1 <- glm(cbind(contra,non.contra)~desire, family=binomial, data=fs)
summary(fer1)



###### Solving problem as 2 var regression ######
fer2 <- glm(cbind(contra,non.contra)~educ+desire, family=binomial, data=fs)
summary(fer2)



###### Solving problem as 3 var regression ######
fer3 <- glm(cbind(contra,non.contra)~educ+desire+age, family=binomial, data=fs)
summary(fer3)

# Now education seems to be working as explanatory variable, but...
anova(fer3, test="Chisq")

# If we run the anova, Education is again no significative


###### Solving problem with interaction of 3 vars ######

# Vector with the right-hand-side formulas
rhs <- c("1","age", "educ", "desire",
         "age + educ", "age + desire", "educ + desire",
         "age * educ", "age * desire", "educ * desire",
         "age + educ + desire",
         "age * educ + desire", "age * desire + educ",
         "age + educ * desire", "age * (educ + desire)",
         "educ * (age + desire)", "(age + educ) * desire",
         "age*educ*desire - age:educ:desire",
         "age*educ*desire")

# Then we include the response for each value in rhs
formulas <- paste("cbind(contra,non.contra)~", rhs)

# We enforce the elements in the previous vectors to be encoded as formula
formulae <- lapply(formulas, formula)

# The object formulae is a list containing all possible model formulation and therefore we (l)-apply glm
# to each element of it:
models <- lapply(formulae, glm, family=binomial, data=fs)

# Finally we extract a table with deviance and degree of freedom:
anovatab <- data.frame(model = rhs,
                       deviance = round(unlist(lapply(models,deviance)),2),
                       df = unlist(lapply(models,df.residual)) )

anovatab



#### Clinical Trial on Polyposis ####

polyps <- read.table("polyps.txt", sep="", header=T)

glimpse(polyps)


###### Solving problem with Null model ######
m0 <- glm(number ~ 1, data = polyps, family = poisson(link=log))

# We estimate a model log-linear model with only a mean, we thus estimated the log of the mean number of polyps:
m0$coef
log(sum(polyps$number)/nrow(polyps))


###### Solving problem with 1 var ######
mA <- glm(number ~ age, data = polyps, family = poisson(link=log))
summary(mA)

# An additional year in age will result to a decrease in the number of polyps equal to about 3%
# $ 1 - \exp{`r mA$mA$coefficients[2]`} = 0.031 $

# Calculation with some fitted values for age 20 and 21
pA <- predict(mA, data.frame(age=c(20, 21)), type="response")
pA
# The estimated coefficient inform about the (constant) difference in the logarithm of the estimated polyps


###### Solving problem with complete model ######
mAT <- glm(number ~ age + treatment + age*treatment, data = polyps, family = poisson(link=log))
summary(mAT)

anova(mAT)
# No interaction is needed

mAT <- glm(number ~ age + treatment, data = polyps, family = poisson(link=log))
summary(mAT)

# A final plot could help to check the outcomes of the model, and we do in a log-scale to appreciate the
# linearity and the parallelism between the placebo and treated group:
ages <- seq(min(polyps$age), max(polyps$age), length=100)
mydf <- expand.grid(age=ages, treatment=c(0,1))
y.hat <- matrix(predict(mAT, newdata = mydf),100)
plot(polyps$age, log(polyps$number), col=polyps$treatment+1, xlab="age", ylab="number of polyps, log-scale")
matlines(ages, y.hat, col=1:2, lty=1)
legend("topright", c("Placebo", "Treatment"), col=1:2, lty=1, pch=1)



###### Solving problem with aggregate data ######

mytab <- table(polyps$age, polyps$treatment)
mydf <- expand.grid(age=as.numeric(dimnames(mytab)[[1]]),
                    treat=as.numeric(dimnames(mytab)[[2]]))
mydf$e <- c(mytab)

sumnum <- tapply(polyps$number, list(polyps$age, polyps$treatment), sum)
mydf$y <- as.vector(sumnum)
mydf <- mydf[mydf$e!=0, ]

# The dataframe mydf contains the same information of the original dataset and fitting the same model will lead to equal outcomes
# OFFSET!! is necessary when using aggregate data
mATagg <- glm(y~age+treat, offset=log(e), family=poisson(), data=mydf)
mATagg$coef
# vs.
mAT$coef




