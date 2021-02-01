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
