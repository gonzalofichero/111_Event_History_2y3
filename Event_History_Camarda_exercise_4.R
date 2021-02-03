#### Importing packages ####
library(tidyverse)
library(wesanderson)
library(ggcats)

#### 1. Load the data and quickly present a possible research questions with the help of descriptive statistics ####
cancer <- read.table("CancerBishopEtAl1975.txt", header=TRUE, sep=",")

glimpse(cancer)


# Create Survivor % rate for each group
# Colors by type of Center, split into 2 by type of tumor
# Calculate Avg Rate of Survival between type of tumor and add dashed line to facet_wrap
# Use Wes Anderson's palette for filling the bars

cancer %>% 
  mutate(survivor = yes/(yes+no)) %>%
  group_by(Malignant) %>% 
  summarise(MN = mean(survivor)) -> mean_mal

cancer %>% 
  mutate(survivor = yes/(yes+no)) %>% 
  ggplot(aes(x=age, y = survivor, fill = center)) + geom_bar(position="dodge", stat="identity") +
  geom_hline(data = mean_mal, aes(yintercept = MN), linetype="dashed", color = "black", size=1.1) +
  facet_wrap(~ Malignant) +
  scale_fill_manual(values=wes_palette(n=3, name="Darjeeling1"))


# We have some mixed results:
# a1) When Tumor is NO malignant, age and survivor are negatively correlated
# a2) When Tumor is YES malignant, we have mixed results: I believe there is another dimension we are missing
# b) Centers have different survivor rates. Tokyo is always better...
# c) Type of Tumor is definitely something to take into account (dashed line)

# We'll try using all dimensions to model the survivorship: YES vs NO => logistic


#### 2. Model these data to answer the question (do not use any interaction) ####

# I'm taking out the intercept so reading the betas is easier

##### Models with 1 variable #####
m.cancer.age <- glm(cbind(yes,no) ~ age - 1, family=binomial(link=logit), data=cancer)
m.cancer.mal <- glm(cbind(yes,no) ~ Malignant - 1, family=binomial(link=logit), data=cancer)
m.cancer.center <- glm(cbind(yes,no) ~ center - 1, family=binomial(link=logit), data=cancer)


summary(m.cancer.age)
summary(m.cancer.mal)
summary(m.cancer.center)
# Everything looks significative in a 1by1 basis


##### Models with 2 variables #####
m.cancer.age.mal <- glm(cbind(yes,no) ~ age +  Malignant - 1, family=binomial(link=logit), data=cancer)
m.cancer.mal.center <- glm(cbind(yes,no) ~ Malignant + center - 1, family=binomial(link=logit), data=cancer)
m.cancer.center.age <- glm(cbind(yes,no) ~ center + age - 1, family=binomial(link=logit), data=cancer)

summary(m.cancer.age.mal)
summary(m.cancer.mal.center)
summary(m.cancer.center.age)

# Combination of Center + Age makes ages look significance

##### Full model #####
# The base category is NO malignant tumor
m.cancer.full <- glm(cbind(yes,no) ~ age + Malignant + center - 1, family=binomial(link=logit), data=cancer)
summary(m.cancer.full)

# All p-values below 5%, except for Glamorgan. GOOD


#### 3. Select the most suitable model. Why do you select it? ####
#### 4. Present a table with the Akaike Information Criterion (AIC) for all estimated models. ####

# Using AIC to compare between the models: the lower the better

AIC(m.cancer.age, m.cancer.mal, m.cancer.center, 
    m.cancer.age.mal, m.cancer.mal.center, m.cancer.center.age,
    m.cancer.full)
# The full model and the one with type of tumor and center performs better
# By parsimony, selecting the model with less parameters: m.cancer.mal.center



#### 5. Plot the probability of dying from the selected model by the significant covariates. Add 95% confidence intervals ####


# Vector with all possible values:
cases <- data.frame(Malignant=c("no", "yes","no", "yes","no", "yes"), 
                    center=c("Boston", "Boston", "Glamorgan", "Glamorgan", "Tokyo", "Tokyo"))


# Saving predictions inside vector
p.cancer <- predict(m.cancer.mal.center, cases, type="response", se.fit = TRUE)

# All together now...
prediction <- cbind(cases, p.cancer$fit, p.cancer$se.fit)
names(prediction)[3] <- "model.prob"
names(prediction)[4] <- "prob.se"

# Calculating lower and upper limit of prediction
prediction %>% 
  mutate(low_pred = model.prob - 1.96 * prob.se,
         high_pred = model.prob + 1.96 * prob.se) -> prediction


# Putting everything together inside the same plot
cancer %>% 
  mutate(survivor = yes/(yes+no)) %>% 
  select(center, Malignant, survivor) %>%
  group_by(center, Malignant) %>% 
  summarise(avg_survivor = mean(survivor)) %>% 
  left_join(prediction, by=c("Malignant", "center")) %>% 
  ggplot(aes(x=center, y = avg_survivor, fill = Malignant)) + geom_bar(position="dodge", stat="identity") +
  #geom_point(aes(x=center, y = model.prob, color = Malignant), position = position_dodge(width = .9)) +
  geom_cat(aes(x=center, y = model.prob), cat = "grumpy", size = 1.5, position = position_dodge(width = .9)) +
  geom_cat(aes(x=center, y = low_pred), cat = "pusheen", size = 0.75, position = position_dodge(width = .9)) +
  geom_cat(aes(x=center, y = high_pred), cat = "pusheen", size = 0.75, position = position_dodge(width = .9)) +
  #scale_color_manual(values=wes_palette(n=2, name="IsleofDogs1")) +
  scale_fill_manual(values=wes_palette(n=3, name="Darjeeling1"))





