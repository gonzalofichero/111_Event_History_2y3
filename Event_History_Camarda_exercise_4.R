#### Importing packages ####
library(tidyverse)
library(wesanderson)

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
# The base category is NO malignant tumor
m.cancer <- glm(cbind(yes,no) ~ age + Malignant + center - 1, family=binomial(link=logit), data=cancer)
summary(m.cancer)

# All p-values below 5%. GOOD



#### 3. Select the most suitable model. Why do you select it? ####

# Let's run model only with age, then only with malignancy and the compared these 2 with the complete model
# Using AIC to compare between the models: the lower the better

m.cancer.age <- glm(cbind(yes,no) ~ age - 1, family=binomial(link=logit), data=cancer)
m.cancer.mal <- glm(cbind(yes,no) ~ Malignant - 1, family=binomial(link=logit), data=cancer)
m.cancer.center <- glm(cbind(yes,no) ~ center - 1, family=binomial(link=logit), data=cancer)


AIC(m.cancer.age, m.cancer.mal, m.cancer.center, m.cancer)
# The complete model performs better





