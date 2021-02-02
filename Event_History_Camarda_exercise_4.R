#### Importing packages ####
library(tidyverse)


#### 1. Load the data and quickly present a possible research questions with the help of descriptive statistics ####
cancer <- read.table("cancer.txt",header=TRUE)

glimpse(cancer)

# Create Survivor % rate for each group
# Pivot_wider to have 2 columns for malignant dimension
cancer %>% 
  mutate(survivor = yes/(yes+no)) %>%
  select(-c(yes,no)) %>% 
  pivot_wider(names_from=malignant, names_prefix="malignant_", values_from=survivor) -> cancer2

#Create a Matrix which will help in creating the plot
value_matrix = matrix(, nrow = 2, ncol = 3)

#An empty matrix is a necessary requirement prior to copying data
value_matrix[1,] = cancer2$malignant_no 
value_matrix[2,] = cancer2$malignant_yes

# Note that the "beside" argument has to be kept "TRUE" in order to place the bars side by side
barplot(value_matrix, names.arg = cancer2$age, beside = TRUE, col = c("peachpuff", "skyblue"), legend.text = c("No malignant", "Malignant"))


# From plot we see that Survivor seems to be negatively correlated to both Age and the "malignity" of the tumor
# We'll try using both dimensions to model the survivorship: YES vs NO => logistic


#### 2. Model these data to answer the question (do not use any interaction) ####

# I'm taking out the intercept so reading the betas is easier
# The base category is NO malignant tumor
m.cancer <- glm(cbind(yes,no) ~ age + malignant - 1, family=binomial(link=logit), data=cancer)
summary(m.cancer)

# All p-values below 5%. GOOD



#### 3. Select the most suitable model. Why do you select it? ####

# Let's run model only with age, then only with malignancy and the compared these 2 with the complete model
# Using AIC to compare between the models: the lower the better

m.cancer.age <- glm(cbind(yes,no) ~ age - 1, family=binomial(link=logit), data=cancer)
m.cancer.mal <- glm(cbind(yes,no) ~ malignant - 1, family=binomial(link=logit), data=cancer)


AIC(m.cancer.age, m.cancer.mal, m.cancer)
# The complete model performs better





