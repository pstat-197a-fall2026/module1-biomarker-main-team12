# Task 3: Experimental analysis with several modifications


# Loading data and libraries
library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(glmnet)

load('data/biomarker-clean.RData') 

set.seed(101422) #same seed as the in-class analysis

# Helper methods for t-test analysis
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)


# Modification 1: Creating a Training/Testing Split
# Here, we are doing an 80% training, 20% testing split of the data
biomarker_split_m1 <- initial_split(biomarker_clean, prop = 0.8)
train_data_m1 <- training(biomarker_split_m1)
test_data_m1  <- testing(biomarker_split_m1)

# Multiple Testing
ttests_out_m1 <- train_data_m1 %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = 'protein', values_to = 'level') %>%
  nest(data = c(level, group)) %>% 
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(), hm = log(m) + 1/(2*m) - digamma(1), rank = row_number(), p.adj = m*hm*p_value/rank)

proteins_s1_m1 <- ttests_out_m1 %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

proteins_s1_m1

# Random Forest
predictors_m1 <- train_data_m1 %>%
  select(-c(group, ados))
response_m1 <- train_data_m1 %>%
  pull(group) %>%
  factor()

rf_out_m1 <- randomForest(x = predictors_m1, y = response_m1, ntree = 1000, importance = T)
proteins_s2_m1 <- rf_out_m1$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out_m1$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

proteins_s2_m1

# Logistic Regression
proteins_sstar_m1 <- intersect(proteins_s1_m1, proteins_s2_m1)

train_sstar_m1 <- train_data_m1 %>%
  select(group, any_of(proteins_sstar_m1)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)
test_sstar_m1  <- test_data_m1  %>%
  select(group, any_of(proteins_sstar_m1)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

fit_m1 <- glm(class ~ ., data = train_sstar_m1, family = 'binomial')

proteins_sstar_m1

# Metrics
test_sstar_m1 %>%
  add_predictions(fit_m1, type = 'response') %>%
  mutate(truth_f = factor(class, levels = c(FALSE, TRUE)),
         est_f = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
  class_metrics(truth = truth_f,
                estimate = est_f,
                pred,
                event_level = 'second') %>%
  print()

# Comparing the results to that of the in-class analysis, with this seed we notice a slight decrease in sensitivity of the model, but specificity increased dramatically.
# Additionally, there was a minor improvement to accuracy, but nothing substantial.
# However, we also note an increase of roughly 0.2 to the AUC statistic, meaning it did meaningfully increase the discriminative ability of the model.

# Modification 2: Larger Amount of Top Predictive Proteins
# Here, we will use K = 25

K <- 15

# Multiple Testing
ttests_out_m2 <- biomarker_clean %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = 'protein', values_to = 'level') %>%
  nest(data = c(level, group)) %>% 
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(), hm = log(m) + 1/(2*m) - digamma(1), rank = row_number(), p.adj = m*hm*p_value/rank)

proteins_s1_m2 <- ttests_out_m2 %>% slice_min(p.adj, n = K) %>% pull(protein) #pulling 15

proteins_s1_m2

# Random Forest
predictors_m2 <- biomarker_clean %>% select(-c(group, ados))
response_m2 <- biomarker_clean %>% pull(group) %>% factor()

rf_out_m2 <- randomForest(x = predictors_m2, y = response_m2, ntree = 1000, importance = T)
proteins_s2_m2 <- rf_out_m2$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out_m2$importance)) %>%
  slice_max(MeanDecreaseGini, n = K) %>% #slicing 15
  pull(protein)

proteins_s2_m2

# Logistic Regression
proteins_sstar_m2 <- intersect(proteins_s1_m2, proteins_s2_m2)

biomarker_sstar_m2 <- biomarker_clean %>%
  select(group, any_of(proteins_sstar_m2)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_split_m2 <- biomarker_sstar_m2 %>%
  initial_split(prop = 0.8)
fit_m2 <- glm(class ~ ., data = training(biomarker_split_m2), family = 'binomial')

proteins_sstar_m2

# Metrics
testing(biomarker_split_m2) %>%
  add_predictions(fit_m2, type = 'response') %>%
  mutate(truth_f = factor(class, levels = c(FALSE, TRUE)), est_f = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
  class_metrics(truth = truth_f, estimate = est_f, pred, event_level = 'second') %>%
  print()

# As a result of having K = 15 predictive proteins, the number of intersecting proteins between the random forest and multiple testing increases. 
# 'DERM' and 'IgD' both appear in this intersection and the one for the training/testing split change, but the rest are different for both.
# As for the metrics of accuracy, the AUC and accuracy for this modification were about the same as the first one.
# However, increasing K led to a model with higher sensitivity at the cost of lower specificity, meaning it was better for prediction of ASD cases.

# Modification 3: Fuzzy Intersection

# Multiple Testing
ttests_out_m3 <- biomarker_clean %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = 'protein', values_to = 'level') %>%
  nest(data = c(level, group)) %>% 
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(), hm = log(m) + 1/(2*m) - digamma(1), rank = row_number(), p.adj = m*hm*p_value/rank)

proteins_s1_m3 <- ttests_out_m3 %>% slice_min(p.adj, n = 10) %>% pull(protein)

proteins_s1_m3

# 2. Random Forest (Top 10 proteins)
predictors_m3 <- biomarker_clean %>%
  select(-c(group, ados))
response_m3 <- biomarker_clean %>%
  pull(group) %>%
  factor()

rf_out_m3 <- randomForest(x = predictors_m3, y = response_m3, ntree = 1000, importance = T)
proteins_s2_m3 <- rf_out_m3$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out_m3$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

proteins_s2_m3

# Logistic Regression
# Using fuzzy intersection, i.e. union
proteins_sstar_m3 <- union(proteins_s1_m3, proteins_s2_m3)

biomarker_sstar_m3 <- biomarker_clean %>%
  select(group, any_of(proteins_sstar_m3)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_split_m3 <- biomarker_sstar_m3 %>%
  initial_split(prop = 0.8)
fit_m3 <- glm(class ~ ., data = training(biomarker_split_m3), family = 'binomial')

proteins_sstar_m3

# Metrics
testing(biomarker_split_m3) %>%
  add_predictions(fit_m3, type = 'response') %>%
  mutate(truth_f = factor(class, levels = c(FALSE, TRUE)), est_f = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
  class_metrics(truth = truth_f, estimate = est_f, pred, event_level = 'second') %>%
  print()

# In this case, sensitivity drops considerably while specificity becomes much higher.
# This means the model gets much more reserved about classifying ASD cases, but identifies TD cases with high accuracy.
# Moreover, AUC increases to the highest out of all three modifications, so the model overall performs well.











