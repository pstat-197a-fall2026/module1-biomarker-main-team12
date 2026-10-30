# analysis-main.R
# Group 12: module 1, proteomic biomarkers of ASD
#
# Primary script for the report (results/report.qmd). It pools the draft
# scripts in scripts/drafts/ (tasks1and2.R, task3.R, task4.R).
#
# How to run: open the project in RStudio so the working directory is the
# repo root (the folder that contains data/ and results/). Run
# scripts/preprocessing.R first if data/biomarker-clean.RData does not exist.
# The in-class code is in scripts/preliminary-analysis.R and is not repeated here.
#
# Each section sets its own seed, so sections can be run on their own.

library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(glmnet)
library(pROC)
library(rsample)


# ============================================================
# TASK 1 AND 2: log transformation and outliers (uses raw data)
# ============================================================

# read in the raw data with the same column naming as preprocessing.R
var_names <- read_csv("data/biomarker-raw.csv",
                      col_names = FALSE, n_max = 2,
                      col_select = -(1:2), show_col_types = FALSE) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, abbreviation = V2) %>%
  na.omit()

biomarker_raw <- read_csv("data/biomarker-raw.csv",
                          skip = 2, col_select = -2L,
                          col_names = c("group", "empty",
                                        pull(var_names, abbreviation), "ados"),
                          na = c("-", ""), show_col_types = FALSE) %>%
  filter(!is.na(group)) %>%
  select(group, ados, everything())

proteins <- setdiff(names(biomarker_raw), c("group", "ados"))

# ------------------------------------------------------------
# Task 1: why log-transform the protein levels?
# ------------------------------------------------------------
# Compare raw and log10 distributions for four example proteins
# (the first four columns: CEBPB, CHIP, NSE, PIAS4).
sample_proteins <- head(proteins, 4)

sample_distributions <- biomarker_raw %>%
  select(all_of(sample_proteins)) %>%
  pivot_longer(everything(), names_to = "protein", values_to = "raw_level") %>%
  mutate(log_level = log10(raw_level)) %>%
  pivot_longer(c(raw_level, log_level), names_to = "version", values_to = "level") %>%
  mutate(version = recode(version, raw_level = "Raw", log_level = "Log10"))

p_raw_log <- ggplot(sample_distributions, aes(x = level)) +
  geom_histogram(bins = 25, fill = "steelblue", color = "white") +
  facet_wrap(~ protein + version, scales = "free_x", ncol = 4) +
  labs(
    title = "Protein distributions before and after log10 transformation",
    x = "Protein measurement (scale shown in panel)",
    y = "Number of subjects"
  ) +
  theme_minimal()
print(p_raw_log)
ggsave("results/task1-distributions.png", plot = p_raw_log,
       width = 10, height = 6, dpi = 300)

# The raw values are right-skewed: most measurements are low and a few are
# extremely high, especially for CHIP and NSE. Log10 makes the distributions
# more symmetric, reduces the pull of the extreme values, and helps stabilize
# the variability, which suits the t-tests and models used later.

# ------------------------------------------------------------
# Task 2: effect of outlier trimming
# ------------------------------------------------------------
# Log-transform and standardize WITHOUT trimming, and keep a subject ID.
biomarker_untrimmed <- biomarker_raw %>%
  mutate(subject_id = row_number()) %>%
  mutate(across(all_of(proteins),
                ~ as.numeric(scale(log10(.x)))))

# count outlying values (|z| > 3) for each subject
subject_outliers <- biomarker_untrimmed %>%
  pivot_longer(cols = all_of(proteins),
               names_to = "protein",
               values_to = "z_score") %>%
  group_by(subject_id, group) %>%
  summarise(n_outliers = sum(abs(z_score) > 3, na.rm = TRUE),
            .groups = "drop")

# subjects with the most outlying measurements
top_subjects <- subject_outliers %>%
  arrange(desc(n_outliers)) %>%
  slice_head(n = 10)
print(top_subjects)

# subjects in the highest 5% of outlier counts
cutoff <- quantile(subject_outliers$n_outliers, 0.95)

high_outlier_subjects <- subject_outliers %>%
  filter(n_outliers >= cutoff) %>%
  arrange(desc(n_outliers))
print(high_outlier_subjects)

# compare outlier counts between groups
group_summary <- subject_outliers %>%
  group_by(group) %>%
  summarise(
    n_subjects = n(),
    mean_outliers = mean(n_outliers),
    median_outliers = median(n_outliers),
    max_outliers = max(n_outliers),
    n_high_outlier_subjects = sum(n_outliers >= cutoff),
    .groups = "drop"
  )
print(group_summary)

p_outliers <- ggplot(subject_outliers,
                     aes(x = group, y = n_outliers, color = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  geom_hline(yintercept = cutoff, linetype = "dashed") +
  labs(
    title = "Outlying protein measurements per subject",
    x = "Group",
    y = "Number of outlying protein measurements"
  ) +
  theme_minimal() +
  guides(color = "none")
print(p_outliers)
ggsave("results/task2-outliers-by-group.png", plot = p_outliers,
       width = 6, height = 4, dpi = 300)

p_top_subjects <- ggplot(high_outlier_subjects,
                         aes(x = reorder(factor(subject_id), n_outliers),
                             y = n_outliers,
                             fill = group)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Subjects with the highest outlier counts",
    x = "Subject ID",
    y = "Number of outlying protein measurements"
  ) +
  theme_minimal()
print(p_top_subjects)
ggsave("results/task2-top-subjects.png", plot = p_top_subjects,
       width = 6, height = 4, dpi = 300)

# Outlying values are concentrated in a few subjects (154, 108, 9, 121, 52, 77
# each have 114 to 157), not spread evenly. TD has the higher mean (17.6 vs 13.2)
# but the medians are the same (8.5). Five of the eight highest-count subjects
# are TD and three are ASD, so there is no clear evidence that outliers are more
# frequent in one group. Trimming limits the influence of these subjects.


# ============================================================
# TASK 3: methodological variations (uses the cleaned data)
# ============================================================
load("data/biomarker-clean.RData")

# same seed as the in-class analysis
set.seed(101422)

# helper for the t-tests
test_fn <- function(.df){
  t_test(.df,
         formula = level ~ group,
         order = c("ASD", "TD"),
         alternative = "two-sided",
         var.equal = F)
}
class_metrics <- metric_set(sensitivity,
                            specificity,
                            accuracy,
                            roc_auc)

# In-class results for comparison (5-protein panel DERM, RELT, MRC2, IgD,
# Cadherin-5; test sensitivity 0.812, specificity 0.733, accuracy 0.774,
# AUROC 0.883).
#
# Note: every modification below uses ONE random split with about 31 test
# subjects, so results are noisy (one subject moves sensitivity by about 0.06).
# The splits are also not the same across modifications (seed is only set once
# and the split is not stratified), so differences between them are within noise.

# ------------------------------------------------------------
# Modification 1: carry out the whole selection on a training partition
# ------------------------------------------------------------
# 80% training, 20% testing split made at the very beginning
biomarker_split_m1 <- initial_split(biomarker_clean, prop = 0.8)
train_data_m1 <- training(biomarker_split_m1)
test_data_m1  <- testing(biomarker_split_m1)

# multiple testing (training data only)
ttests_out_m1 <- train_data_m1 %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = "protein", values_to = "level") %>%
  nest(data = c(level, group)) %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(), hm = log(m) + 1/(2*m) - digamma(1), rank = row_number(),
         p.adj = m*hm*p_value/rank)

proteins_s1_m1 <- ttests_out_m1 %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)
proteins_s1_m1

# random forest (training data only)
predictors_m1 <- train_data_m1 %>%
  select(-c(group, ados))
response_m1 <- train_data_m1 %>%
  pull(group) %>%
  factor()

rf_out_m1 <- randomForest(x = predictors_m1, y = response_m1,
                          ntree = 1000, importance = T)
proteins_s2_m1 <- rf_out_m1$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out_m1$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)
proteins_s2_m1

# logistic regression on the intersection
proteins_sstar_m1 <- intersect(proteins_s1_m1, proteins_s2_m1)
proteins_sstar_m1   # DERM, IgD, TSP4, MRC2

train_sstar_m1 <- train_data_m1 %>%
  select(group, any_of(proteins_sstar_m1)) %>%
  mutate(class = (group == "ASD")) %>%
  select(-group)
test_sstar_m1 <- test_data_m1 %>%
  select(group, any_of(proteins_sstar_m1)) %>%
  mutate(class = (group == "ASD")) %>%
  select(-group)

fit_m1 <- glm(class ~ ., data = train_sstar_m1, family = "binomial")

test_sstar_m1 %>%
  add_predictions(fit_m1, type = "response") %>%
  mutate(truth_f = factor(class, levels = c(FALSE, TRUE)),
         est_f = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
  class_metrics(truth = truth_f,
                estimate = est_f,
                pred,
                event_level = "second") %>%
  print()

# Result: sensitivity 0.688, specificity 0.733, accuracy 0.710, AUROC 0.846.
# Compared with the in-class results, sensitivity, accuracy and AUROC are lower
# and specificity is the same. This is the expected direction: in class the test
# subjects helped choose the proteins, which makes those numbers optimistic.
# The panel shares only some proteins with the in-class one (DERM, IgD, MRC2).

# ------------------------------------------------------------
# Modification 2: a larger number of top proteins per method (K = 15)
# ------------------------------------------------------------
K <- 15

ttests_out_m2 <- biomarker_clean %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = "protein", values_to = "level") %>%
  nest(data = c(level, group)) %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(), hm = log(m) + 1/(2*m) - digamma(1), rank = row_number(),
         p.adj = m*hm*p_value/rank)

proteins_s1_m2 <- ttests_out_m2 %>% slice_min(p.adj, n = K) %>% pull(protein)
proteins_s1_m2

predictors_m2 <- biomarker_clean %>% select(-c(group, ados))
response_m2 <- biomarker_clean %>% pull(group) %>% factor()

rf_out_m2 <- randomForest(x = predictors_m2, y = response_m2,
                          ntree = 1000, importance = T)
proteins_s2_m2 <- rf_out_m2$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out_m2$importance)) %>%
  slice_max(MeanDecreaseGini, n = K) %>%
  pull(protein)
proteins_s2_m2

proteins_sstar_m2 <- intersect(proteins_s1_m2, proteins_s2_m2)
proteins_sstar_m2   # 8 proteins

biomarker_sstar_m2 <- biomarker_clean %>%
  select(group, any_of(proteins_sstar_m2)) %>%
  mutate(class = (group == "ASD")) %>%
  select(-group)

biomarker_split_m2 <- biomarker_sstar_m2 %>%
  initial_split(prop = 0.8)
fit_m2 <- glm(class ~ ., data = training(biomarker_split_m2), family = "binomial")

testing(biomarker_split_m2) %>%
  add_predictions(fit_m2, type = "response") %>%
  mutate(truth_f = factor(class, levels = c(FALSE, TRUE)),
         est_f = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
  class_metrics(truth = truth_f, estimate = est_f, pred, event_level = "second") %>%
  print()

# Result: sensitivity 0.647, specificity 0.786, accuracy 0.710, AUROC 0.857.
# With K = 15 the two lists overlap in 8 proteins instead of a handful, but test
# performance did not improve on this split (AUROC and accuracy are about the
# same as Modification 1 and below the in-class results). DERM and IgD are in
# every panel so far.

# ------------------------------------------------------------
# Modification 3: fuzzy intersection (union of the two top 10 lists)
# ------------------------------------------------------------
ttests_out_m3 <- biomarker_clean %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = "protein", values_to = "level") %>%
  nest(data = c(level, group)) %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(), hm = log(m) + 1/(2*m) - digamma(1), rank = row_number(),
         p.adj = m*hm*p_value/rank)

proteins_s1_m3 <- ttests_out_m3 %>% slice_min(p.adj, n = 10) %>% pull(protein)
proteins_s1_m3

predictors_m3 <- biomarker_clean %>%
  select(-c(group, ados))
response_m3 <- biomarker_clean %>%
  pull(group) %>%
  factor()

rf_out_m3 <- randomForest(x = predictors_m3, y = response_m3,
                          ntree = 1000, importance = T)
proteins_s2_m3 <- rf_out_m3$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out_m3$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)
proteins_s2_m3

# fuzzy intersection used here: the union of the two lists
proteins_sstar_m3 <- union(proteins_s1_m3, proteins_s2_m3)
proteins_sstar_m3   # 16 proteins

biomarker_sstar_m3 <- biomarker_clean %>%
  select(group, any_of(proteins_sstar_m3)) %>%
  mutate(class = (group == "ASD")) %>%
  select(-group)

biomarker_split_m3 <- biomarker_sstar_m3 %>%
  initial_split(prop = 0.8)
fit_m3 <- glm(class ~ ., data = training(biomarker_split_m3), family = "binomial")

testing(biomarker_split_m3) %>%
  add_predictions(fit_m3, type = "response") %>%
  mutate(truth_f = factor(class, levels = c(FALSE, TRUE)),
         est_f = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
  class_metrics(truth = truth_f, estimate = est_f, pred, event_level = "second") %>%
  print()

# Result: sensitivity 0.667, specificity 0.789, accuracy 0.742, AUROC 0.781.
# The union gives 16 proteins. AUROC is the lowest of the three modifications,
# and with only about 123 training subjects a 16-protein model likely overfits.


# ============================================================
# TASK 4: a simpler panel with comparable accuracy
# ============================================================
# Benchmark: the in-class panel. Note that the random forest top 10 is not
# perfectly stable across runs, so this list is hard-coded to match the in-class
# results we benchmark against.
baseline_panel <- c("DERM", "RELT", "MRC2", "IgD", "Cadherin-5")

dat <- biomarker_clean %>%
  select(-ados) %>%
  mutate(y = as.numeric(group == "ASD"),
         strat = factor(group)) %>%
  select(-group)

# one split: returns test AUROC, panel size, and panel for each method
run_once <- function(seed) {
  set.seed(seed)
  sp <- initial_split(dat, prop = 0.75, strata = strat)
  train <- training(sp); test <- testing(sp)
  
  x_tr <- as.matrix(select(train, -y, -strat)); y_tr <- train$y
  x_te <- as.matrix(select(test,  -y, -strat)); y_te <- test$y
  
  # LASSO with cross-validation on the training data only
  cv <- cv.glmnet(x_tr, y_tr, family = "binomial",
                  alpha = 1, type.measure = "auc", nfolds = 5)
  
  # LASSO methods use the LASSO's own predictions
  lasso_res <- map_dfr(c("lambda.min", "lambda.1se"), function(s) {
    b <- coef(cv, s = s)
    panel <- setdiff(rownames(b)[as.vector(b != 0)], "(Intercept)")
    pr <- as.vector(predict(cv, newx = x_te, s = s, type = "response"))
    tibble(method = paste0("lasso_", sub("lambda.", "", s)),
           auc = as.numeric(auc(roc(y_te, pr, quiet = TRUE))),
           size = length(panel),
           panel = list(panel))
  })
  
  # baseline: refit a logistic regression on the in-class panel
  fit <- glm(y ~ ., data = train[, c("y", baseline_panel)], family = binomial)
  pr <- predict(fit, test[, baseline_panel], type = "response")
  base_res <- tibble(method = "baseline",
                     auc = as.numeric(auc(roc(y_te, pr, quiet = TRUE))),
                     size = length(baseline_panel),
                     panel = list(baseline_panel))
  
  bind_rows(lasso_res, base_res) %>% mutate(seed = seed)
}

# repeat over 50 random splits (takes a minute or two)
res <- map_dfr(1:50, run_once)

# summary table
print(res %>%
        group_by(method) %>%
        summarize(mean_auc = mean(auc), sd_auc = sd(auc), mean_size = mean(size)))

# most frequently selected proteins (1se LASSO)
print(res %>%
        filter(method == "lasso_1se") %>%
        unnest(panel) %>%
        count(panel, sort = TRUE) %>%
        head(15))

# Neither LASSO panel was simpler and better: baseline mean AUROC 0.806 (5
# proteins), LASSO 1se 0.784 (about 15 proteins), LASSO min 0.819 (about 40).
# DERM, IgD and MAPK2 were chosen most often by the 1se LASSO.

# fixed 3-protein panel vs baseline on the same 50 splits
# (this panel was suggested by the LASSO counts above, which came from these
# same splits, so the comparison favors it somewhat)
alt_panel <- c("DERM", "IgD", "MAPK2")

res2 <- map_dfr(1:50, function(seed) {
  set.seed(seed)
  sp <- initial_split(dat, prop = 0.75, strata = strat)
  train <- training(sp); test <- testing(sp)
  map_dfr(list(baseline = baseline_panel, alt = alt_panel), function(p) {
    fit <- glm(y ~ ., data = train[, c("y", p)], family = binomial)
    pr <- predict(fit, test[, p], type = "response")
    tibble(auc = as.numeric(auc(roc(test$y, pr, quiet = TRUE))))
  }, .id = "method") %>% mutate(seed = seed)
})

print(res2 %>%
        group_by(method) %>%
        summarize(mean_auc = mean(auc), sd_auc = sd(auc)))

# paired difference on the same splits
print(res2 %>%
        pivot_wider(names_from = method, values_from = auc) %>%
        summarize(mean_diff = mean(alt - baseline),
                  sd_diff = sd(alt - baseline)))

# Result: the 3-protein panel has mean AUROC 0.848 (sd 0.046) vs 0.806 for the
# in-class panel; the mean paired difference is +0.042 (sd 0.041). It is a
# simpler panel with comparable accuracy, but the gain is less than one sd of
# the baseline (0.062), so we do not call it an improvement.

# boxplots saved for the report
p4a <- ggplot(res, aes(method, auc)) +
  geom_boxplot() +
  scale_x_discrete(labels = c(baseline  = "In-class\n(5 proteins)",
                              lasso_1se = "LASSO 1se",
                              lasso_min = "LASSO min")) +
  labs(y = "Test AUROC", x = NULL,
       title = "Test AUROC across 50 random splits")
ggsave("results/task4-boxplot1.png", plot = p4a, width = 6, height = 4, dpi = 300)

p4b <- ggplot(res2, aes(method, auc)) +
  geom_boxplot() +
  scale_x_discrete(labels = c(baseline = "In-class\n(5 proteins)",
                              alt      = "DERM + IgD + MAPK2\n(3 proteins)")) +
  labs(y = "Test AUROC", x = NULL,
       title = "In-class panel vs. 3-protein panel (50 splits)")
ggsave("results/task4-boxplot2.png", plot = p4b, width = 6, height = 4, dpi = 300)