library(tidyverse)
library(glmnet)
library(pROC)
library(rsample)

# load data
load('data/biomarker-clean.RData')

dat <- biomarker_clean %>%
  select(-ados) %>%
  mutate(y = as.numeric(group == 'ASD'),
         strat = factor(group)) %>%
  select(-group)

# baseline panel from in-class analysis
baseline_panel <- c('DERM', 'RELT', 'MRC2', 'IgD', 'Cadherin-5')

# one split: returns test AUC, panel size, and panel
run_once <- function(seed) {
  set.seed(seed)
  sp <- initial_split(dat, prop = 0.75, strata = strat)
  train <- training(sp); test <- testing(sp)
  
  x_tr <- as.matrix(select(train, -y, -strat)); y_tr <- train$y
  x_te <- as.matrix(select(test,  -y, -strat)); y_te <- test$y
  
  # LASSO with CV on training data only
  cv <- cv.glmnet(x_tr, y_tr, family = 'binomial',
                  alpha = 1, type.measure = 'auc', nfolds = 5)
  
  # LASSO methods: use the LASSO's own predictions
  lasso_res <- map_dfr(c('lambda.min', 'lambda.1se'), function(s) {
    b <- coef(cv, s = s)
    panel <- setdiff(rownames(b)[as.vector(b != 0)], '(Intercept)')
    pr <- as.vector(predict(cv, newx = x_te, s = s, type = 'response'))
    tibble(method = paste0('lasso_', sub('lambda.', '', s)),
           auc = as.numeric(auc(roc(y_te, pr, quiet = TRUE))),
           size = length(panel),
           panel = list(panel))
  })
  
  # Baseline: refit glm on the in-class panel using training data
  fit <- glm(y ~ ., data = train[, c('y', baseline_panel)], family = binomial)
  pr <- predict(fit, test[, baseline_panel], type = 'response')
  base_res <- tibble(method = 'baseline',
                     auc = as.numeric(auc(roc(y_te, pr, quiet = TRUE))),
                     size = length(baseline_panel),
                     panel = list(baseline_panel))
  
  bind_rows(lasso_res, base_res) %>% mutate(seed = seed)
}

# repeat over 50 splits (takes a minute or two)
res <- map_dfr(1:50, run_once)

# summary table
res %>%
  group_by(method) %>%
  summarize(mean_auc  = mean(auc),
            sd_auc    = sd(auc),
            mean_size = mean(size))

# most frequently selected proteins (1se LASSO)
res %>%
  filter(method == 'lasso_1se') %>%
  unnest(panel) %>%
  count(panel, sort = TRUE) %>%
  head(15)

# boxplot
ggplot(res, aes(method, auc)) +
  geom_boxplot() +
  labs(y = 'Test AUROC', x = NULL)

print(res %>%
        group_by(method) %>%
        summarize(mean_auc = mean(auc), sd_auc = sd(auc), mean_size = mean(size)))


print(res %>%
        filter(method == 'lasso_1se') %>%
        unnest(panel) %>%
        count(panel, sort = TRUE) %>%
        head(15))



#fixed 3-protein panel vs baseline (same 50 splits) 
alt_panel <- c('DERM', 'IgD', 'MAPK2')

res2 <- map_dfr(1:50, function(seed) {
  set.seed(seed)
  sp <- initial_split(dat, prop = 0.75, strata = strat)
  train <- training(sp); test <- testing(sp)
  map_dfr(list(baseline = baseline_panel, alt = alt_panel), function(p) {
    fit <- glm(y ~ ., data = train[, c('y', p)], family = binomial)
    pr <- predict(fit, test[, p], type = 'response')
    tibble(auc = as.numeric(auc(roc(test$y, pr, quiet = TRUE))))
  }, .id = 'method') %>% mutate(seed = seed)
})

print(res2 %>%
        group_by(method) %>%
        summarize(mean_auc = mean(auc), sd_auc = sd(auc)))

# paired difference on the same splits
print(res2 %>%
        pivot_wider(names_from = method, values_from = auc) %>%
        summarize(mean_diff = mean(alt - baseline),
                  sd_diff = sd(alt - baseline)))

#Load boxplot images and save to results for final report. 
p <- ggplot(res, aes(method, auc)) +
  geom_boxplot() +
  scale_x_discrete(labels = c(baseline  = 'In-class\n(5 proteins)',
                              lasso_1se = 'LASSO 1se',
                              lasso_min = 'LASSO min')) +
  labs(y = 'Test AUROC', x = NULL,
       title = 'Test AUROC across 50 random splits')

ggsave('results/task4-boxplot1.png', plot = p, width = 6, height = 4, dpi = 300)

p2 <- ggplot(res2, aes(method, auc)) +
  geom_boxplot() +
  scale_x_discrete(labels = c(baseline = 'In-class\n(5 proteins)',
                              alt      = 'DERM + IgD + MAPK2\n(3 proteins)')) +
  labs(y = 'Test AUROC', x = NULL,
       title = 'In-class panel vs. 3-protein panel (50 splits)')

ggsave('results/task4-boxplot2.png', plot = p2, width = 6, height = 4, dpi = 300)


