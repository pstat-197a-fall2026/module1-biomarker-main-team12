# Tasks 1 and 2: Impact of preprocessing and outliers
# Run from the main project folder (the folder containing 'data/').

library(tidyverse)
library(randomForest)

# read in the data
# use same column naming method as preprocessing.R.
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
# TASK 1: WHY LOG-TRANSFORM PROTEIN LEVELS?
# Question: What do you imagine is the reason for log-transforming the protein levels in biomarker-raw.csv?
# ------------------------------------------------------------
# Compare raw and log10 distributions for four example proteins.
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

# The protein levels were likely log-transformed because the raw data are right-skewed, with most measurements clustered at lower values and a few extremely high values. 
# This is especially noticeable for `CHIP` and `NSE`. 
# Taking the log10 of the measurements reduces the influence of these extreme values and makes the distributions more symmetric and easier to compare. 
# Log transformation can also help stabilize variability, making the data more suitable for statistical analysis.

