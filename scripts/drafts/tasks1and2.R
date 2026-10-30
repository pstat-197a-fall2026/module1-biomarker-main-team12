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

# ---------------------------------
# TASK 2: EFFECT OF OUTLIER TRIMMING
# Temporarily remove the outlier trimming from preprocessing and do some exploratory analysis of the outlying values. Are there specific subjects (not values) that seem to be outliers? If so, are outliers more frequent in one group or the other? (Hint: consider tabulating the number of outlying values per subject.)
# ---------------------------------

# Log-transform and standardize protein levels W/O trimming
# Add a subject ID to track each individual
biomarker_untrimmed <- biomarker_raw %>%
  mutate(subject_id = row_number()) %>%
  mutate(across(all_of(proteins),
                ~ as.numeric(scale(log10(.x)))))

# Identify outlying values (3+ stdev from mean) & count them for each subject
subject_outliers <- biomarker_untrimmed %>%
  pivot_longer(cols = all_of(proteins),
               names_to = "protein",
               values_to = "z_score") %>%
  group_by(subject_id, group) %>%
  summarise(
    n_outliers = sum(abs(z_score) > 3, na.rm = TRUE),
    .groups = "drop")

# View subjects w/ the most outlying protein measurements
top_subjects <- subject_outliers %>%
  arrange(desc(n_outliers)) %>%
  slice_head(n = 10)

print(top_subjects)

# Identify subjects in the highest 5% of outlier counts
cutoff <- quantile(subject_outliers$n_outliers, 0.95)

high_outlier_subjects <- subject_outliers %>%
  filter(n_outliers >= cutoff) %>%
  arrange(desc(n_outliers))

print(high_outlier_subjects)

# Compare outlier counts between ASD & TD groups
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

# Plot # of outlying values per subject by group
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

# Plot the subjects w/ unusually high outlier counts
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

# Outlying measurements appear to be concentrated among certain subjects rather than evenly distributed across the dataset.
# Although TD subjects had more outliers on average, both TD and ASD had similar medians, suggesting that the typical subject in each group had similar numbers of outlying measurements. 
# These findings help explain why trimming may be useful for reducing the influence of extreme measurements on subsequent analyses.