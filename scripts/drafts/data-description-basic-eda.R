library(tidyverse)



load('data/biomarker-clean.RData')

biomarker_clean %>%
  group_by(group) %>%
  summarize(n = n(),
            ados_n = sum(!is.na(ados)),
            ados_mean = mean(ados, na.rm = TRUE),
            ados_sd = sd(ados, na.rm = TRUE))

n_proteins <- ncol(biomarker_clean) - 2
n_proteins

ncol(biomarker_clean)