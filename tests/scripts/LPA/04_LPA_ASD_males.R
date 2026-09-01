# ============================================================
# 04_LPA_ASD_males.R
# Purpose: LPA model comparison + 3-class solution for ASD males
# ============================================================

# 04_LPA_ASD_males.R
sex_label  <- "males"
result_dir <- here::here("tests", "results", "LPA", "04_ASD_males")
source(here::here("tests", "scripts", "LPA", "config.R"))
source(here::here("tests", "scripts", "LPA", "LPA_ASD_template.R"))