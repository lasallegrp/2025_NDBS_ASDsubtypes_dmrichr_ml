# ============================================================
# config.R
# Purpose: Central configuration — all paths defined here
# Usage:   source(here::here("scripts", "LPA", "config.R"))
#          at the top of every LPA script
# ============================================================
# ============================================================
# Package availability check
# Runs automatically when config.R is sourced
# ============================================================

required_packages <- c(
  "tidyLPA", "here", "tidyverse", "writexl", "lubridate",
  "ggpubr", "rstatix", "FSA", "patchwork", "car", "rcompanion",
  "ggplot2", "readxl", "openxlsx", "magrittr"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(paste0(
    "\n================================================================\n",
    "ERROR: The following required packages are not available:\n",
    paste(" -", missing_packages, collapse = "\n"), "\n",
    "----------------------------------------------------------------\n",
    "If using pixi, run:\n",
    "  pixi install\n",
    "Then verify with:\n",
    "  pixi run check-lpa-pkgs\n",
    "If a package is still missing, add it to pixi.toml under\n",
    "the LPA analysis packages section and re-run pixi install.\n",
    "================================================================\n"
  ))
} else {
  cat("All required packages available.\n")
}

library(here)

# --- Root (resolves to 2025_NDBS_ASDsubtypes_dmrichr_ml_George/tests) ---
# ROOT <- here::here() for main root but I wanted to resolve to the tests directory for testing purposes
ROOT <- file.path(here::here(), "tests")
# --- Data paths ---
PATHS <- list(
  raw_data   = file.path(ROOT, "data", "metadata", "raw"),
  processed  = file.path(ROOT, "data", "processed"),
  males      = file.path(ROOT, "data", "processed", "males"),
  females    = file.path(ROOT, "data", "processed", "females")
)

# --- Results paths (mirrors script numbering) ---
RESULTS <- list(
  # --- LPA results ---
  cleaning         = file.path(ROOT, "results", "LPA", "01_cleaning"),
  recoding         = file.path(ROOT, "results", "LPA", "02_recoding"),
  subsets          = file.path(ROOT, "results", "LPA", "03_subsets"),
  asd_males        = file.path(ROOT, "results", "LPA", "04_ASD_males"),
  td_males         = file.path(ROOT, "results", "LPA", "05_TD_males"),
  combined_males   = file.path(ROOT, "results", "LPA", "06_combined_males"),
  asd_females      = file.path(ROOT, "results", "LPA", "07_ASD_females"),
  td_females       = file.path(ROOT, "results", "LPA", "08_TD_females"),
  combined_females = file.path(ROOT, "results", "LPA", "09_combined_females"),
  dmrichr_inputs   = file.path(ROOT, "results", "LPA", "10_DMRichR_inputs")
)

# --- Raw data filename ---
RAW_FILENAME <- "ASDandTDChildrenwithDNAfromDBS_M5514_13NOV23.csv"

# --- Sample IDs to remove (failed library prep) ---
FAILED_LIBRARY_PREP_IDS <- c(128203, 138483, 126305, 128165,
                              138200, 138326, 146703, 146704)

# --- Sex-specific NA sample IDs in ADOS/ADI-R ---
# Males:   confirmed from data exploration
# Females: update after inspecting female ASD data
NA_SAMPLE_IDS <- list(
  males   = c(125864, 138315),
  females = c()               # confirmed no NAs in female ADOS/ADI-R
)

# --- LPA settings ---
# UPDATE — run models separately with different K limits
# --- LPA settings ---
LPA_CONFIG <- list(
  models = c(1, 2),

  # --- K limits by sex ---
  # Males (n=132): can support more classes
  # Females (n=27): capped lower to avoid convergence failures
  max_K_ASD_m1 = list(males = 6, females = 3),
  max_K_ASD_m2 = list(males = 4, females = 2),
  max_K_TD_m1 = list(males = 4, females = 3),
  max_K_TD_m2 = list(males = 4, females = 2),  # females cap at 2 for Model 2
  
  # --- Chosen solutions by sex ---
  # Update females after reviewing fit table from script 07
  best_K_ASD   = list(males = 3, females = 2),
  best_K_TD    = list(males = 2, females = 1),

  # --- Model number chosen as best ---
  best_model_ASD = list(males = 1, females = 1),
  best_model_TD  = list(males = 1, females = 1),

  # --- File names ---
  asd_input_file   = "ASD_lpa_input_ADOS_ADIR_Mullen_Vine.rds",
  asd_classed_file = list(
    males   = "ASD_lpa_classed_Mullen_Vine_k3.rds",
    females = "ASD_lpa_classed_Mullen_Vine_k2.rds"
  ),
  td_classed_file = list(
  males   = "TD_lpa_classed_k2.rds",
  females = "TD_lpa_classed_k1.rds"   # all females, no subgrouping
  )
)

#skip LPA
# Add a flag to skip TD LPA for females
skip_TD_LPA = list(males = FALSE, females = TRUE)

# --- Create all directories ---
invisible(lapply(c(PATHS, RESULTS), function(p) {
  dir.create(p, showWarnings = FALSE, recursive = TRUE)
}))

cat("Config loaded. Project root:", ROOT, "\n")