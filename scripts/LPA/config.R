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

# --- Root (resolves to 2025_NDBS_ASDsubtypes_dmrichr_ml_George/) ---
ROOT <- here::here()

# --- Data paths ---
PATHS <- list(
  raw_data   = file.path(ROOT, "data", "metadata", "raw"),
  processed  = file.path(ROOT, "data", "processed"),
  males      = file.path(ROOT, "data", "processed", "males"),
  females    = file.path(ROOT, "data", "processed", "females")
)

# --- Results paths (mirrors script numbering) ---
RESULTS <- list(
  cleaning          = file.path(ROOT, "results", "LPA", "01_cleaning"),
  recoding          = file.path(ROOT, "results", "LPA", "02_recoding"),
  subsets           = file.path(ROOT, "results", "LPA", "03_subsets"),
  asd_males         = file.path(ROOT, "results", "LPA", "04_ASD_males"),
  td_males          = file.path(ROOT, "results", "LPA", "05_TD_males"),
  combined_males    = file.path(ROOT, "results", "LPA", "06_combined_males"),
  asd_females       = file.path(ROOT, "results", "LPA", "07_ASD_females"),
  td_females        = file.path(ROOT, "results", "LPA", "08_TD_females"),
  combined_females  = file.path(ROOT, "results", "LPA", "09_combined_females")
)

# --- Raw data filename ---
RAW_FILENAME <- "ASDandTDChildrenwithDNAfromDBS_M5514_13NOV23 (3) (1)(ASDandTDChildrenwithDNAfromDBS_) (1).csv"

# --- Sample IDs to remove (failed library prep) ---
FAILED_LIBRARY_PREP_IDS <- c(128203, 138483, 126305, 128165,
                              138200, 138326, 146703, 146704)

# --- Sex-specific NA sample IDs in ADOS/ADI-R ---
# Males:   confirmed from data exploration
# Females: update after inspecting female ASD data
NA_SAMPLE_IDS <- list(
  males   = c(125864, 138315),
  females = NA_character_     # <-- UPDATE BEFORE RUNNING FEMALES
)

# --- LPA settings ---
LPA_CONFIG <- list(
  models     = c(1, 2),
  max_K_ASD  = 6,
  max_K_TD   = 4,
  best_K_ASD = 3,
  best_K_TD  = 2
)

# --- Create all directories ---
invisible(lapply(c(PATHS, RESULTS), function(p) {
  dir.create(p, showWarnings = FALSE, recursive = TRUE)
}))

cat("Config loaded. Project root:", ROOT, "\n")