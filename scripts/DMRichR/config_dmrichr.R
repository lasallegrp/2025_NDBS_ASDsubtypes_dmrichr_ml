# ============================================================
# config_dmrichr.R
# Purpose: DMRichR-specific configuration
# Usage:   source(here::here("scripts", "DMRichR", "config_dmrichr.R"))
# ============================================================

library(here)

ROOT <- here::here()

# --- DMRichR result paths ---
DMRICHR_RESULTS <- list(
  base    = file.path(ROOT, "results", "DMRichR"),
  males   = file.path(ROOT, "results", "DMRichR", "males"),
  females = file.path(ROOT, "results", "DMRichR", "females")
)

# --- Bismark raw data location ---
BISMARK_RAW_DIR <- "/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_epigenerator_George/cytosine_reports_combined"
BISMARK_EXT     <- "_merged_name_sorted.deduplicated.bismark.cov.gz.CpG_report.txt.gz"

# --- Sample ID inputs from LPA pipeline ---
DMRICHR_INPUTS <- list(
  males   = file.path(ROOT, "results", "LPA", "10_DMRichR_inputs", "males"),
  females = file.path(ROOT, "results", "LPA", "10_DMRichR_inputs", "females")
)

# --- Create directories ---
invisible(lapply(DMRICHR_RESULTS, function(p) {
  dir.create(p, showWarnings = FALSE, recursive = TRUE)
}))

cat("DMRichR config loaded. Project root:", ROOT, "\n")