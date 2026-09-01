# ============================================================
# 10_prepare_DMRichR_inputs.R
# Purpose: Generate sample_info.xlsx and sample ID text files
#          for all pairwise DMRichR comparisons
#          ASD subtypes vs each other + ASD subtypes vs TD
# Input:   data/processed/males/ and data/processed/females/
# Output:  results/LPA/10_DMRichR_inputs/males/
#          results/LPA/10_DMRichR_inputs/females/
# ============================================================

library(tidyverse)
library(writexl)
source(here::here("tests", "scripts", "LPA", "config.R"))

# --- Add output paths ---
dmrichr_out <- list(
  males   = file.path(ROOT, "results", "LPA", "10_DMRichR_inputs", "males"),
  females = file.path(ROOT, "results", "LPA", "10_DMRichR_inputs", "females")
)
invisible(lapply(dmrichr_out, function(p) {
  dir.create(file.path(p, "sample_info"),    showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(p, "sample_ids"),     showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(p, "comparison_log"), showWarnings = FALSE, recursive = TRUE)
}))

# ============================================================
# Helper functions
# ============================================================

# --- Validate sample names: no underscores allowed ---
validate_names <- function(names, context) {
  bad <- names[grepl("_", names)]
  if (length(bad) > 0) {
    stop(paste0(
      "\n================================================================\n",
      "ERROR: Sample names contain underscores in ", context, ".\n",
      "DMRichR matches on text before the first underscore.\n",
      "Affected names:\n",
      paste(" -", bad, collapse = "\n"), "\n",
      "----------------------------------------------------------------\n",
      "Remove underscores from sample names before proceeding.\n",
      "================================================================\n"
    ))
  }
}

# --- Validate alphabetical ordering ---
# DMRichR requires control label to come BEFORE experimental alphabetically
# to get results in the correct direction (experimental vs control)
validate_alpha_order <- function(control_label, experimental_label, context) {
  if (control_label >= experimental_label) {
    stop(paste0(
      "\n================================================================\n",
      "ERROR: Label ordering violation in ", context, ".\n",
      "----------------------------------------------------------------\n",
      "DMRichR requires the control label to come BEFORE the\n",
      "experimental label alphabetically.\n",
      "Control label:      ", control_label, "\n",
      "Experimental label: ", experimental_label, "\n",
      "----------------------------------------------------------------\n",
      "Update TD_LABEL or ASD_SUBTYPE_LABELS in script 10 so that\n",
      "the TD label comes before the ASD label alphabetically.\n",
      "Current setup uses Neurotypical (N) before Subtype (S).\n",
      "================================================================\n"
    ))
  }
  invisible(TRUE)
}

# --- Check minimum sample size ---
validate_min_samples <- function(df, comparison_name, min_per_group = 2) {
  counts <- table(df$Diagnosis)
  below_min <- counts[counts < min_per_group]
  if (length(below_min) > 0) {
    warning(paste0(
      "\nWARNING: Some groups have fewer than ", min_per_group,
      " samples in comparison: ", comparison_name, "\n",
      "DMRichR requires minimum 4 samples total (2 per group).\n",
      "Groups below minimum:\n",
      paste(names(below_min), "n =", below_min, collapse = "\n"), "\n",
      "This comparison will be skipped.\n"
    ))
    return(FALSE)
  }
  return(TRUE)
}

# --- Write one comparison ---
write_comparison <- function(df, comparison_name, out_dir) {

  if (!validate_min_samples(df, comparison_name)) return(invisible(NULL))

  sample_info <- df %>%
    dplyr::select(Name, Diagnosis) %>%
    arrange(Diagnosis, Name)

  write_xlsx(sample_info,
             file.path(out_dir, "sample_info",
                       paste0(comparison_name, "_sample_info.xlsx")))

  write_lines(sample_info$Name,
              file.path(out_dir, "sample_ids",
                        paste0(comparison_name, "_sample_ids.txt")))

  cat("  Written:", comparison_name,
      "| n =", nrow(sample_info),
      "| groups:", paste(names(table(sample_info$Diagnosis)),
                         table(sample_info$Diagnosis),
                         sep = "=", collapse = ", "), "\n")

  invisible(sample_info)
}

# ============================================================
# Main function — runs for one sex
# ============================================================

prepare_dmrichr_inputs <- function(sex_label) {

  cat("\n================================================================\n")
  cat("Preparing DMRichR inputs for:", sex_label, "\n")
  cat("================================================================\n")

  out_dir <- dmrichr_out[[sex_label]]

  # --- Label configuration ---
  # Neurotypical (N) comes before Subtype (S) alphabetically
  # This ensures DMRichR results are in the correct direction:
  # Subtype_ASD vs Neurotypical (experimental vs control)

  TD_LABEL <- "Neurotypical"

  # ASD subtype labels — automatically adjusts for females (2 classes)
  ASD_SUBTYPE_LABELS <- if (sex_label == "females") {
    list(
      "1" = "Subtype1_ASD",
      "2" = "Subtype2_ASD"
    )
  } else {
    list(
      "1" = "Subtype1_ASD",
      "2" = "Subtype2_ASD",
      "3" = "Subtype3_ASD"
    )
  }

  cat("TD label:  ", TD_LABEL, "\n")
  cat("ASD labels:", paste(unlist(ASD_SUBTYPE_LABELS), collapse = ", "), "\n")

  # --- Load ASD classed data ---
  asd_classed_path <- file.path(PATHS[[sex_label]],
                                LPA_CONFIG$asd_classed_file[[sex_label]])

  if (!file.exists(asd_classed_path)) {
    stop(paste0(
      "\n================================================================\n",
      "ERROR: ASD classed file not found for ", sex_label, ":\n",
      asd_classed_path, "\n",
      "Have you run script 04/07 for ", sex_label, "?\n",
      "================================================================\n"
    ))
  }

  ASD_classed <- readRDS(asd_classed_path)
  cat("ASD", sex_label, "loaded | n =", nrow(ASD_classed), "\n")
  cat("ASD class distribution:\n")
  print(table(ASD_classed$Class))

  # --- Load TD refined comparison group ---
  td_refined_path <- file.path(PATHS[[sex_label]],
                               "TD_refined_comparison_group.rds")

  if (!file.exists(td_refined_path)) {
    stop(paste0(
      "\n================================================================\n",
      "ERROR: TD refined comparison group not found for ", sex_label, ":\n",
      td_refined_path, "\n",
      "Have you run script 05/08 for ", sex_label, "?\n",
      "================================================================\n"
    ))
  }

  TD_refined <- readRDS(td_refined_path)
  cat("TD", sex_label, "refined group loaded | n =", nrow(TD_refined), "\n")

  # --- Build ASD sample dataframe ---
  ASD_classed$Name <- as.character(rownames(ASD_classed))
  validate_names(ASD_classed$Name, paste("ASD", sex_label))

  ASD_classed$Diagnosis <- unlist(
    ASD_SUBTYPE_LABELS[as.character(ASD_classed$Class)]
  )

  cat("ASD diagnosis label mapping:\n")
  print(table(ASD_classed$Diagnosis))

  # --- Build TD sample dataframe ---
  TD_refined$Name <- as.character(rownames(TD_refined))
  validate_names(TD_refined$Name, paste("TD", sex_label))
  TD_refined$Diagnosis <- TD_LABEL

  cat("TD label:", TD_LABEL, "| n =", nrow(TD_refined), "\n")

  # --- Validate alphabetical ordering for all comparisons ---
  for (asd_label in ASD_SUBTYPE_LABELS) {
    validate_alpha_order(TD_LABEL, asd_label,
                         paste(TD_LABEL, "vs", asd_label))
  }
  asd_labels <- unlist(ASD_SUBTYPE_LABELS)
  if (length(asd_labels) > 1) {
    asd_pairs <- combn(asd_labels, 2, simplify = FALSE)
    for (pair in asd_pairs) {
      validate_alpha_order(pair[[1]], pair[[2]],
                           paste(pair[[1]], "vs", pair[[2]]))
    }
  }

  # ============================================================
  # Generate comparisons
  # ============================================================

  comparison_log <- list()
  # --- 0. All ASD participants vs Neurotypical ---
  cat("\n--- Total ASD vs Neurotypical ---\n")

  TOTAL_ASD_LABEL <- "Total_ASD"

  # Confirm the labels give the intended DMRichR direction:
  # Neurotypical is control and Total_ASD is experimental.
  validate_alpha_order(
    control_label      = TD_LABEL,
    experimental_label = TOTAL_ASD_LABEL,
    context            = paste(TD_LABEL, "vs", TOTAL_ASD_LABEL)
  )

  # Combine every ASD subtype into one ASD group
  asd_total_subset <- ASD_classed %>%
    dplyr::transmute(
      Name      = Name,
      Diagnosis = TOTAL_ASD_LABEL
    )

  # Use the same sex-specific refined TD comparison group
  td_total_subset <- TD_refined %>%
    dplyr::transmute(
      Name      = Name,
      Diagnosis = TD_LABEL
    )

  total_comparison_df <- bind_rows(
    asd_total_subset,
    td_total_subset
  )

  total_comparison_name <- paste0(
    TOTAL_ASD_LABEL,
    "_vs_",
    TD_LABEL
  )

  total_result <- write_comparison(
    df              = total_comparison_df,
    comparison_name = total_comparison_name,
    out_dir         = out_dir
  )

  if (!is.null(total_result)) {
    comparison_log[[total_comparison_name]] <- data.frame(
      comparison = total_comparison_name,
      group1     = TOTAL_ASD_LABEL,
      group2     = TD_LABEL,
      n_group1   = sum(
        total_comparison_df$Diagnosis == TOTAL_ASD_LABEL
      ),
      n_group2   = sum(
        total_comparison_df$Diagnosis == TD_LABEL
      ),
      n_total    = nrow(total_comparison_df),
      sex        = sex_label,
      type       = "Total_ASD_vs_TD"
    )
  }
  # --- 1. Each ASD subtype vs TD ---
  cat("\n--- ASD subtypes vs Neurotypical ---\n")
  for (class_num in names(ASD_SUBTYPE_LABELS)) {
    asd_label <- ASD_SUBTYPE_LABELS[[class_num]]

    asd_subset <- ASD_classed %>%
      dplyr::filter(Diagnosis == asd_label) %>%
      dplyr::select(Name, Diagnosis)

    td_subset <- TD_refined %>%
      dplyr::select(Name, Diagnosis)

    comparison_df   <- bind_rows(asd_subset, td_subset)
    comparison_name <- paste0(asd_label, "_vs_", TD_LABEL)

    result <- write_comparison(comparison_df, comparison_name, out_dir)

    if (!is.null(result)) {
      comparison_log[[comparison_name]] <- data.frame(
        comparison = comparison_name,
        group1     = asd_label,
        group2     = TD_LABEL,
        n_group1   = sum(comparison_df$Diagnosis == asd_label),
        n_group2   = sum(comparison_df$Diagnosis == TD_LABEL),
        n_total    = nrow(comparison_df),
        sex        = sex_label,
        type       = "ASD_vs_TD"
      )
    }
  }

  # --- 2. ASD subtypes vs each other ---
  cat("\n--- ASD subtypes pairwise ---\n")
  asd_class_nums <- names(ASD_SUBTYPE_LABELS)

  if (length(asd_class_nums) > 1) {
    asd_pairs <- combn(asd_class_nums, 2, simplify = FALSE)

    for (pair in asd_pairs) {
      label1 <- ASD_SUBTYPE_LABELS[[pair[[1]]]]
      label2 <- ASD_SUBTYPE_LABELS[[pair[[2]]]]

      # Ensure alphabetical order — swap if needed
      if (label1 > label2) {
        tmp    <- label1
        label1 <- label2
        label2 <- tmp
      }

      subset1 <- ASD_classed %>%
        dplyr::filter(Diagnosis == label1) %>%
        dplyr::select(Name, Diagnosis)

      subset2 <- ASD_classed %>%
        dplyr::filter(Diagnosis == label2) %>%
        dplyr::select(Name, Diagnosis)

      comparison_df   <- bind_rows(subset1, subset2)
      comparison_name <- paste0(label1, "_vs_", label2)

      result <- write_comparison(comparison_df, comparison_name, out_dir)

      if (!is.null(result)) {
        comparison_log[[comparison_name]] <- data.frame(
          comparison = comparison_name,
          group1     = label1,
          group2     = label2,
          n_group1   = sum(comparison_df$Diagnosis == label1),
          n_group2   = sum(comparison_df$Diagnosis == label2),
          n_total    = nrow(comparison_df),
          sex        = sex_label,
          type       = "ASD_vs_ASD"
        )
      }
    }
  } else {
    cat("Only 1 ASD class for", sex_label,
        "— skipping ASD vs ASD comparisons\n")
  }

  # --- Save comparison log ---
  comparison_log_df <- bind_rows(comparison_log)
  write_xlsx(comparison_log_df,
             file.path(out_dir, "comparison_log",
                       paste0(sex_label, "_comparison_log.xlsx")))
  write.csv(comparison_log_df,
            file.path(out_dir, "comparison_log",
                      paste0(sex_label, "_comparison_log.csv")),
            row.names = FALSE)

  cat("\n--- Summary for", sex_label, "---\n")
  cat("Total comparisons generated:", nrow(comparison_log_df), "\n")
  cat("ASD vs TD:  ",
      sum(comparison_log_df$type == "ASD_vs_TD"), "\n")
  cat("ASD vs ASD: ",
      sum(comparison_log_df$type == "ASD_vs_ASD"), "\n")
  print(comparison_log_df[, c("comparison", "n_group1",
                               "n_group2", "n_total")])

  invisible(comparison_log_df)
}

# ============================================================
# Run for males and females
# ============================================================

males_log   <- prepare_dmrichr_inputs("males")
females_log <- prepare_dmrichr_inputs("females")

# --- Combined log across both sexes ---
combined_log <- bind_rows(males_log, females_log)
write_xlsx(combined_log,
           file.path(ROOT, "results", "LPA", "10_DMRichR_inputs",
                     "all_comparisons_log.xlsx"))

cat("\n================================================================\n")
cat("All DMRichR inputs prepared.\n")
cat("Total comparisons:", nrow(combined_log), "\n")
cat("Males:  ", sum(combined_log$sex == "males"), "\n")
cat("Females:", sum(combined_log$sex == "females"), "\n")
cat("================================================================\n")