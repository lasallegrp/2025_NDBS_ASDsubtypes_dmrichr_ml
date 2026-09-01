# ============================================================
# 03_subset_neuromeasures.R
# Purpose: Split by sex and diagnosis, create named LPA inputs
# Input:   results/LPA/02_recoding/ASDandTD_recoded.rds
# Output:  data/processed/males/ and data/processed/females/
#          results/LPA/03_subsets/
# ============================================================

library(tidyverse)
source(here::here("tests", "scripts", "LPA", "config.R"))

ASDandTD_recoded <- readRDS(file.path(RESULTS$recoding, "ASDandTD_recoded.rds"))

build_subsets <- function(data, sex_value, sex_label) {

  # --- Validate sex_label ---
  if (!sex_label %in% c("males", "females")) {
    stop(paste0(
      "\n================================================================\n",
      "ERROR: sex_label '", sex_label, "' is not recognized.\n",
      "Expected 'males' or 'females'.\n",
      "================================================================\n"
    ))
  }

  sex_data <- data[data$Child_Sex == sex_value, ]

  core_cols <- c(
    "Freezerworks_ID", "family_id", "Child_ASD_Status", "Child_TD_Status",
    "F2_SCQ_Total", "Autism_Severity_Score", "ADOS_Age_in_Months",
    "ADOS2_Social_Affect", "ADOS2_Restricted_Repetitive_Behaviors",
    "ADOS2_Total_Score", "ADI_R_Social_Interaction_Score",
    "ADI_R_Verbal_Communication_Score", "ADI_R_Nonverbal_Communication_Score",
    "ADI_R_Restricted_Behavior_Score", "ADI_R_Onset_Total_by_36m",
    "Mullen_Chronological_Age", "Mullen_Visual_Reception", "Mullen_Fine_Motor",
    "Mullen_Receptive_Language", "Mullen_Expressive_Language",
    "Mullen_Early_Learning_Composite", "Mullen_Verbal_Score",
    "Mullen_Nonverbal_Score", "Vineland_Chronological_Age",
    "Vineland_Communication", "Vineland_Daily_Living_Skills",
    "Vineland_Socialization", "Vineland_Motor_Skills", "Vineland_Composite_Score"
  )

  neuro_data <- sex_data[, core_cols]

  # --- ASD subset ---
  ASD_data <- neuro_data[neuro_data$Child_ASD_Status == 1, ]
  ASD_data <- ASD_data[, !names(ASD_data) %in% "Child_TD_Status"]
  rownames(ASD_data) <- ASD_data$Freezerworks_ID

  # --- NA removal with error checking ---
  # Note: ADI_R_Verbal and Nonverbal excluded from this check
  # because nonverbal kids will only have one or the other —
  # handled downstream by creating ADI_R_Communication_Score
  ados_adir_cols <- c(
    "Autism_Severity_Score", "ADOS_Age_in_Months",
    "ADOS2_Social_Affect", "ADOS2_Restricted_Repetitive_Behaviors",
    "ADOS2_Total_Score", "ADI_R_Social_Interaction_Score",
    "ADI_R_Restricted_Behavior_Score", "ADI_R_Onset_Total_by_36m"
  )

  na_check <- ASD_data %>%
    dplyr::select(Freezerworks_ID, all_of(ados_adir_cols)) %>%
    dplyr::filter(if_any(all_of(ados_adir_cols), is.na))

  ids_with_missing <- NA_SAMPLE_IDS[[sex_label]]

  # --- Female guard ---
  # Triggers only on literal NA_character_ placeholder
  # c() means "confirmed no NAs" and passes through correctly
  # This is necessary because all(is.na(c())) is vacuously TRUE in R
  if (sex_label == "females" &&
      (identical(ids_with_missing, NA_character_) ||
       (length(ids_with_missing) == 1 && is.na(ids_with_missing)))) {
    stop(paste0(
      "\n================================================================\n",
      "ERROR: Female missing sample IDs have not been specified.\n",
      "----------------------------------------------------------------\n",
      "Before running females:\n",
      "  1. Review NAs detected below\n",
      "  2. Open scripts/LPA/config.R\n",
      "  3. Update NA_SAMPLE_IDS$females with confirmed IDs\n",
      "     Or set to c() if none found\n",
      "----------------------------------------------------------------\n",
      if (nrow(na_check) > 0) {
        paste0(
          "DETECTED: These female samples have NAs in ADOS/ADI-R:\n",
          paste(capture.output(print(na_check)), collapse = "\n"), "\n"
        )
      } else {
        "GOOD NEWS: No NAs detected. Set NA_SAMPLE_IDS$females <- c()\n"
      },
      "================================================================\n"
    ))
  }

  # --- Error if NAs exist that are not listed ---
  if (nrow(na_check) > 0) {
    unhandled <- setdiff(na_check$Freezerworks_ID, ids_with_missing)
    if (length(unhandled) > 0) {
      stop(paste0(
        "\n================================================================\n",
        "ERROR: Unhandled NAs in ADOS/ADI-R for ", sex_label, ":\n",
        paste(unhandled, collapse = ", "), "\n",
        "Add these to NA_SAMPLE_IDS$", sex_label, " in config.R\n",
        "================================================================\n"
      ))
    }
  }

  # --- Warn if listed IDs have no NAs (stale entries) ---
  if (!identical(ids_with_missing, NA_character_) &&
      length(ids_with_missing) > 0) {
    overclaimed <- setdiff(ids_with_missing, na_check$Freezerworks_ID)
    if (length(overclaimed) > 0) {
      warning(paste0(
        "\nWARNING: IDs in NA_SAMPLE_IDS$", sex_label, " have no NAs ",
        "in current data: ", paste(overclaimed, collapse = ", "), "\n",
        "Review NA_SAMPLE_IDS in config.R\n"
      ))
    }
    ASD_data <- ASD_data[!ASD_data$Freezerworks_ID %in% ids_with_missing, ]
  }

  cat("ASD", sex_label, "after NA removal: n =", nrow(ASD_data), "\n")

  # --- ADI-R communication score ---
  # Uses verbal score where available, nonverbal otherwise
  ASD_data$ADI_R_Communication_Score <- ifelse(
    is.na(ASD_data$ADI_R_Verbal_Communication_Score),
    ASD_data$ADI_R_Nonverbal_Communication_Score,
    ASD_data$ADI_R_Verbal_Communication_Score
  )

  # --- LPA input: ADOS + ADI-R only ---
  ASD_lpa_input_ADOS_ADIR <- ASD_data[, c(
    "Autism_Severity_Score", "ADOS_Age_in_Months",
    "ADOS2_Social_Affect", "ADOS2_Restricted_Repetitive_Behaviors",
    "ADOS2_Total_Score", "ADI_R_Social_Interaction_Score",
    "ADI_R_Restricted_Behavior_Score", "ADI_R_Onset_Total_by_36m",
    "ADI_R_Communication_Score"
  )]

  # --- LPA input: ADOS + ADI-R + Mullen composite + Vineland composite ---
  ASD_lpa_input_ADOS_ADIR_Mullen_Vine <- ASD_data[, c(
    "Autism_Severity_Score", "ADOS_Age_in_Months",
    "ADOS2_Social_Affect", "ADOS2_Restricted_Repetitive_Behaviors",
    "ADOS2_Total_Score", "ADI_R_Social_Interaction_Score",
    "ADI_R_Restricted_Behavior_Score", "ADI_R_Onset_Total_by_36m",
    "ADI_R_Communication_Score", "Mullen_Early_Learning_Composite",
    "Vineland_Composite_Score"
  )]

  # --- Full Mullen/Vineland for post-hoc plotting (not LPA input) ---
  ASD_mullen_vineland_full <- ASD_data[, c(
    "Freezerworks_ID", "Mullen_Chronological_Age", "Mullen_Visual_Reception",
    "Mullen_Fine_Motor", "Mullen_Receptive_Language", "Mullen_Expressive_Language",
    "Mullen_Early_Learning_Composite", "Mullen_Verbal_Score", "Mullen_Nonverbal_Score",
    "Vineland_Chronological_Age", "Vineland_Communication", "Vineland_Daily_Living_Skills",
    "Vineland_Socialization", "Vineland_Motor_Skills", "Vineland_Composite_Score"
  )]

  # --- TD subset ---
  TD_data <- neuro_data[neuro_data$Child_TD_Status == 1, ]
  TD_data <- TD_data[, !names(TD_data) %in% "Child_ASD_Status"]
  rownames(TD_data) <- TD_data$Freezerworks_ID

  # --- Identify which TD samples are removed by na.omit ---
  TD_scq_check <- TD_data[, c(
    "Freezerworks_ID", "F2_SCQ_Total",
    "Mullen_Early_Learning_Composite",
    "Vineland_Composite_Score"
  )]

  missing_td <- TD_scq_check[!complete.cases(TD_scq_check), ]

  if (nrow(missing_td) > 0) {
    cat("\nTD sample(s) removed due to missing data for", sex_label, ":\n")
    print(missing_td)
    write.csv(missing_td,
              file.path(RESULTS$subsets,
                        paste0(sex_label, "_TD_removed_missing_data.csv")),
              row.names = FALSE)
  } else {
    cat("\nNo TD samples removed for", sex_label, "\n")
  }

  # --- LPA input: Mullen + Vineland + SCQ ---
  TD_lpa_input_Mullen_Vine_SCQ <- TD_data[, c(
    "Mullen_Early_Learning_Composite", "Vineland_Composite_Score", "F2_SCQ_Total"
  )] %>% na.omit()

  # --- Save subsets to data/processed ---
  out <- PATHS[[sex_label]]

  saveRDS(ASD_lpa_input_ADOS_ADIR,
          file.path(out, "ASD_lpa_input_ADOS_ADIR.rds"))
  saveRDS(ASD_lpa_input_ADOS_ADIR_Mullen_Vine,
          file.path(out, "ASD_lpa_input_ADOS_ADIR_Mullen_Vine.rds"))
  saveRDS(ASD_mullen_vineland_full,
          file.path(out, "ASD_mullen_vineland_full.rds"))
  saveRDS(TD_lpa_input_Mullen_Vine_SCQ,
          file.path(out, "TD_lpa_input_Mullen_Vine_SCQ.rds"))

  # --- Save subset summary to results ---
  subset_log <- data.frame(
    group           = c("ASD", "TD"),
    n_before_filter = c(sum(neuro_data$Child_ASD_Status == 1, na.rm = TRUE),
                        sum(neuro_data$Child_TD_Status  == 1, na.rm = TRUE)),
    n_final         = c(nrow(ASD_lpa_input_ADOS_ADIR),
                        nrow(TD_lpa_input_Mullen_Vine_SCQ)),
    n_removed       = c(sum(neuro_data$Child_ASD_Status == 1, na.rm = TRUE) -
                          nrow(ASD_lpa_input_ADOS_ADIR),
                        sum(neuro_data$Child_TD_Status  == 1, na.rm = TRUE) -
                          nrow(TD_lpa_input_Mullen_Vine_SCQ))
  )

  write.csv(subset_log,
            file.path(RESULTS$subsets,
                      paste0(sex_label, "_subset_log.csv")),
            row.names = FALSE)

  cat("Saved subsets for", sex_label,
      "| ASD n =", nrow(ASD_lpa_input_ADOS_ADIR),
      "| TD n =",  nrow(TD_lpa_input_Mullen_Vine_SCQ), "\n")
}

# --- Run for males and females ---
build_subsets(ASDandTD_recoded, sex_value = 1, sex_label = "males")
build_subsets(ASDandTD_recoded, sex_value = 0, sex_label = "females")