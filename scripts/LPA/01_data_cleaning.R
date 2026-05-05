# ============================================================
# 01_data_cleaning.R
# Purpose: Load raw data, remove failed samples, merge batches
# Output:  results/LPA/01_cleaning/ASDandTD_cleaned.rds
#          results/LPA/01_cleaning/
# ============================================================

library(tidyverse)
library(writexl)
source(here::here("scripts", "LPA", "config.R"))

# --- Load raw data ---
ASDandTDGrandChildren <- read_csv(
  file.path(PATHS$raw_data, RAW_FILENAME)
)
cat("Raw data dimensions:", dim(ASDandTDGrandChildren), "\n")

# --- Fix family ID ---
ASDandTDGrandChildren$subject_id <- gsub("PRI", "", ASDandTDGrandChildren$subject_id)
names(ASDandTDGrandChildren)[names(ASDandTDGrandChildren) == "subject_id"] <- "family_id"

# --- Remove failed library prep samples ---
rows_removed <- ASDandTDGrandChildren[
  ASDandTDGrandChildren$Freezerworks_ID %in% FAILED_LIBRARY_PREP_IDS,
  c("Freezerworks_ID", "family_id")
]
write_xlsx(rows_removed,
           file.path(RESULTS$cleaning, "removed_failed_library_prep.xlsx"))

ASDandTDGrandChildren_remove_8 <- ASDandTDGrandChildren[
  !ASDandTDGrandChildren$Freezerworks_ID %in% FAILED_LIBRARY_PREP_IDS, ]
cat("After removing failed samples:", nrow(ASDandTDGrandChildren_remove_8), "rows\n")

# --- Separate George and Jules batches ---
Jules_subset  <- ASDandTDGrandChildren_remove_8[
  ASDandTDGrandChildren_remove_8$Batch == "Jules", ]
head(Jules_subset)
dim(Jules_subset)

George_subset <- ASDandTDGrandChildren_remove_8[
  grepl("George", ASDandTDGrandChildren_remove_8$Batch), ]
head(George_subset)
dim(George_subset)

# --- Fix IBC vs Freezerworks ID ---
Jules_subset_1a <- Jules_subset[, !colnames(Jules_subset) %in% "Freezerworks_ID"]
colnames(Jules_subset_1a)[colnames(Jules_subset_1a) == "IBC"] <- "Freezerworks_ID"

George_subset$Freezerworks_ID[George_subset$IBC == 117981] <- 117981
George_subset_1a <- George_subset[, !colnames(George_subset) %in% "IBC"]

# --- Merge George and Jules ---
merged_data_raw <- merge(George_subset_1a, Jules_subset_1a,
                         by = "Freezerworks_ID", all = TRUE,
                         suffixes = c("_George", "_Jules"))

common_cols <- intersect(
  gsub("_George$", "", grep("_George$", names(merged_data_raw), value = TRUE)),
  gsub("_Jules$",  "", grep("_Jules$",  names(merged_data_raw), value = TRUE))
)
for (col in common_cols) {
  merged_data_raw[[col]] <- coalesce(
    merged_data_raw[[paste0(col, "_George")]],
    merged_data_raw[[paste0(col, "_Jules")]]
  )
}
ASDandTD_cleaned <- merged_data_raw %>%
  dplyr::select(-matches("_George$|_Jules$"))

# --- Remove all-NA columns and unnecessary columns ---
ASDandTD_cleaned <- ASDandTD_cleaned[,
  colSums(is.na(ASDandTD_cleaned)) != nrow(ASDandTD_cleaned)]
ASDandTD_cleaned <- ASDandTD_cleaned[,
  !colnames(ASDandTD_cleaned) %in% "C_NEXPDQ9spec"]

# --- Save cleaning log ---
cleaning_log <- data.frame(
  step = c("Raw data", "After removing failed samples", "After merging batches"),
  n_rows = c(nrow(ASDandTDGrandChildren),
             nrow(ASDandTDGrandChildren_remove_8),
             nrow(ASDandTD_cleaned))
)
write.csv(cleaning_log,
          file.path(RESULTS$cleaning, "cleaning_log.csv"),
          row.names = FALSE)

# --- Save cleaned data ---
saveRDS(ASDandTD_cleaned,
        file.path(RESULTS$cleaning, "ASDandTD_cleaned.rds"))
cat("Saved: ASDandTD_cleaned.rds —", nrow(ASDandTD_cleaned), "rows\n")