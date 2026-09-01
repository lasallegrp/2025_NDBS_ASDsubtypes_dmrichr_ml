# scripts/LPA/README.md

## Run order

Scripts must be run in order. Each script depends on outputs from the previous.

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| 01_data_cleaning.R | Load raw data, remove failed samples, merge batches | data/metadata/raw/ | data/processed/ASDandTD_cleaned.rds |
| 02_recode_variables.R | Rename and recode all variables | ASDandTD_cleaned.rds | ASDandTD_recoded.rds |
| 03_subset_neuromeasures.R | Split by sex/diagnosis, create LPA inputs | ASDandTD_recoded.rds | data/processed/males/ and females/ |
| 04_LPA_ASD_males.R | LPA model comparison + 3-class solution for ASD males | males/ | results/LPA/04_ASD_males/ |
| 05_LPA_TD_males.R | LPA for TD males, identify refined comparison group | males/ | results/LPA/05_TD_males/ |
| 06_LPA_combined_ASD_TD_males.R | Combine and compare ASD + TD male profiles | males/ | results/LPA/06_combined_males/ |
| 07_LPA_ASD_females.R | LPA for ASD females | females/ | results/LPA/07_ASD_females/ |
| 08_LPA_TD_females.R | LPA for TD females | females/ | results/LPA/08_TD_females/ |
| 09_LPA_combined_ASD_TD_females.R | Combine and compare ASD + TD female profiles | females/ | results/LPA/09_combined_females/ |

## Before running females

You must update `NA_SAMPLE_IDS$females` in `config.R` before running
scripts 07-09. See the error message in `03_subset_neuromeasures.R`
for instructions — it will print the specific IDs to add.

## Interactive testing

To test any script line by line in RStudio:
1. Open the project (.Rproj file) so `here::here()` resolves correctly
2. Run `source(here::here("scripts", "LPA", "config.R"))` once
3. All PATHS, RESULTS, and LPA_CONFIG objects are now in your environment
4. Run any line from any script interactively