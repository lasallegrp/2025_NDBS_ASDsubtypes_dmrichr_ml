# ============================================================
# 02_recode_variables.R
# Purpose: Recode and rename all variables
# Input:   results/LPA/01_cleaning/ASDandTD_cleaned.rds
# Output:  results/LPA/02_recoding/ASDandTD_recoded.rds
# ============================================================

library(tidyverse)
source(here::here("tests", "scripts", "LPA", "config.R"))

ASDandTD_cleaned <- readRDS(file.path(RESULTS$cleaning, "ASDandTD_cleaned.rds"))

# --- Recode diagnosis ---
table(ASDandTD_cleaned$dx2)
ASDandTD_cleaned$ASD <- ifelse(ASDandTD_cleaned$dx2 == 2, 1, 0)
ASDandTD_cleaned$TD  <- ifelse(ASDandTD_cleaned$dx2 == 5, 1, 0)

# --- Recode sex (1 = Male, 0 = Female) ---
table(ASDandTD_cleaned$correctsex)
ASDandTD_cleaned$correctsex <- ifelse(ASDandTD_cleaned$correctsex == "M", 1, 0)

# --- Recode race ---
ASDandTD_cleaned <- ASDandTD_cleaned %>%
  mutate(
    F2_Race_White        = case_when(ChildRace_5cat == 1 ~ 1L, !is.na(ChildRace_5cat) ~ 0L, TRUE ~ NA_integer_),
    F2_Race_Black        = case_when(ChildRace_5cat == 2 ~ 1L, !is.na(ChildRace_5cat) ~ 0L, TRUE ~ NA_integer_),
    F2_Race_Asian        = case_when(ChildRace_5cat == 4 ~ 1L, !is.na(ChildRace_5cat) ~ 0L, TRUE ~ NA_integer_),
    F2_Race_HispanicAny  = case_when(ChildRace_5cat == 6 ~ 1L, !is.na(ChildRace_5cat) ~ 0L, TRUE ~ NA_integer_),
    F2_Race_Multi_racial = case_when(ChildRace_5cat == 9 ~ 1L, !is.na(ChildRace_5cat) ~ 0L, TRUE ~ NA_integer_)
  )

# --- Rename variables ---
ASDandTD_recoded <- ASDandTD_cleaned %>%
  dplyr::rename(
    Child_ASD_Status                      = ASD,
    Child_TD_Status                       = TD,
    Child_Sex                             = correctsex,
    F2_BirthMonth                         = BirthMonth,
    F2_BirthYear                          = YOB,
    F2_BirthSeason                        = BirthSeason,
    F2_AgeMonths                          = ChildAgeMos,
    F2_ChildRace_5cat                     = ChildRace_5cat,
    Locale                                = locale,
    F2_PaymentDelivery                    = PaymentDelivery,
    F2_PregnancyWeightGain                = PregWtGainlb,
    F2_C_Section                          = Cesarean,
    F2_BirthLength                        = C_birthlength,
    F2_BirthWeight                        = C_BirthWeight,
    F2_APGAR1                             = C_APGAR1,
    F2_APGAR5                             = C_APGAR5,
    F2_APGAR10                            = C_APGAR10,
    F2_BirthDefect                        = BirthDefect,
    F2_NewbornExamDefect                  = C_NEXPDQ9defect,
    F2_DisorderCode                       = DisorderCode,
    F2_NICU                               = C_NICU,
    Autism_Severity_Score                 = AutismSevScore,
    F2_ADOS_Module                        = ADOSmodule,
    ADOS_Age_in_Months                    = ADOSAGEmos,
    ADOS2_Social_Affect                   = ADOS2_SA,
    ADOS2_Restricted_Repetitive_Behaviors = ADOS2_RRB,
    ADOS2_Total_Score                     = ADOS2_TOT,
    ADI_R_Social_Interaction_Score        = ADRtotalAs,
    ADI_R_Verbal_Communication_Score      = ADRtotalVBs,
    ADI_R_Nonverbal_Communication_Score   = ADRtotalNVBs,
    ADI_R_Restricted_Behavior_Score       = ADRtotalCs,
    ADI_R_Onset_Total_by_36m              = ADRtotalDs,
    Mullen_Chronological_Age              = MSLchronAge,
    Mullen_Visual_Reception               = MSLvrDQ,
    Mullen_Fine_Motor                     = MSLfmDQ,
    Mullen_Receptive_Language             = MSLrlDQ,
    Mullen_Expressive_Language            = MSLelDQ,
    Mullen_Early_Learning_Composite       = MSLelcDQ,
    Mullen_Verbal_Score                   = MSLverbDQ,
    Mullen_Nonverbal_Score                = MSLnoverbDQ,
    Vineland_Chronological_Age            = VINchronAge,
    Vineland_Communication                = VINcommDQ,
    Vineland_Daily_Living_Skills          = VINdlivDQ,
    Vineland_Socialization                = VINsociDQ,
    Vineland_Motor_Skills                 = VINmotoDQ,
    Vineland_Composite_Score              = VINcompDQ,
    F2_SCQ_Total                          = SCQtotal
  ) %>%
  dplyr::select(-dx2)

# --- Save recoding summary ---
recode_summary <- data.frame(
  variable = c("Child_ASD_Status", "Child_TD_Status", "Child_Sex"),
  n_1      = c(sum(ASDandTD_recoded$Child_ASD_Status == 1, na.rm = TRUE),
               sum(ASDandTD_recoded$Child_TD_Status  == 1, na.rm = TRUE),
               sum(ASDandTD_recoded$Child_Sex        == 1, na.rm = TRUE)),
  n_0      = c(sum(ASDandTD_recoded$Child_ASD_Status == 0, na.rm = TRUE),
               sum(ASDandTD_recoded$Child_TD_Status  == 0, na.rm = TRUE),
               sum(ASDandTD_recoded$Child_Sex        == 0, na.rm = TRUE))
)
print(recode_summary)
write.csv(recode_summary,
          file.path(RESULTS$recoding, "recoding_summary.csv"),
          row.names = FALSE)

saveRDS(ASDandTD_recoded,
        file.path(RESULTS$recoding, "ASDandTD_recoded.rds"))
cat("Saved: ASDandTD_recoded.rds —", nrow(ASDandTD_recoded), "rows\n")