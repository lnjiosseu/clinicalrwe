
# =============================================================================
# ClinicalRWE: 07_table1_baseline.R
# Baseline Characteristics Table
# =============================================================================

library(tidyverse)

cohort <- read_csv("data/cohort.csv", show_col_types = FALSE)

table1_tbl <- cohort %>%
  mutate(
    therapy_group = ifelse(
      preventive_therapy == 1,
      "Preventive Therapy",
      "Control"
    )
  ) %>%
  group_by(therapy_group) %>%
  summarise(
    N = n(),
    Mean_Age = round(mean(age, na.rm = TRUE), 1),
    Female_Pct = round(100 * mean(sex == "Female"), 1),
    Mean_BMI = round(mean(bmi, na.rm = TRUE), 1),
    ECOG_High_Pct = round(100 * mean(ecog_score >= 2), 1),
    Biomarker_Positive_Pct = round(
      100 * mean(biomarker_positive == 1, na.rm = TRUE),
      1
    ),
    Charlson_Mean = round(
      mean(charlson_index, na.rm = TRUE),
      2
    ),
    Current_Smoker_Pct = round(
      100 * mean(smoking_status == "Current"),
      1
    ),
    Medicaid_Pct = round(
      100 * mean(insurance_type == "Medicaid"),
      1
    ),
    Medicare_Pct = round(
      100 * mean(insurance_type == "Medicare"),
      1
    ),
    .groups = "drop"
  )

write_csv(
  table1_tbl,
  "outputs/table1_baseline.csv"
)

print(table1_tbl)
