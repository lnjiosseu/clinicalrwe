
# =============================================================================
# ClinicalRWE: 02_disparity_analysis.R
# HEOR access disparities
# =============================================================================

library(tidyverse)

cohort <- read_csv("data/cohort.csv", show_col_types = FALSE)
claims <- read_csv("data/claims.csv", show_col_types = FALSE)

therapy_disp <- cohort %>%
  group_by(race_ethnicity, insurance_type) %>%
  summarise(
    n = n(),
    therapy_rate = round(100 * mean(preventive_therapy), 2),
    progression_rate = round(100 * mean(disease_progression), 2),
    .groups = "drop"
  )

write_csv(therapy_disp, "outputs/therapy_disparities.csv")

therapy_model <- glm(
  preventive_therapy ~ race_ethnicity +
    insurance_type + poverty_pct +
    age + ecog_score,
  data = cohort,
  family = binomial()
)

coef_tbl <- tibble(
  variable = names(coef(therapy_model)),
  odds_ratio = round(exp(coef(therapy_model)), 3)
)

write_csv(coef_tbl, "outputs/therapy_access_model.csv")

p1 <- therapy_disp %>%
  ggplot(aes(x = race_ethnicity,
             y = therapy_rate,
             fill = insurance_type)) +
  geom_col(position = "dodge") +
  theme_minimal(base_size = 13)

ggsave("plots/therapy_disparities.png", p1, width = 7, height = 4)

cost_tbl <- claims %>%
  left_join(cohort %>% select(patient_id, insurance_type),
            by = "patient_id") %>%
  group_by(insurance_type) %>%
  summarise(mean_cost = round(mean(total_cost), 2),
            .groups = "drop")

write_csv(cost_tbl, "outputs/cost_disparities.csv")

p2 <- cost_tbl %>%
  ggplot(aes(x = insurance_type,
             y = mean_cost,
             fill = insurance_type)) +
  geom_col() +
  theme_minimal(base_size = 13)

ggsave("plots/cost_disparities.png", p2, width = 6, height = 4)
