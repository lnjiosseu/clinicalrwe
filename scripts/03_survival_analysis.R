
# =============================================================================
# ClinicalRWE: 03_survival_analysis.R
# Survival analysis
# =============================================================================

library(tidyverse)
library(survival)
library(broom)

cohort <- read_csv("data/cohort.csv", show_col_types = FALSE)

surv_obj <- Surv(cohort$observed_time,
                 cohort$disease_progression)

km_therapy <- survfit(surv_obj ~ preventive_therapy,
                      data = cohort)

km_tbl <- broom::tidy(km_therapy)

write_csv(km_tbl, "outputs/km_intervention.csv")

cox_model <- coxph(
  Surv(observed_time, disease_progression) ~
    preventive_therapy + age + ecog_score +
    biomarker_positive + charlson_index + bmi,
  data = cohort
)

cox_sum <- summary(cox_model)

cox_tbl <- tibble(
  variable = rownames(cox_sum$coefficients),
  hr = round(exp(coef(cox_model)), 3),
  hr_lower95 = round(cox_sum$conf.int[, "lower .95"], 3),
  hr_upper95 = round(cox_sum$conf.int[, "upper .95"], 3),
  p_value = round(cox_sum$coefficients[, "Pr(>|z|)"], 4)
)

write_csv(cox_tbl, "outputs/cox_results.csv")

ph_test <- cox.zph(cox_model)

ph_tbl <- tibble(
  variable = rownames(ph_test$table),
  p_value = round(ph_test$table[, "p"], 4)
)

write_csv(ph_tbl, "outputs/ph_test.csv")

p1 <- km_tbl %>%
  ggplot(aes(x = time,
             y = estimate,
             color = strata)) +
  geom_step(linewidth = 1) +
  theme_minimal(base_size = 13)

ggsave("plots/km_therapy.png", p1, width = 7, height = 5)

p2 <- cox_tbl %>%
  mutate(variable = forcats::fct_reorder(variable, hr)) %>%
  ggplot(aes(x = variable,
             y = hr,
             ymin = hr_lower95,
             ymax = hr_upper95)) +
  geom_pointrange(color = "#2980b9") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  theme_minimal(base_size = 13)

ggsave("plots/cox_forest.png", p2, width = 7, height = 5)
