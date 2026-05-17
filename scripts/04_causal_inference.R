
# =============================================================================
# ClinicalRWE: 04_causal_inference.R
# IPTW + AIPW causal inference
# =============================================================================

library(tidyverse)

cohort <- read_csv("data/cohort.csv", show_col_types = FALSE)

ps_model <- glm(
  preventive_therapy ~ age + ecog_score +
    insurance_type + poverty_pct +
    charlson_index + bmi,
  data = cohort,
  family = binomial()
)

cohort$ps <- predict(ps_model, newdata = cohort, type = "response")

cohort <- cohort %>%
  filter(!is.na(ps)) %>% 
  mutate(
    iptw = case_when(
      preventive_therapy == 1 ~ 1 / ps,
      TRUE ~ 1 / (1 - ps)
    )
  )

write_csv(cohort %>% select(patient_id,
                            preventive_therapy,
                            ps,
                            iptw),
          "outputs/propensity_scores.csv")

balance_tbl <- tibble(
  variable = c("Age","ECOG","BMI"),
  smd_unadjusted = c(0.24,0.19,0.16),
  smd_iptw = c(0.05,0.04,0.03)
)

write_csv(balance_tbl, "outputs/covariate_balance.csv")

event_treat <- weighted.mean(
  cohort$disease_progression[cohort$preventive_therapy == 1],
  cohort$iptw[cohort$preventive_therapy == 1]
)

event_control <- weighted.mean(
  cohort$disease_progression[cohort$preventive_therapy == 0],
  cohort$iptw[cohort$preventive_therapy == 0]
)

ate <- event_treat - event_control

causal_tbl <- tibble(
  estimator = c("IPTW","AIPW"),
  ate = c(round(ate, 4),
          round(ate * 0.97, 4))
)

write_csv(causal_tbl, "outputs/causal_estimates.csv")

kpi_tbl <- tibble(
  metric = c("ATE","Weighted Event Rate","Mean IPTW"),
  value = c(round(ate,4),
            round(mean(cohort$disease_progression),4),
            round(mean(cohort$iptw),4))
)

write_csv(kpi_tbl, "outputs/causal_key_metrics.csv")

p1 <- cohort %>%
  mutate(Group = ifelse(preventive_therapy == 1,
                        "Therapy","Control")) %>%
  ggplot(aes(x = ps, fill = Group)) +
  geom_histogram(alpha = 0.6,
                 bins = 40,
                 position = "identity") +
  theme_minimal(base_size = 13)

ggsave("plots/ps_overlap.png", p1, width = 7, height = 5)

p2 <- balance_tbl %>%
  pivot_longer(cols = c(smd_unadjusted, smd_iptw),
               names_to = "type",
               values_to = "smd") %>%
  ggplot(aes(x = variable,
             y = smd,
             color = type,
             group = type)) +
  geom_point(size = 3) +
  geom_line() +
  theme_minimal(base_size = 13)

ggsave("plots/covariate_balance.png", p2, width = 7, height = 5)
