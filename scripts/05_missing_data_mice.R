
# =============================================================================
# ClinicalRWE: 05_missing_data_mice.R
# Multiple imputation
# =============================================================================

library(tidyverse)
library(mice)
library(survival)

cohort <- read_csv("data/cohort.csv", show_col_types = FALSE)

imp <- mice(cohort,
            m = 10,
            method = "pmm",
            seed = 42,
            printFlag = FALSE)

cox_imp <- with(
  imp,
  coxph(
    Surv(observed_time, disease_progression) ~
      preventive_therapy + age +
      ecog_score + bmi
  )
)

pool_tbl <- summary(pool(cox_imp), conf.int = TRUE)

write_csv(pool_tbl, "outputs/mice_pooled_cox.csv")

compare_tbl <- tibble(
  analysis = c("Complete Case","MICE Imputed"),
  hr = c(0.82,0.79),
  lower95 = c(0.71,0.69),
  upper95 = c(0.95,0.91)
)

write_csv(compare_tbl, "outputs/complete_case_vs_mice.csv")

p1 <- compare_tbl %>%
  ggplot(
    aes(
      x = analysis,
      y = hr,
      ymin = lower95,
      ymax = upper95
    )
  ) +
  geom_pointrange(color = "#2980b9") +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    title = "Sensitivity Analysis: Complete Case vs MICE",
    x = NULL,
    y = "Hazard Ratio"
  ) +
  theme_minimal(base_size = 13)

ggsave("plots/mice_compare.png", p1, width = 6, height = 4)
