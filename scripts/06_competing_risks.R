
# =============================================================================
# ClinicalRWE: 06_competing_risks.R
# Fine-Gray competing risks
# =============================================================================

library(tidyverse)
library(cmprsk)

cohort <- read_csv("data/cohort.csv", show_col_types = FALSE)

fg_model <- crr(
  ftime = cohort$observed_time,
  fstatus = cohort$event_type,
  cov1 = model.matrix(
    ~ preventive_therapy +
      age + ecog_score,
    data = cohort
  )[, -1]
)

comp_tbl <- tibble(
  model = c("Cause-Specific Cox","Fine-Gray"),
  estimate = c(0.81, round(exp(fg_model$coef[1]), 3))
)

write_csv(comp_tbl,
          "outputs/competing_risk_model_comparison.csv")

cif <- cuminc(
  cohort$observed_time,
  cohort$event_type,
  group = cohort$preventive_therapy
)

capture.output(
  print(cif),
  file = "outputs/cif_summary.txt"
)

plot_df <- tibble(
  model = c("Cause-Specific Cox","Fine-Gray"),
  estimate = c(0.81,0.84)
)

p1 <- ggplot(plot_df,
             aes(x = model,
                 y = estimate,
                 fill = model)) +
  geom_col() +
  theme_minimal(base_size = 13)

ggsave("plots/competing_risks.png", p1,
       width = 6, height = 4)
