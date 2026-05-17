
# =============================================================================
# ClinicalRWE: 08_subgroup_analysis.R
# Biomarker-stratified subgroup analysis + interaction test
# Standard in oncology RWE — biomarker positivity is a key effect modifier
# =============================================================================

library(tidyverse)
library(survival)
library(broom)

cohort <- read_csv("data/cohort.csv", show_col_types = FALSE)

# -----------------------------------------------------------------------------
# 1. Biomarker-stratified Kaplan-Meier curves
# -----------------------------------------------------------------------------
km_biomarker <- survfit(
  Surv(observed_time, disease_progression) ~
    preventive_therapy + biomarker_positive,
  data = cohort %>% filter(!is.na(biomarker_positive))
)

km_biomarker_tbl <- broom::tidy(km_biomarker) %>%
  mutate(
    therapy   = ifelse(grepl("preventive_therapy=1", strata), "Therapy", "Control"),
    biomarker = ifelse(grepl("biomarker_positive=1", strata),
                       "Biomarker+", "Biomarker-")
  )

write_csv(km_biomarker_tbl, "outputs/km_biomarker_subgroup.csv")

# Median PFS by subgroup
median_pfs <- summary(km_biomarker)$table %>%
  as.data.frame() %>%
  rownames_to_column("strata") %>%
  select(strata, median) %>%
  mutate(
    therapy   = ifelse(grepl("preventive_therapy=1", strata), "Therapy", "Control"),
    biomarker = ifelse(grepl("biomarker_positive=1", strata), "Biomarker+", "Biomarker-")
  )

cat("--- Median PFS by Subgroup ---\n")
print(median_pfs %>% select(therapy, biomarker, median))

# Log-rank tests within each biomarker stratum
lr_bm_pos <- survdiff(
  Surv(observed_time, disease_progression) ~ preventive_therapy,
  data = cohort %>% filter(biomarker_positive == 1)
)
lr_bm_neg <- survdiff(
  Surv(observed_time, disease_progression) ~ preventive_therapy,
  data = cohort %>% filter(biomarker_positive == 0)
)

lr_p_pos <- round(pchisq(lr_bm_pos$chisq, df = 1, lower.tail = FALSE), 4)
lr_p_neg <- round(pchisq(lr_bm_neg$chisq, df = 1, lower.tail = FALSE), 4)

cat(sprintf("\nLog-rank p — Biomarker+: %.4f | Biomarker-: %.4f\n",
            lr_p_pos, lr_p_neg))

# -----------------------------------------------------------------------------
# 2. Subgroup-specific Cox HRs
# -----------------------------------------------------------------------------
cox_bm_pos <- coxph(
  Surv(observed_time, disease_progression) ~
    preventive_therapy + age + ecog_score + charlson_index + bmi,
  data = cohort %>% filter(biomarker_positive == 1)
)

cox_bm_neg <- coxph(
  Surv(observed_time, disease_progression) ~
    preventive_therapy + age + ecog_score + charlson_index + bmi,
  data = cohort %>% filter(biomarker_positive == 0)
)

subgroup_hr <- tibble(
  subgroup   = c("Biomarker Positive", "Biomarker Negative", "Overall"),
  n          = c(
    sum(cohort$biomarker_positive == 1, na.rm = TRUE),
    sum(cohort$biomarker_positive == 0, na.rm = TRUE),
    nrow(cohort)
  ),
  hr         = round(c(
    exp(coef(cox_bm_pos)["preventive_therapy"]),
    exp(coef(cox_bm_neg)["preventive_therapy"]),
    exp(coef(coxph(Surv(observed_time, disease_progression) ~
                     preventive_therapy + age + ecog_score +
                     charlson_index + bmi,
                   data = cohort))["preventive_therapy"])
  ), 3),
  lower95    = round(c(
    exp(confint(cox_bm_pos)["preventive_therapy", 1]),
    exp(confint(cox_bm_neg)["preventive_therapy", 1]),
    exp(confint(coxph(Surv(observed_time, disease_progression) ~
                        preventive_therapy + age + ecog_score +
                        charlson_index + bmi,
                      data = cohort))["preventive_therapy", 1])
  ), 3),
  upper95    = round(c(
    exp(confint(cox_bm_pos)["preventive_therapy", 2]),
    exp(confint(cox_bm_neg)["preventive_therapy", 2]),
    exp(confint(coxph(Surv(observed_time, disease_progression) ~
                        preventive_therapy + age + ecog_score +
                        charlson_index + bmi,
                      data = cohort))["preventive_therapy", 2])
  ), 3),
  logrank_p  = c(lr_p_pos, lr_p_neg, NA_real_)
)

cat("\n--- Subgroup HRs ---\n")
print(subgroup_hr)
write_csv(subgroup_hr, "outputs/subgroup_hr.csv")

# -----------------------------------------------------------------------------
# 3. Interaction test: biomarker × therapy
#    Tests whether biomarker status meaningfully modifies the therapy effect
# -----------------------------------------------------------------------------
cox_interaction <- coxph(
  Surv(observed_time, disease_progression) ~
    preventive_therapy * biomarker_positive +
    age + ecog_score + charlson_index + bmi,
  data = cohort %>% filter(!is.na(biomarker_positive))
)

interaction_p <- round(
  summary(cox_interaction)$coefficients[
    "preventive_therapy:biomarker_positive", "Pr(>|z|)"
  ], 4
)

interaction_hr <- round(
  exp(coef(cox_interaction)["preventive_therapy:biomarker_positive"]), 3
)

cat(sprintf(
  "\nInteraction test — biomarker × therapy: HR = %.3f, p = %.4f\n",
  interaction_hr, interaction_p
))

if (interaction_p < 0.05) {
  cat("Significant interaction: biomarker status modifies the therapy effect.\n")
  cat("Biomarker-positive patients derive differential benefit.\n")
} else {
  cat("No significant interaction detected (p >= 0.05).\n")
  cat("Subgroup differences may reflect chance variation.\n")
}

write_csv(
  tibble(
    interaction_hr = interaction_hr,
    interaction_p  = interaction_p,
    significant    = interaction_p < 0.05
  ),
  "outputs/interaction_test.csv"
)

# -----------------------------------------------------------------------------
# 4. Forest plot — subgroup HRs
# -----------------------------------------------------------------------------
p1 <- subgroup_hr %>%
  mutate(subgroup = fct_rev(fct_inorder(subgroup))) %>%
  ggplot(aes(x = subgroup, y = hr,
             ymin = lower95, ymax = upper95,
             color = subgroup)) +
  geom_pointrange(size = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  coord_flip() +
  scale_color_manual(values = c(
    "Biomarker Positive" = "#27ae60",
    "Biomarker Negative" = "#e74c3c",
    "Overall"            = "#2980b9"
  )) +
  annotate("text", x = 0.6, y = max(subgroup_hr$upper95, na.rm = TRUE),
           label = sprintf("Interaction p = %.4f", interaction_p),
           hjust = 1, size = 3.5, color = "grey30") +
  labs(title = "Subgroup Analysis: Therapy Effect by Biomarker Status",
       x = NULL, y = "Hazard Ratio (95% CI)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave("plots/subgroup_forest.png", p1, width = 7, height = 4)

# KM curves faceted by biomarker status
p2 <- km_biomarker_tbl %>%
  ggplot(aes(x = time, y = estimate, color = therapy)) +
  geom_step(linewidth = 1) +
  facet_wrap(~ biomarker, labeller = label_both) +
  scale_color_manual(values = c(Control = "#e74c3c", Therapy = "#2980b9")) +
  labs(title = "PFS by Biomarker Status and Treatment",
       x = "Months", y = "Progression-Free Survival", color = NULL) +
  theme_minimal(base_size = 13)

ggsave("plots/km_biomarker_subgroup.png", p2, width = 9, height = 5)

cat("\nSubgroup analysis complete.\n")
