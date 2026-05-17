
# =============================================================================
# ClinicalRWE: 04_causal_inference.R
# IPTW + AIPW causal inference — with bootstrap CIs and E-value
# =============================================================================

library(tidyverse)
library(survival)

cohort <- read_csv("data/cohort.csv", show_col_types = FALSE)

# -----------------------------------------------------------------------------
# 1. Propensity score model
# -----------------------------------------------------------------------------
ps_model <- glm(
  preventive_therapy ~ age + ecog_score +
    insurance_type + poverty_pct +
    charlson_index + bmi,
  data   = cohort,
  family = binomial()
)

cohort$ps <- predict(ps_model, newdata = cohort, type = "response")

cohort <- cohort %>%
  filter(!is.na(ps)) %>%
  mutate(
    p_treat = mean(preventive_therapy),
    # Stabilized weights — reduces variance vs unstabilized 1/ps
    iptw = case_when(
      preventive_therapy == 1 ~ p_treat / ps,
      TRUE                    ~ (1 - p_treat) / (1 - ps)
    ),
    # Trim at 99th percentile to reduce influence of extreme weights
    iptw = pmin(iptw, quantile(iptw, 0.99))
  )

write_csv(
  cohort %>% select(patient_id, preventive_therapy, ps, iptw),
  "outputs/propensity_scores.csv"
)

# -----------------------------------------------------------------------------
# 2. Covariate balance — computed from data, not hardcoded
# -----------------------------------------------------------------------------
smd_fn <- function(x, treat, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  x  <- as.numeric(x)
  m1 <- sum(x[treat == 1] * w[treat == 1], na.rm = TRUE) /
    sum(w[treat == 1], na.rm = TRUE)
  m0 <- sum(x[treat == 0] * w[treat == 0], na.rm = TRUE) /
    sum(w[treat == 0], na.rm = TRUE)
  v1 <- sum(w[treat == 1] * (x[treat == 1] - m1)^2, na.rm = TRUE) /
    sum(w[treat == 1], na.rm = TRUE)
  v0 <- sum(w[treat == 0] * (x[treat == 0] - m0)^2, na.rm = TRUE) /
    sum(w[treat == 0], na.rm = TRUE)
  abs((m1 - m0) / sqrt((v1 + v0) / 2))
}

# Build balance dataset — impute NAs, include iptw directly
cohort_bal <- cohort %>%
  mutate(
    medicaid           = as.integer(insurance_type == "Medicaid"),
    medicare           = as.integer(insurance_type == "Medicare"),
    smoker             = as.integer(smoking_status  == "Current"),
    lot3               = as.integer(as.character(line_of_therapy) == "3L+"),
    bmi                = ifelse(is.na(bmi),
                                median(bmi,            na.rm = TRUE), bmi),
    charlson_index     = ifelse(is.na(charlson_index),
                                median(charlson_index, na.rm = TRUE), charlson_index),
    poverty_pct        = ifelse(is.na(poverty_pct),
                                median(poverty_pct,    na.rm = TRUE), poverty_pct),
    biomarker_positive = ifelse(is.na(biomarker_positive),
                                round(mean(biomarker_positive, na.rm = TRUE)),
                                biomarker_positive)
    # iptw already present from cohort %>% mutate(iptw = ...)
  )

balance_vars <- list(
  "Age"                 = "age",
  "ECOG Score"          = "ecog_score",
  "BMI"                 = "bmi",
  "Charlson Index"      = "charlson_index",
  "Poverty %"           = "poverty_pct",
  "Biomarker Positive"  = "biomarker_positive",
  "Line of Therapy 3L+" = "lot3",
  "Current Smoker"      = "smoker",
  "Medicaid"            = "medicaid",
  "Medicare"            = "medicare"
)

balance_tbl <- tibble(
  variable       = names(balance_vars),
  smd_unadjusted = sapply(balance_vars, function(v)
    round(smd_fn(as.numeric(cohort_bal[[v]]),
                 cohort_bal$preventive_therapy), 3)),
  smd_iptw       = sapply(balance_vars, function(v)
    round(smd_fn(as.numeric(cohort_bal[[v]]),
                 cohort_bal$preventive_therapy,
                 cohort_bal$iptw), 3))
)

write_csv(balance_tbl, "outputs/covariate_balance.csv")

cat(sprintf("Variables balanced (|SMD| < 0.10 after IPTW): %d / %d\n",
            sum(balance_tbl$smd_iptw < 0.10), nrow(balance_tbl)))

# -----------------------------------------------------------------------------
# 3. IPTW-weighted ATE with bootstrap 95% CI
# -----------------------------------------------------------------------------
iptw_ate_fn <- function(df) {
  wt <- sum(df$iptw[df$preventive_therapy == 1] *
              df$disease_progression[df$preventive_therapy == 1]) /
    sum(df$iptw[df$preventive_therapy == 1])
  wc <- sum(df$iptw[df$preventive_therapy == 0] *
              df$disease_progression[df$preventive_therapy == 0]) /
    sum(df$iptw[df$preventive_therapy == 0])
  wt - wc
}

ate_iptw <- iptw_ate_fn(cohort)

set.seed(42)
boot_iptw <- replicate(500, {
  idx <- sample(nrow(cohort), replace = TRUE)
  iptw_ate_fn(cohort[idx, ])
})
ci_iptw <- quantile(boot_iptw, c(0.025, 0.975))

cat(sprintf("IPTW ATE: %.4f | 95%% CI: [%.4f, %.4f]\n",
            ate_iptw, ci_iptw[1], ci_iptw[2]))

# -----------------------------------------------------------------------------
# 4. Doubly Robust AIPW estimator (real implementation)
# -----------------------------------------------------------------------------
# Impute NAs for outcome model covariates (median imputation)
cohort_om <- cohort %>%
  mutate(
    bmi            = ifelse(is.na(bmi),            median(bmi,            na.rm = TRUE), bmi),
    charlson_index = ifelse(is.na(charlson_index), median(charlson_index, na.rm = TRUE), charlson_index),
    poverty_pct    = ifelse(is.na(poverty_pct),    median(poverty_pct,    na.rm = TRUE), poverty_pct),
    biomarker_positive = ifelse(is.na(biomarker_positive),
                                round(mean(biomarker_positive, na.rm = TRUE)),
                                biomarker_positive)
  )

outcome_model <- glm(
  disease_progression ~ preventive_therapy + age + ecog_score +
    bmi + charlson_index + poverty_pct + biomarker_positive +
    smoking_status + line_of_therapy,
  data   = cohort_om,
  family = binomial()
)

cohort_om$mu1 <- predict(outcome_model,
                         newdata = mutate(cohort_om, preventive_therapy = 1),
                         type = "response")
cohort_om$mu0 <- predict(outcome_model,
                         newdata = mutate(cohort_om, preventive_therapy = 0),
                         type = "response")

# Merge mu1/mu0 back to main cohort for AIPW
cohort <- cohort %>%
  mutate(mu1 = cohort_om$mu1,
         mu0 = cohort_om$mu0)

aipw_scores <- with(cohort, {
  (preventive_therapy / ps) * (disease_progression - mu1) + mu1 -
    ((1 - preventive_therapy) / (1 - ps)) * (disease_progression - mu0) - mu0
})

ate_aipw <- mean(aipw_scores, na.rm = TRUE)
se_aipw  <- sd(aipw_scores,   na.rm = TRUE) / sqrt(sum(!is.na(aipw_scores)))
ci_aipw  <- ate_aipw + c(-1.96, 1.96) * se_aipw

cat(sprintf("AIPW ATE: %.4f | 95%% CI: [%.4f, %.4f]\n",
            ate_aipw, ci_aipw[1], ci_aipw[2]))

# -----------------------------------------------------------------------------
# 5. E-value sensitivity analysis
# -----------------------------------------------------------------------------
base_rate <- mean(cohort$disease_progression, na.rm = TRUE)
# RR for protective effect: use ratio of counterfactual rates
rr_treated  <- mean(cohort$mu1, na.rm = TRUE)
rr_control  <- mean(cohort$mu0, na.rm = TRUE)
rr_est      <- rr_treated / rr_control
# E-value defined for RR >= 1; for protective effects use 1/RR
rr_for_eval <- ifelse(rr_est >= 1, rr_est, 1 / rr_est)
e_value     <- rr_for_eval + sqrt(rr_for_eval * (rr_for_eval - 1))

cat(sprintf("E-value: %.3f\n", e_value))

# -----------------------------------------------------------------------------
# 6. Save outputs
# -----------------------------------------------------------------------------
causal_tbl <- tibble(
  estimator = c("IPTW", "AIPW (Doubly Robust)"),
  ate       = round(c(ate_iptw, ate_aipw), 4),
  ci_lower  = round(c(ci_iptw[1], ci_aipw[1]), 4),
  ci_upper  = round(c(ci_iptw[2], ci_aipw[2]), 4)
)
write_csv(causal_tbl, "outputs/causal_estimates.csv")

kpi_tbl <- tibble(
  metric = c("IPTW ATE", "IPTW 95% CI", "AIPW ATE",
             "AIPW 95% CI", "E-value",
             "Covariates balanced (|SMD|<0.10)"),
  value  = c(round(ate_iptw, 4),
             sprintf("[%.4f, %.4f]", ci_iptw[1], ci_iptw[2]),
             round(ate_aipw, 4),
             sprintf("[%.4f, %.4f]", ci_aipw[1], ci_aipw[2]),
             round(e_value, 3),
             paste0(sum(balance_tbl$smd_iptw < 0.10), "/",
                    nrow(balance_tbl)))
)
write_csv(kpi_tbl, "outputs/causal_key_metrics.csv")

# Plots
p1 <- cohort %>%
  mutate(Group = ifelse(preventive_therapy == 1, "Therapy", "Control")) %>%
  ggplot(aes(x = ps, fill = Group)) +
  geom_histogram(alpha = 0.6, bins = 40, position = "identity") +
  scale_fill_manual(values = c(Therapy = "#2980b9", Control = "#e74c3c")) +
  labs(title = "Propensity Score Overlap",
       x = "Propensity Score", y = "Count", fill = NULL) +
  theme_minimal(base_size = 13)
ggsave("plots/ps_overlap.png", p1, width = 7, height = 5)

p2 <- balance_tbl %>%
  pivot_longer(cols = c(smd_unadjusted, smd_iptw),
               names_to = "type", values_to = "smd") %>%
  mutate(
    type     = recode(type,
                      smd_unadjusted = "Unadjusted",
                      smd_iptw       = "IPTW Weighted"),
    smd      = ifelse(is.na(smd), 0, smd),
    variable = forcats::fct_reorder(variable, smd)
  ) %>%
  ggplot(aes(x = smd, y = variable, color = type, shape = type)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0.10, linetype = "dashed", color = "red") +
  labs(title = "Covariate Balance Before and After IPTW",
       x = "Absolute Standardized Mean Difference",
       y = NULL, color = NULL, shape = NULL) +
  theme_minimal(base_size = 13)
ggsave("plots/covariate_balance.png", p2, width = 7, height = 5)

p3 <- causal_tbl %>%
  mutate(estimator = fct_inorder(estimator)) %>%
  ggplot(aes(x = estimator, y = ate,
             ymin = ci_lower, ymax = ci_upper,
             color = estimator)) +
  geom_pointrange(size = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  scale_color_manual(values = c("#2980b9", "#27ae60")) +
  labs(title = "Causal Estimates: IPTW vs AIPW (with 95% CI)",
       x = NULL, y = "Average Treatment Effect (ATE)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")
ggsave("plots/causal_comparison.png", p3, width = 7, height = 4)
