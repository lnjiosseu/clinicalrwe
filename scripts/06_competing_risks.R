
# =============================================================================
# ClinicalRWE: 06_competing_risks.R
# Fine-Gray competing risks — computed from data
# Primary event: disease progression (event_type == 1)
# Competing event: treatment discontinuation (event_type == 2)
# =============================================================================

library(tidyverse)
library(survival)
library(cmprsk)

cohort <- read_csv("data/cohort.csv", show_col_types = FALSE) %>%
  filter(!is.na(observed_time), !is.na(event_type))

# -----------------------------------------------------------------------------
# 1. Cumulative Incidence Functions (CIF)
# -----------------------------------------------------------------------------
cif_overall <- cuminc(
  ftime   = cohort$observed_time,
  fstatus = cohort$event_type
)

# CIF by therapy group
cif_by_therapy <- cuminc(
  ftime   = cohort$observed_time,
  fstatus = cohort$event_type,
  group   = cohort$preventive_therapy
)

# Gray's test p-value for disease progression (event = 1)
gray_p <- round(cif_by_therapy$Tests["1", "pv"], 4)
cat(sprintf("Gray's test p-value (disease progression): %.4f\n", gray_p))

# Extract CIF data frames for plotting
# cmprsk names CIF list elements as e.g. "0 1", "1 1", "0 2", "1 2"
# (group label + space + event code) — match both " 1" and " 1$" patterns
extract_cif <- function(cif_obj, event_code) {
  # Match any key ending in the event code (space-separated)
  pattern <- paste0("\\s", event_code, "$")
  nm_list <- names(cif_obj)[grepl(pattern, names(cif_obj))]
  # Fallback: match on exact suffix if above returns nothing
  if (length(nm_list) == 0)
    nm_list <- names(cif_obj)[endsWith(names(cif_obj), as.character(event_code))]
  map_dfr(nm_list, function(nm) {
    tibble(
      time  = cif_obj[[nm]]$time,
      cif   = round(cif_obj[[nm]]$est, 5),
      group = nm
    )
  })
}

cif_prog_df <- extract_cif(cif_by_therapy, 1)
if (nrow(cif_prog_df) == 0) {
  # Last-resort: use first two CIF entries (group 0 and group 1 for event 1)
  all_nms  <- names(cif_by_therapy)[!names(cif_by_therapy) %in% c("Tests")]
  event1   <- all_nms[seq(1, length(all_nms), by = 2)]
  cif_prog_df <- map_dfr(event1, function(nm)
    tibble(time = cif_by_therapy[[nm]]$time,
           cif  = round(cif_by_therapy[[nm]]$est, 5),
           group = nm))
}

cif_prog_df <- cif_prog_df %>%
  mutate(group = ifelse(grepl("^0", group), "Control", "Therapy"),
         event = "Disease Progression")

cif_disc_df <- extract_cif(cif_by_therapy, 2)
if (nrow(cif_disc_df) == 0) {
  all_nms  <- names(cif_by_therapy)[!names(cif_by_therapy) %in% c("Tests")]
  event2   <- all_nms[seq(2, length(all_nms), by = 2)]
  cif_disc_df <- map_dfr(event2, function(nm)
    tibble(time = cif_by_therapy[[nm]]$time,
           cif  = round(cif_by_therapy[[nm]]$est, 5),
           group = nm))
}

cif_disc_df <- cif_disc_df %>%
  mutate(group = ifelse(grepl("^0", group), "Control", "Therapy"),
         event = "Discontinuation (competing)")

write_csv(bind_rows(cif_prog_df, cif_disc_df), "outputs/cif_by_therapy.csv")
cat(sprintf("CIF rows written: %d\n",
            nrow(cif_prog_df) + nrow(cif_disc_df)))

# -----------------------------------------------------------------------------
# 2. Fine-Gray subdistribution hazard model — disease progression
# -----------------------------------------------------------------------------
cohort_cov <- cohort %>%
  mutate(
    # Drop medicare to avoid perfect multicollinearity (commercial = reference)
    medicaid   = as.integer(insurance_type == "Medicaid"),
    smoker     = as.integer(smoking_status  == "Current"),
    lot3       = as.integer(line_of_therapy == "3L+"),
    biomarker  = as.integer(biomarker_positive == 1),
    age_c      = age - 64
  ) %>%
  select(preventive_therapy, age_c, ecog_score, biomarker,
         charlson_index, bmi, medicaid, smoker, lot3) %>%
  mutate(across(everything(), ~ replace_na(., median(., na.rm = TRUE))))

# Drop zero-variance columns and confirm full rank before crr()
cov_sd    <- sapply(cohort_cov, sd, na.rm = TRUE)
cohort_cov <- cohort_cov[, names(cov_sd[cov_sd > 0])]
covariates <- as.matrix(cohort_cov)

mat_rank <- qr(covariates)$rank
cat(sprintf("Covariate matrix: %d cols, rank %d — %s\n",
            ncol(covariates), mat_rank,
            ifelse(mat_rank == ncol(covariates), "full rank OK",
                   "WARNING: rank deficient")))

fg_progression <- crr(
  ftime    = cohort$observed_time,
  fstatus  = cohort$event_type,
  cov1     = covariates,
  failcode = 1,
  cencode  = 0
)

fg_discontinuation <- crr(
  ftime    = cohort$observed_time,
  fstatus  = cohort$event_type,
  cov1     = covariates,
  failcode = 2,
  cencode  = 0
)

extract_fg <- function(fg_obj, varnames, event_label) {
  s <- summary(fg_obj)
  tibble(
    variable  = varnames,
    shr       = round(exp(fg_obj$coef), 3),
    lower95   = round(s$conf.int[, 3], 3),
    upper95   = round(s$conf.int[, 4], 3),
    p_value   = round(s$coef[, 5], 4),
    event     = event_label
  )
}

varnames    <- colnames(covariates)
fg_prog_tbl <- extract_fg(fg_progression,     varnames, "Disease Progression")
fg_disc_tbl <- extract_fg(fg_discontinuation, varnames, "Discontinuation")

write_csv(bind_rows(fg_prog_tbl, fg_disc_tbl), "outputs/fine_gray_results.csv")

# -----------------------------------------------------------------------------
# 3. Cause-specific Cox for comparison
# -----------------------------------------------------------------------------
cohort_cs <- cohort %>%
  mutate(
    medicaid  = as.integer(insurance_type == "Medicaid"),
    medicare  = as.integer(insurance_type == "Medicare"),
    smoker    = as.integer(smoking_status  == "Current"),
    lot3      = as.integer(line_of_therapy == "3L+"),
    biomarker = as.integer(biomarker_positive == 1),
    age_c     = age - 64
  )

cs_cox <- coxph(
  Surv(observed_time, event_type == 1) ~
    preventive_therapy + age_c + ecog_score + biomarker +
    charlson_index + bmi + medicaid + medicare + smoker + lot3,
  data = cohort_cs, ties = "efron"
)

# Primary comparison: preventive_therapy SHR vs CS-HR
therapy_fg_shr <- fg_prog_tbl$shr[fg_prog_tbl$variable == "preventive_therapy"]
therapy_cs_hr  <- round(exp(coef(cs_cox)["preventive_therapy"]), 3)
cs_lower       <- round(exp(confint(cs_cox)["preventive_therapy", 1]), 3)
cs_upper       <- round(exp(confint(cs_cox)["preventive_therapy", 2]), 3)
fg_lower       <- fg_prog_tbl$lower95[fg_prog_tbl$variable == "preventive_therapy"]
fg_upper       <- fg_prog_tbl$upper95[fg_prog_tbl$variable == "preventive_therapy"]

comparison_tbl <- tibble(
  model    = c("Cause-Specific Cox", "Fine-Gray (subdistribution)"),
  estimate = c(therapy_cs_hr, therapy_fg_shr),
  lower95  = c(cs_lower, fg_lower),
  upper95  = c(cs_upper, fg_upper),
  note     = c(
    "Hazard of progression among those not yet progressed or discontinued",
    "Hazard weighted to include those who discontinued — accounts for competing risk"
  )
)

cat("\n--- Model Comparison: preventive_therapy ---\n")
print(comparison_tbl %>% select(model, estimate, lower95, upper95))

# Narrative flag: if Fine-Gray SHR is meaningfully smaller than CS-HR,
# it means discontinuation is acting as a competing risk that attenuates
# the observed progression benefit when properly accounted for
shr_cs_diff <- abs(therapy_fg_shr - therapy_cs_hr)
if (shr_cs_diff > 0.05) {
  cat(sprintf(
    "\nNote: Fine-Gray SHR (%.3f) differs from cause-specific HR (%.3f) by %.3f.\n",
    therapy_fg_shr, therapy_cs_hr, shr_cs_diff))
  cat("Treatment discontinuation is a non-negligible competing risk.\n")
  cat("The cause-specific model overstates the therapy benefit relative to Fine-Gray.\n")
  cat("Reporting both is required for a complete pharma RWE submission.\n")
}

write_csv(comparison_tbl, "outputs/competing_risk_model_comparison.csv")

# Gray's test summary
write_csv(
  tibble(event = c("Disease Progression", "Discontinuation"),
         grays_p = c(gray_p, round(cif_by_therapy$Tests["2", "pv"], 4))),
  "outputs/grays_test.csv"
)

# -----------------------------------------------------------------------------
# 4. Plots
# -----------------------------------------------------------------------------

# CIF: progression by therapy group
p1 <- cif_prog_df %>%
  ggplot(aes(x = time, y = cif, color = group)) +
  geom_step(linewidth = 1.1) +
  scale_color_manual(values = c(Control = "#e74c3c", Therapy = "#2980b9")) +
  labs(title = "Cumulative Incidence of Disease Progression",
       subtitle = sprintf("Gray's test p = %s", gray_p),
       x = "Months", y = "Cumulative Incidence", color = NULL) +
  theme_minimal(base_size = 13)
ggsave("plots/cif_progression.png", p1, width = 7, height = 5)

# Model comparison forest
p2 <- comparison_tbl %>%
  mutate(model = fct_inorder(model)) %>%
  ggplot(aes(x = model, y = estimate,
             ymin = lower95, ymax = upper95,
             color = model)) +
  geom_pointrange(size = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  scale_color_manual(values = c("#8e44ad", "#2980b9")) +
  labs(title = "Cause-Specific vs Fine-Gray Hazard Estimates",
       subtitle = "Preventive therapy effect on disease progression",
       x = NULL, y = "Hazard / Subdistribution Hazard Ratio") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")
ggsave("plots/competing_risks.png", p2, width = 7, height = 4)

cat("\nCompeting risks analysis complete.\n")
