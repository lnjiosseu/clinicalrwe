
# =============================================================================
# ClinicalRWE: 01_simulate_data.R
# Simulated oncology RWE cohort generation
# =============================================================================

library(tidyverse)

dir.create("data", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

set.seed(42)

N <- 8000
max_followup <- 48

cohort <- tibble(
  patient_id = paste0("PT", str_pad(1:N, 5, pad = "0")),
  age = round(rnorm(N, 64, 10)),
  sex = sample(c("Male","Female"), N, TRUE),
  race_ethnicity = sample(c("White","Black","Hispanic","Asian","Other"),
                          N, TRUE, prob = c(0.58,0.16,0.14,0.08,0.04)),
  insurance_type = sample(c("Commercial","Medicare","Medicaid"),
                           N, TRUE, prob = c(0.42,0.44,0.14)),
  ecog_score = sample(0:3, N, TRUE, prob = c(0.35,0.40,0.20,0.05)),
  biomarker_positive = rbinom(N, 1, 0.38),
  line_of_therapy = sample(c("1L","2L","3L+"),
                            N, TRUE, prob = c(0.58,0.29,0.13)),
  smoking_status = sample(c("Never","Former","Current"),
                           N, TRUE, prob = c(0.52,0.33,0.15)),
  bmi = round(rnorm(N, 28, 5), 1),
  poverty_pct = round(runif(N, 0, 1), 3),
  charlson_index = rpois(N, 2)
)

therapy_logit <- with(cohort,
  -0.8 +
  0.6 * biomarker_positive -
  0.25 * ecog_score -
  0.20 * (line_of_therapy == "3L+") -
  0.12 * poverty_pct
)

cohort$preventive_therapy <- rbinom(N, 1, plogis(therapy_logit))

prog_hazard <- with(cohort,
  exp(
    -0.42 * preventive_therapy +
    0.20 * ecog_score +
    0.08 * charlson_index +
    0.15 * (smoking_status == "Current")
  )
)

tox_hazard <- with(cohort,
  exp(
    0.25 * preventive_therapy +
    0.18 * ecog_score
  )
)

progression_time <- rweibull(N, shape = 1.3,
                             scale = 30 / prog_hazard)

discontinuation_time <- rweibull(N, shape = 1.1,
                                 scale = 42 / tox_hazard)

cohort <- cohort %>%
  mutate(
    observed_time = pmin(progression_time,
                         discontinuation_time,
                         max_followup),

    event_type = case_when(
      progression_time <= discontinuation_time &
      progression_time <= max_followup ~ 1,

      discontinuation_time < progression_time &
      discontinuation_time <= max_followup ~ 2,

      TRUE ~ 0
    ),

    disease_progression = as.integer(event_type == 1),
    treatment_discontinuation = as.integer(event_type == 2)
  )

# Missingness
set.seed(99)

cohort <- cohort %>%
  mutate(
    biomarker_positive = ifelse(runif(N) < 0.10, NA, biomarker_positive),
    bmi = ifelse(runif(N) < 0.08, NA, bmi),
    charlson_index = ifelse(runif(N) < 0.07, NA, charlson_index),
    poverty_pct = ifelse(runif(N) < 0.06, NA, poverty_pct)
  )

claims <- cohort %>%
  transmute(
    patient_id,
    ed_visits = rpois(N, 1.6 + 0.3 * ecog_score),
    outpatient_visits = rpois(N, 5 + preventive_therapy),
    total_cost = round(rlnorm(N, 9.5, 0.6), 2)
  )

write_csv(cohort, "data/cohort.csv")
write_csv(claims, "data/claims.csv")

p1 <- cohort %>%
  count(line_of_therapy) %>%
  ggplot(aes(x = line_of_therapy, y = n, fill = line_of_therapy)) +
  geom_col() +
  theme_minimal(base_size = 13)

ggsave("plots/line_of_therapy.png", p1, width = 6, height = 4)

p2 <- cohort %>%
  count(ecog_score) %>%
  ggplot(aes(x = factor(ecog_score), y = n, fill = factor(ecog_score))) +
  geom_col() +
  theme_minimal(base_size = 13)

ggsave("plots/ecog_distribution.png", p2, width = 6, height = 4)
