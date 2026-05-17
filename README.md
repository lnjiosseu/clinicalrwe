
# ClinicalRWE: Oncology Survival, Causal Inference & Competing Risks Pipeline

An end-to-end pharmaceutical real-world evidence (RWE) pipeline on a simulated longitudinal oncology cohort of 8,000 patients — built in R, covering HEOR disparity analysis, Kaplan-Meier and Cox PH survival analysis, IPTW / doubly robust causal inference (AIPW), MICE multiple imputation, and Fine-Gray competing risks modeling.

Directly mirrors workflows used in oncology RWE, HEOR, and regulatory-grade post-market effectiveness studies across pharmaceutical and biotech organizations.

---

## Research Question

Does receipt of a preventive maintenance biologic therapy delay time to disease progression, and does this effect persist after rigorously adjusting for confounding, treatment-selection bias, missing data, and competing discontinuation events?

---

## Why This Project Exists

Pharmaceutical real-world evidence requires a methodological toolkit that goes beyond standard ML workflows: longitudinal EHR/claims linkage, comparative effectiveness research, time-to-event modeling, competing risks methodology, causal inference under non-random treatment assignment, and regulatory-style reproducibility.

This project builds that full stack — from synthetic cohort generation through interactive dashboards and reproducible reporting — in the same methodological register as FDA-aligned oncology observational studies and HEOR evidence-generation pipelines.

---

## What It Does

**Module 1 — Clinical Cohort Construction**
- 8,000-patient longitudinal oncology cohort with demographics, insurance, socioeconomic variables, smoking history, BMI, and follow-up time
- ECOG performance status, biomarker positivity, and line-of-therapy assignment
- Preventive biologic therapy assignment with realistic confounding structures
- Weibull-distributed progression outcomes with known hazard multipliers
- Competing treatment discontinuation events
- Claims-style utilization and cost dataset generation

**Module 2 — HEOR & Access Disparities**
- Therapy-access disparities by race/ethnicity and insurance type
- Comparative progression burden across demographic groups
- Logistic regression for treatment-access inequities
- Cost disparity summaries across payer categories

**Module 3 — Survival Analysis**
- Kaplan-Meier progression-free survival curves
- Log-rank tests for treatment-group comparison
- Cox proportional hazards model with multivariable adjustment
- Schoenfeld residual diagnostics for PH assumption testing
- Forest plot of hazard ratios with 95% confidence intervals

**Module 4 — Causal Inference**
- Propensity score estimation with positivity assessment
- Stabilized IPTW weighting
- Standardized mean difference (SMD) balance assessment
- IPTW-weighted average treatment effect (ATE)
- Doubly robust AIPW estimator
- E-value sensitivity analysis for unmeasured confounding

**Module 5 — Missing Data**
- Structured MCAR/MAR missingness generation
- Multiple imputation using MICE (m = 10)
- Rubin pooled inference
- Complete-case vs imputed estimate comparison

**Module 6 — Competing Risks**
- Cause-specific Cox regression
- Fine-Gray subdistribution hazard modeling
- Cumulative incidence functions (CIF)
- Disease progression vs treatment discontinuation comparison

**Deliverables**
- Interactive Shiny dashboard (7 tabs: Overview, HEOR Disparities, Survival Analysis, Causal Inference, Missing Data, Competing Risks, Patient Explorer)
- Quarto HTML report with reproducible code, figures, diagnostics, and methods narrative

---

## Key Results

| Module | Finding |
|---|---|
| HEOR disparities | Medicaid and uninsured patients have lower preventive-therapy uptake |
| Survival | Preventive biologic therapy associated with longer progression-free survival |
| Cox model | ECOG score, smoking, and Charlson index are strongest progression predictors |
| Causal (IPTW) | Negative ATE: therapy reduces progression probability |
| Causal (AIPW) | Doubly robust estimate consistent with IPTW findings |
| Missing data | MICE pooled estimates remain directionally consistent with complete-case analysis |
| Competing risks | Fine-Gray estimates differ from cause-specific hazards after accounting for discontinuation |
| Balance | All major covariates achieve \|SMD\| < 0.10 after weighting |

---

## Methods

| Component | Method |
|---|---|
| HEOR profiling | Stratified rates, logistic regression, odds ratios |
| Survival analysis | Kaplan-Meier estimator, log-rank test |
| Multivariable survival | Cox proportional hazards model (Efron tie-handling) |
| PH diagnostics | Schoenfeld residuals (cox.zph) |
| Competing risks | Fine-Gray subdistribution hazards |
| Missing data | MICE multiple imputation |
| Pooling | Rubin's Rules |
| Causal identification | Conditional exchangeability + positivity |
| Propensity score | Logistic regression |
| Weighting | Stabilized IPTW |
| Balance assessment | Standardized mean differences (SMD) |
| Causal estimator 1 | IPTW-weighted outcome model |
| Causal estimator 2 | Augmented IPW (AIPW) — doubly robust |
| Sensitivity analysis | E-value (VanderWeele & Ding, 2017) |
| Visualization | ggplot2, Shiny, bslib |
| Standards alignment | ICH E9(R1), FDA RWE guidance |

---

## Project Structure

```text
clinicalrwe/
├── scripts/
│   ├── 01_simulate_data.R
│   ├── 02_disparity_analysis.R
│   ├── 03_survival_analysis.R
│   ├── 04_causal_inference.R
│   ├── 05_missing_data_mice.R
│   └── 06_competing_risks.R
├── data/                         # Generated datasets (git-ignored)
├── outputs/                      # CSV/model artifacts (git-ignored)
├── plots/                        # Publication/dashboard figures
├── shiny/
│   └── app.R
├── pubhealth_report.qmd
├── README.md
└── clinicalrwe.Rproj
```

---

## Reproducing the Project

```r
# Install dependencies
install.packages(c(
  "tidyverse","survival","cmprsk","mice",
  "shiny","bslib","bsicons","DT","broom"
))

# Run scripts sequentially
source("scripts/01_simulate_data.R")
source("scripts/02_disparity_analysis.R")
source("scripts/03_survival_analysis.R")
source("scripts/04_causal_inference.R")
source("scripts/05_missing_data_mice.R")
source("scripts/06_competing_risks.R")

# Launch dashboard
shiny::runApp("shiny/")

# Render Quarto report
quarto::quarto_render("pubhealth_report.qmd")
```

---

## Connection to Real-World Work

The methods in this project directly mirror workflows used in:
- oncology RWE
- HEOR comparative effectiveness research
- observational biostatistics
- pharmacoepidemiology
- post-market evidence generation

and align closely with my prior work across:

- **Pfizer RWE (2021–2023):** End-to-end RWE study design, patient-level EHR/claims analysis, comparative effectiveness workflows, and causal inference across multiple therapeutic areas
- **FDNY/Montefiore (2025–2026):** Survival analysis, IPTW pipelines, longitudinal cohort analytics, and regulatory-grade reproducibility under ICH E9(R1)
- **NYU Biostatistics (2019–2020):** Claims-based utilization analysis, regression modeling, and ETL pipeline development on healthcare datasets

---

## Author

**Ludovic Njiosseu** | MS Biostatistics, NYU | Data Scientist & Statistician  
5+ years across healthcare, pharma, public health, and consumer products  
Open to roles in RWE, HEOR, biostatistics, epidemiology, and healthcare data science
