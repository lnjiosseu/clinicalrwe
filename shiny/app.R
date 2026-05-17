# =============================================================================
# ClinicalRWE: shiny/app.R — Oncology RWE Dashboard
# =============================================================================

library(shiny)
library(bslib)
library(tidyverse)
library(survival)
library(DT)

load_data <- function() {
  list(
    cohort        = read_csv("../data/cohort.csv", show_col_types = FALSE),
    claims        = read_csv("../data/claims.csv", show_col_types = FALSE),
    therapy_disp  = read_csv("../outputs/therapy_disparities.csv", show_col_types = FALSE),
    therapy_model = read_csv("../outputs/therapy_access_model.csv", show_col_types = FALSE),
    cost_disp     = read_csv("../outputs/cost_disparities.csv", show_col_types = FALSE),
    km_interv     = read_csv("../outputs/km_intervention.csv", show_col_types = FALSE),
    cox_results   = read_csv("../outputs/cox_results.csv", show_col_types = FALSE),
    ph_test       = read_csv("../outputs/ph_test.csv", show_col_types = FALSE),
    balance       = read_csv("../outputs/covariate_balance.csv", show_col_types = FALSE),
    ps_scores     = read_csv("../outputs/propensity_scores.csv", show_col_types = FALSE),
    causal        = read_csv("../outputs/causal_estimates.csv", show_col_types = FALSE),
    causal_kpi    = read_csv("../outputs/causal_key_metrics.csv", show_col_types = FALSE),
    mice_compare  = read_csv("../outputs/complete_case_vs_mice.csv", show_col_types = FALSE),
    mice_cox      = read_csv("../outputs/mice_pooled_cox.csv", show_col_types = FALSE),
    comp_risk     = read_csv("../outputs/competing_risk_model_comparison.csv", show_col_types = FALSE)
  )
}

d <- load_data()

PAL <- list(
  treat   = "#2980b9",
  control = "#e74c3c",
  green   = "#27ae60",
  gold    = "#f39c12",
  purple  = "#8e44ad",
  dark    = "#2c3e50",
  grey    = "#7f8c8d"
)

ui <- page_navbar(
  title = "ClinicalRWE — Oncology Analytics",
  theme = bs_theme(bootswatch = "flatly", primary = PAL$dark),
  
  # ── Tab 1: Overview ────────────────────────────────────────────────────────
  nav_panel("Overview",
            layout_columns(
              col_widths = c(3,3,3,3),
              
              value_box("Cohort Size",
                        scales::comma(nrow(d$cohort)),
                        showcase = bsicons::bs_icon("people-fill"),
                        theme = "primary"),
              
              value_box("Progression Rate",
                        paste0(round(100 * mean(d$cohort$disease_progression),1),"%"),
                        showcase = bsicons::bs_icon("activity"),
                        theme = "danger"),
              
              value_box("Therapy Exposure",
                        paste0(round(100 * mean(d$cohort$preventive_therapy),1),"%"),
                        showcase = bsicons::bs_icon("capsule"),
                        theme = "success"),
              
              value_box("Biomarker Positive",
                        paste0(round(100 * mean(d$cohort$biomarker_positive,
                                                na.rm = TRUE),1),"%"),
                        showcase = bsicons::bs_icon("clipboard2-pulse"),
                        theme = "warning")
            ),
            
            layout_columns(
              col_widths = c(6,6),
              
              card(
                card_header("ECOG Distribution"),
                plotOutput("plot_ecog", height = "300px")
              ),
              
              card(
                card_header("Line of Therapy"),
                plotOutput("plot_lot", height = "300px")
              )
            )
  ),
  
  # ── Tab 2: HEOR Disparities ───────────────────────────────────────────────
  nav_panel("HEOR Disparities",
            layout_columns(
              col_widths = c(6,6),
              
              card(
                card_header("Therapy Access Disparities"),
                plotOutput("plot_disp", height = "320px")
              ),
              
              card(
                card_header("Mean Cost by Insurance"),
                plotOutput("plot_cost", height = "320px")
              )
            ),
            
            card(
              card_header("Therapy Access Model"),
              DT::dataTableOutput("tbl_model")
            )
  ),
  
  # ── Tab 3: Survival Analysis ──────────────────────────────────────────────
  nav_panel("Survival Analysis",
            layout_columns(
              col_widths = c(6,6),
              
              card(
                card_header("Kaplan-Meier Curves"),
                plotOutput("plot_km", height = "350px")
              ),
              
              card(
                card_header("Hazard Ratios"),
                plotOutput("plot_cox", height = "350px")
              )
            ),
            
            card(
              card_header("PH Diagnostics"),
              tableOutput("tbl_ph")
            )
  ),
  
  # ── Tab 4: Causal Inference ───────────────────────────────────────────────
  nav_panel("Causal Inference",
            layout_columns(
              col_widths = c(6,6),
              
              card(
                card_header("Propensity Score Overlap"),
                plotOutput("plot_ps", height = "320px")
              ),
              
              card(
                card_header("Covariate Balance"),
                plotOutput("plot_balance", height = "320px")
              )
            ),
            
            card(
              card_header("ATE Estimates"),
              tableOutput("tbl_causal")
            )
  ),
  
  # ── Tab 5: Missing Data ───────────────────────────────────────────────────
  nav_panel("Missing Data",
            layout_columns(
              col_widths = c(6,6),
              
              card(
                card_header("Complete Case vs MICE"),
                plotOutput("plot_mice", height = "320px")
              ),
              
              card(
                card_header("Pooled Cox Results"),
                DT::dataTableOutput("tbl_mice")
              )
            )
  ),
  
  # ── Tab 6: Competing Risks ────────────────────────────────────────────────
  nav_panel("Competing Risks",
            layout_columns(
              col_widths = c(6,6),
              
              card(
                card_header("Cause-Specific vs Fine-Gray"),
                plotOutput("plot_comp", height = "320px")
              ),
              
              card(
                card_header("Competing Risk Summary"),
                DT::dataTableOutput("tbl_comp")
              )
            )
  ),
  
  # ── Tab 7: Patient Explorer ───────────────────────────────────────────────
  nav_panel("Patient Explorer",
            DT::dataTableOutput("tbl_patients")
  )
)

server <- function(input, output, session) {
  
  # Overview
  output$plot_ecog <- renderPlot({
    cohort <- d$cohort %>%
      count(ecog_score)
    
    ggplot(cohort,
           aes(x = factor(ecog_score),
               y = n,
               fill = factor(ecog_score))) +
      geom_col() +
      theme_minimal(base_size = 13)
  })
  
  output$plot_lot <- renderPlot({
    cohort <- d$cohort %>%
      count(line_of_therapy)
    
    ggplot(cohort,
           aes(x = line_of_therapy,
               y = n,
               fill = line_of_therapy)) +
      geom_col() +
      theme_minimal(base_size = 13)
  })
  
  # HEOR Disparities
  output$plot_disp <- renderPlot({
    ggplot(d$therapy_disp,
           aes(x = race_ethnicity,
               y = therapy_rate,
               fill = insurance_type)) +
      geom_col(position = "dodge") +
      theme_minimal(base_size = 13)
  })
  
  output$plot_cost <- renderPlot({
    ggplot(d$cost_disp,
           aes(x = insurance_type,
               y = mean_cost,
               fill = insurance_type)) +
      geom_col() +
      theme_minimal(base_size = 13)
  })
  
  output$tbl_model <- DT::renderDataTable({
    DT::datatable(d$therapy_model,
                  options = list(pageLength = 10),
                  rownames = FALSE)
  })
  
  # Survival Analysis
  output$plot_km <- renderPlot({
    ggplot(d$km_interv,
           aes(x = time,
               y = estimate,
               color = strata)) +
      geom_step(size = 1) +
      theme_minimal(base_size = 13)
  })
  
  output$plot_cox <- renderPlot({
    cox <- d$cox_results %>%
      mutate(variable = forcats::fct_reorder(variable, hr))
    
    ggplot(cox,
           aes(x = variable,
               y = hr,
               ymin = hr_lower95,
               ymax = hr_upper95)) +
      geom_pointrange(color = PAL$treat) +
      geom_hline(yintercept = 1,
                 linetype = "dashed") +
      coord_flip() +
      theme_minimal(base_size = 13)
  })
  
  output$tbl_ph <- renderTable({
    d$ph_test
  })
  
  # Causal Inference
  output$plot_ps <- renderPlot({
    ps <- d$ps_scores %>%
      mutate(Group = ifelse(preventive_therapy == 1,
                            "Therapy","Control"))
    
    ggplot(ps,
           aes(x = ps,
               fill = Group)) +
      geom_histogram(alpha = 0.6,
                     bins = 40,
                     position = "identity") +
      theme_minimal(base_size = 13)
  })
  
  output$plot_balance <- renderPlot({
    bal <- d$balance %>%
      pivot_longer(cols = c(smd_unadjusted, smd_iptw),
                   names_to = "type",
                   values_to = "smd")
    
    ggplot(bal,
           aes(x = variable,
               y = smd,
               color = type,
               group = type)) +
      geom_point(size = 3) +
      geom_line() +
      theme_minimal(base_size = 13)
  })
  
  output$tbl_causal <- renderTable({
    d$causal
  })
  
  # Missing Data
  output$plot_mice <- renderPlot({
    ggplot(d$mice_compare,
           aes(x = analysis,
               y = hr,
               fill = analysis)) +
      geom_col() +
      theme_minimal(base_size = 13)
  })
  
  output$tbl_mice <- DT::renderDataTable({
    DT::datatable(d$mice_cox,
                  options = list(pageLength = 8),
                  rownames = FALSE)
  })
  
  # Competing Risks
  output$plot_comp <- renderPlot({
    ggplot(d$comp_risk,
           aes(x = model,
               y = estimate,
               fill = model)) +
      geom_col() +
      theme_minimal(base_size = 13)
  })
  
  output$tbl_comp <- DT::renderDataTable({
    DT::datatable(d$comp_risk,
                  options = list(pageLength = 5),
                  rownames = FALSE)
  })
  
  # Patient Explorer
  output$tbl_patients <- DT::renderDataTable({
    DT::datatable(
      d$cohort %>%
        select(patient_id,
               age,
               sex,
               ecog_score,
               biomarker_positive,
               line_of_therapy,
               preventive_therapy,
               disease_progression),
      options = list(pageLength = 15),
      rownames = FALSE
    )
  })
}

shinyApp(ui, server)