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
    cohort         = read_csv("../data/cohort.csv", show_col_types = FALSE,
                              col_types = cols(line_of_therapy = col_character())),
    claims         = read_csv("../data/claims.csv",                              show_col_types = FALSE),
    therapy_disp   = read_csv("../outputs/therapy_disparities.csv",              show_col_types = FALSE),
    therapy_model  = read_csv("../outputs/therapy_access_model.csv",             show_col_types = FALSE),
    cost_disp      = read_csv("../outputs/cost_disparities.csv",                 show_col_types = FALSE),
    km_interv      = read_csv("../outputs/km_intervention.csv",                  show_col_types = FALSE),
    cox_results    = read_csv("../outputs/cox_results.csv",                      show_col_types = FALSE),
    ph_test        = read_csv("../outputs/ph_test.csv",                          show_col_types = FALSE),
    balance        = read_csv("../outputs/covariate_balance.csv",                show_col_types = FALSE),
    ps_scores      = read_csv("../outputs/propensity_scores.csv",                show_col_types = FALSE),
    causal         = read_csv("../outputs/causal_estimates.csv",                 show_col_types = FALSE),
    causal_kpi     = read_csv("../outputs/causal_key_metrics.csv",               show_col_types = FALSE),
    mice_compare   = read_csv("../outputs/complete_case_vs_mice.csv",            show_col_types = FALSE),
    mice_cox       = read_csv("../outputs/mice_pooled_cox.csv",                  show_col_types = FALSE),
    comp_risk      = read_csv("../outputs/competing_risk_model_comparison.csv",  show_col_types = FALSE),
    cif_therapy    = read_csv("../outputs/cif_by_therapy.csv",                   show_col_types = FALSE),
    fg_results     = read_csv("../outputs/fine_gray_results.csv",                show_col_types = FALSE),
    grays_test     = read_csv("../outputs/grays_test.csv",                       show_col_types = FALSE),
    subgroup_hr    = read_csv("../outputs/subgroup_hr.csv",                      show_col_types = FALSE),
    interaction    = read_csv("../outputs/interaction_test.csv",                 show_col_types = FALSE),
    km_biomarker   = read_csv("../outputs/km_biomarker_subgroup.csv",            show_col_types = FALSE),
    table1         = read_csv("../outputs/table1_baseline.csv",                  show_col_types = FALSE)
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
  
  # ── Tab 1: Overview ──────────────────────────────────────────────────────
  nav_panel("Overview",
            layout_columns(
              col_widths = c(3, 3, 3, 3),
              value_box("Cohort Size",
                        scales::comma(nrow(d$cohort)),
                        showcase = bsicons::bs_icon("people-fill"), theme = "primary"),
              value_box("Progression Rate",
                        paste0(round(100 * mean(d$cohort$disease_progression), 1), "%"),
                        showcase = bsicons::bs_icon("activity"), theme = "danger"),
              value_box("Therapy Exposure",
                        paste0(round(100 * mean(d$cohort$preventive_therapy), 1), "%"),
                        showcase = bsicons::bs_icon("capsule"), theme = "success"),
              value_box("Median Follow-Up",
                        paste0(round(median(d$cohort$observed_time), 1), " months"),
                        showcase = bsicons::bs_icon("calendar2-week"), theme = "warning")
            ),
            layout_columns(
              col_widths = c(6, 6),
              card(
                card_header("Progression Rate by ECOG Score"),
                plotOutput("plot_ecog_prog", height = "300px")
              ),
              card(
                card_header("Therapy Exposure by Line of Therapy"),
                plotOutput("plot_lot_therapy", height = "300px")
              )
            )
  ),
  
  # ── Tab 2: Baseline Characteristics ─────────────────────────────────────
  nav_panel("Baseline",
            card(
              card_header("Table 1: Baseline Characteristics by Therapy Group"),
              DT::dataTableOutput("tbl_table1")
            )
  ),
  
  # ── Tab 3: HEOR Disparities ──────────────────────────────────────────────
  nav_panel("HEOR Disparities",
            layout_columns(
              col_widths = c(6, 6),
              card(
                card_header("Therapy Access by Race/Ethnicity and Insurance"),
                plotOutput("plot_therapy_disp", height = "360px")
              ),
              card(
                card_header("Mean Cost by Insurance Type"),
                plotOutput("plot_cost_disp", height = "360px")
              )
            )
  ),
  
  # ── Tab 4: Survival Analysis ─────────────────────────────────────────────
  nav_panel("Survival Analysis",
            layout_columns(
              col_widths = c(7, 5),
              card(
                card_header("Kaplan-Meier: Progression-Free Survival by Therapy"),
                plotOutput("plot_km", height = "360px")
              ),
              card(
                card_header("Cox PH Model: Hazard Ratios (95% CI)"),
                plotOutput("plot_cox_forest", height = "360px")
              )
            ),
            card(
              card_header("Proportional Hazards Assumption (Schoenfeld)"),
              DT::dataTableOutput("tbl_ph")
            )
  ),
  
  # ── Tab 5: Causal Inference ──────────────────────────────────────────────
  nav_panel("Causal Inference",
            layout_columns(
              col_widths = c(6, 6),
              card(
                card_header("Propensity Score Overlap"),
                plotOutput("plot_ps", height = "300px")
              ),
              card(
                card_header("Covariate Balance (SMD)"),
                plotOutput("plot_balance", height = "300px")
              )
            ),
            layout_columns(
              col_widths = c(6, 6),
              card(
                card_header("IPTW vs AIPW: ATE with 95% CI"),
                plotOutput("plot_causal_compare", height = "280px")
              ),
              card(
                card_header("Key Causal Metrics"),
                DT::dataTableOutput("tbl_causal_kpi")
              )
            )
  ),
  
  # ── Tab 6: Missing Data ──────────────────────────────────────────────────
  nav_panel("Missing Data",
            layout_columns(
              col_widths = c(6, 6),
              card(
                card_header("Complete Case vs MICE: Preventive Therapy HR"),
                plotOutput("plot_mice", height = "300px")
              ),
              card(
                card_header("MICE Pooled Cox Results"),
                DT::dataTableOutput("tbl_mice")
              )
            )
  ),
  
  # ── Tab 7: Competing Risks ───────────────────────────────────────────────
  nav_panel("Competing Risks",
            layout_columns(
              col_widths = c(6, 6),
              card(
                card_header("Cumulative Incidence: Progression by Therapy Group"),
                plotOutput("plot_cif", height = "320px")
              ),
              card(
                card_header("Cause-Specific vs Fine-Gray HR"),
                plotOutput("plot_comp", height = "320px")
              )
            ),
            layout_columns(
              col_widths = c(6, 6),
              card(
                card_header("Model Comparison Table"),
                DT::dataTableOutput("tbl_comp")
              ),
              card(
                card_header("Gray's Test"),
                DT::dataTableOutput("tbl_grays")
              )
            )
  ),
  
  # ── Tab 8: Subgroup Analysis ─────────────────────────────────────────────
  nav_panel("Subgroup Analysis",
            layout_columns(
              col_widths = c(7, 5),
              card(
                card_header("KM Curves by Biomarker Status and Therapy"),
                plotOutput("plot_km_biomarker", height = "360px")
              ),
              card(
                card_header("Subgroup Forest Plot"),
                plotOutput("plot_subgroup_forest", height = "360px")
              )
            ),
            layout_columns(
              col_widths = c(6, 6),
              card(
                card_header("Subgroup Hazard Ratios"),
                DT::dataTableOutput("tbl_subgroup")
              ),
              card(
                card_header("Interaction Test: Therapy × Biomarker"),
                DT::dataTableOutput("tbl_interaction")
              )
            )
  )
)

server <- function(input, output, session) {
  
  # ── Overview ────────────────────────────────────────────────────────────
  output$plot_ecog_prog <- renderPlot({
    d$cohort %>%
      group_by(ecog_score) %>%
      summarise(prog_rate = round(100 * mean(disease_progression), 1),
                .groups = "drop") %>%
      ggplot(aes(x = factor(ecog_score), y = prog_rate, fill = factor(ecog_score))) +
      geom_col() +
      scale_fill_brewer(palette = "Reds") +
      labs(x = "ECOG Score", y = "Progression Rate (%)", fill = NULL) +
      theme_minimal(base_size = 13) + theme(legend.position = "none")
  })
  
  output$plot_lot_therapy <- renderPlot({
    lot_order <- c("1L", "2L", "3L+")
    lot_colors <- c("1L" = PAL$green, "2L" = PAL$gold, "3L+" = PAL$control)
    df <- d$cohort %>%
      filter(!is.na(line_of_therapy)) %>%
      mutate(lot = factor(as.character(line_of_therapy), levels = lot_order)) %>%
      group_by(lot, .drop = FALSE) %>%
      summarise(
        therapy_rate = round(100 * mean(preventive_therapy, na.rm = TRUE), 1),
        .groups = "drop"
      ) %>%
      filter(!is.na(lot))
    validate(need(nrow(df) > 0, "line_of_therapy not available"))
    ggplot(df, aes(x = lot, y = therapy_rate, fill = lot)) +
      geom_col() +
      scale_fill_manual(values = lot_colors, drop = FALSE) +
      labs(x = "Line of Therapy", y = "Therapy Exposure (%)", fill = NULL) +
      theme_minimal(base_size = 13) + theme(legend.position = "none")
  })
  
  # ── Baseline ────────────────────────────────────────────────────────────
  output$tbl_table1 <- DT::renderDataTable({
    DT::datatable(d$table1, options = list(pageLength = 5), rownames = FALSE)
  })
  
  # ── HEOR ────────────────────────────────────────────────────────────────
  output$plot_therapy_disp <- renderPlot({
    d$therapy_disp %>%
      ggplot(aes(x = race_ethnicity, y = therapy_rate, fill = insurance_type)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c(Commercial = PAL$green,
                                   Medicaid   = PAL$gold,
                                   Medicare   = PAL$treat)) +
      labs(x = NULL, y = "Therapy Access Rate (%)", fill = "Insurance") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  })
  
  output$plot_cost_disp <- renderPlot({
    d$cost_disp %>%
      ggplot(aes(x = insurance_type, y = mean_cost, fill = insurance_type)) +
      geom_col() +
      scale_fill_manual(values = c(Commercial = PAL$green,
                                   Medicaid   = PAL$gold,
                                   Medicare   = PAL$treat)) +
      labs(x = NULL, y = "Mean Cost (USD)", fill = NULL) +
      theme_minimal(base_size = 13) + theme(legend.position = "none")
  })
  
  # ── Survival ────────────────────────────────────────────────────────────
  output$plot_km <- renderPlot({
    d$km_interv %>%
      mutate(Group = ifelse(grepl("=1", strata), "Therapy", "Control")) %>%
      ggplot(aes(x = time, y = estimate, color = Group)) +
      geom_step(linewidth = 1.1) +
      scale_color_manual(values = c(Control = PAL$control, Therapy = PAL$treat)) +
      labs(x = "Months", y = "Progression-Free Survival", color = NULL) +
      theme_minimal(base_size = 13)
  })
  
  output$plot_cox_forest <- renderPlot({
    d$cox_results %>%
      mutate(variable  = forcats::fct_reorder(variable, hr),
             sig       = ifelse(p_value < 0.05, "p < 0.05", "p ≥ 0.05")) %>%
      ggplot(aes(x = variable, y = hr,
                 ymin = hr_lower95, ymax = hr_upper95,
                 color = sig)) +
      geom_pointrange(size = 0.7) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      coord_flip() +
      scale_color_manual(values = c("p < 0.05" = PAL$dark, "p ≥ 0.05" = PAL$grey)) +
      labs(x = NULL, y = "Hazard Ratio (95% CI)", color = NULL) +
      theme_minimal(base_size = 13)
  })
  
  output$tbl_ph <- DT::renderDataTable({
    DT::datatable(d$ph_test, options = list(pageLength = 8), rownames = FALSE)
  })
  
  # ── Causal Inference ────────────────────────────────────────────────────
  output$plot_ps <- renderPlot({
    d$ps_scores %>%
      mutate(Group = ifelse(preventive_therapy == 1, "Therapy", "Control")) %>%
      ggplot(aes(x = ps, fill = Group)) +
      geom_histogram(alpha = 0.6, bins = 40, position = "identity") +
      scale_fill_manual(values = c(Therapy = PAL$treat, Control = PAL$control)) +
      labs(x = "Propensity Score", y = "Count", fill = NULL) +
      theme_minimal(base_size = 13)
  })
  
  output$plot_balance <- renderPlot({
    d$balance %>%
      pivot_longer(cols = c(smd_unadjusted, smd_iptw),
                   names_to = "type", values_to = "smd") %>%
      mutate(
        type = recode(type,
                      smd_unadjusted = "Unadjusted",
                      smd_iptw       = "IPTW Weighted"),
        smd  = ifelse(is.na(smd), 0, smd)
      ) %>%
      mutate(variable = forcats::fct_reorder(variable, smd)) %>%
      ggplot(aes(x = smd, y = variable, color = type, shape = type)) +
      geom_point(size = 3) +
      geom_vline(xintercept = 0.10, linetype = "dashed", color = "red") +
      labs(x = "Absolute SMD", y = NULL, color = NULL, shape = NULL) +
      theme_minimal(base_size = 13)
  })
  
  output$plot_causal_compare <- renderPlot({
    causal_df <- d$causal %>%
      mutate(estimator = fct_inorder(estimator))
    # Handle both old CSV (no CIs) and new CSV (with CIs)
    has_ci <- all(c("ci_lower", "ci_upper") %in% names(causal_df))
    if (has_ci) {
      p <- causal_df %>%
        ggplot(aes(x = estimator, y = ate,
                   ymin = ci_lower, ymax = ci_upper,
                   color = estimator)) +
        geom_pointrange(size = 0.9)
    } else {
      p <- causal_df %>%
        ggplot(aes(x = estimator, y = ate, color = estimator)) +
        geom_point(size = 3)
    }
    p +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      coord_flip() +
      scale_color_manual(values = c(PAL$treat, PAL$green)) +
      labs(x = NULL, y = "Average Treatment Effect (ATE)") +
      theme_minimal(base_size = 13) + theme(legend.position = "none")
  })
  
  output$tbl_causal_kpi <- DT::renderDataTable({
    DT::datatable(d$causal_kpi, options = list(pageLength = 8), rownames = FALSE)
  })
  
  # ── Missing Data ────────────────────────────────────────────────────────
  output$plot_mice <- renderPlot({
    d$mice_compare %>%
      ggplot(aes(x = analysis, y = hr, ymin = lower95, ymax = upper95)) +
      geom_pointrange(color = PAL$treat, size = 0.9) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      coord_flip() +
      labs(title = "Sensitivity Analysis: Complete Case vs MICE",
           x = NULL, y = "Hazard Ratio (95% CI)") +
      theme_minimal(base_size = 13)
  })
  
  output$tbl_mice <- DT::renderDataTable({
    DT::datatable(d$mice_cox, options = list(pageLength = 8), rownames = FALSE)
  })
  
  # ── Competing Risks ─────────────────────────────────────────────────────
  output$plot_cif <- renderPlot({
    d$cif_therapy %>%
      filter(event == "Disease Progression") %>%
      ggplot(aes(x = time, y = cif, color = group)) +
      geom_step(linewidth = 1.1) +
      scale_color_manual(values = c(Control = PAL$control, Therapy = PAL$treat)) +
      labs(x = "Months", y = "Cumulative Incidence of Progression", color = NULL) +
      theme_minimal(base_size = 13)
  })
  
  output$plot_comp <- renderPlot({
    d$comp_risk %>%
      mutate(model = fct_inorder(model)) %>%
      ggplot(aes(x = model, y = estimate, ymin = lower95, ymax = upper95,
                 color = model)) +
      geom_pointrange(size = 0.9) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      coord_flip() +
      scale_color_manual(values = c(PAL$purple, PAL$treat)) +
      labs(x = NULL, y = "Hazard / Subdistribution Hazard Ratio") +
      theme_minimal(base_size = 13) + theme(legend.position = "none")
  })
  
  output$tbl_comp <- DT::renderDataTable({
    DT::datatable(
      d$comp_risk %>% select(model, estimate, lower95, upper95),
      options = list(pageLength = 5), rownames = FALSE
    )
  })
  
  output$tbl_grays <- DT::renderDataTable({
    DT::datatable(d$grays_test, options = list(pageLength = 5), rownames = FALSE)
  })
  
  # ── Subgroup Analysis ───────────────────────────────────────────────────
  output$plot_km_biomarker <- renderPlot({
    d$km_biomarker %>%
      mutate(
        therapy   = ifelse(grepl("preventive_therapy=1", strata), "Therapy", "Control"),
        biomarker = ifelse(grepl("biomarker_positive=1", strata), "Biomarker+", "Biomarker-")
      ) %>%
      ggplot(aes(x = time, y = estimate, color = therapy)) +
      geom_step(linewidth = 1) +
      facet_wrap(~ biomarker) +
      scale_color_manual(values = c(Control = PAL$control, Therapy = PAL$treat)) +
      labs(x = "Months", y = "Progression-Free Survival", color = NULL) +
      theme_minimal(base_size = 13)
  })
  
  output$plot_subgroup_forest <- renderPlot({
    int_p <- d$interaction$interaction_p[1]
    d$subgroup_hr %>%
      mutate(subgroup = fct_rev(fct_inorder(subgroup))) %>%
      ggplot(aes(x = subgroup, y = hr, ymin = lower95, ymax = upper95,
                 color = subgroup)) +
      geom_pointrange(size = 0.9) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      coord_flip() +
      scale_color_manual(values = c(
        "Biomarker Positive" = PAL$green,
        "Biomarker Negative" = PAL$control,
        "Overall"            = PAL$treat
      )) +
      annotate("text",
               x    = 0.6,
               y    = max(d$subgroup_hr$upper95, na.rm = TRUE),
               label = sprintf("Interaction p = %.4f", int_p),
               hjust = 1, size = 3.5, color = "grey30") +
      labs(x = NULL, y = "Hazard Ratio (95% CI)") +
      theme_minimal(base_size = 13) + theme(legend.position = "none")
  })
  
  output$tbl_subgroup <- DT::renderDataTable({
    DT::datatable(d$subgroup_hr, options = list(pageLength = 5), rownames = FALSE)
  })
  
  output$tbl_interaction <- DT::renderDataTable({
    DT::datatable(d$interaction, options = list(pageLength = 5), rownames = FALSE)
  })
}

shinyApp(ui, server)