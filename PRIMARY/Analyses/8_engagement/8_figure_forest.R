# 8_figure_forest.R
# SAP §8 forest plots: WHOOP-engagement predictors of BP-responder status.
#   primary   = 7 separate single-predictor adjusted logistic models per outcome.
#   secondary = 1 joint mutually-adjusted model per outcome (~7.5k complete-case).
#
# House-style pass mirroring §6.1 / §6.2: publication title, short methods
# subtitle, no caption; boxed panel, no y-gridlines, bold x-title; SBP blue /
# DBP green dodged; filled marker = significant at α = 0.005. Predictor labels
# carry no units; the OR unit (+10 pp for the % metrics, per-log-unit for the two
# log metrics) is documented in Methods / the table.
#
# Standalone: reads the frozen §8 coefficient CSVs, writes only the two PNGs.
#
# Inputs:  results/8_coefficients_separate_primary.csv, results/8_coefficients_joint_primary.csv
# Outputs: figures/8_forest_separate_primary.png, figures/8_forest_joint_secondary.png

library(tidyverse)

project_dir <- "./PRIMARY"
results_dir <- file.path(project_dir, "Analyses", "8_engagement", "derived", "results")
figures_dir <- file.path(project_dir, "Analyses", "8_engagement", "derived", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Clean labels keyed by the model term (per-10pp / log predictor names). OR units
# (+10 pp for the % metrics, per-log-unit for the two log metrics) are documented
# in Methods / the table, not on the figure.
forest_labels <- c(
  LOG_AVG_DAILY_TIME_IN_APP_MIN        = "Daily App Time",
  PCT_DAYS_JOURNALED_per10             = "Journal Completion",
  PCT_DAYS_ENGAGED_NOTIFICATIONS_per10 = "Notification Engagement",
  PCT_DAYS_USED_COACH_per10            = "AI Coach Engagement",
  AVG_SLEEP_NEED_PCT_per10             = "Sleep Need Achieved",
  PCT_DAYS_MET_OPTIMAL_STRAIN_per10    = "Physical Activity Recommendation",
  LOG_N_BP_READINGS_STUDY              = "Blood Pressure Readings")
order_top_to_bottom <- unname(forest_labels)   # y-axis top→bottom = this order

outcome_cols <- c("SBP" = "#1f78b4", "DBP" = "#33a02c")

build_forest <- function(df, title_txt, subtitle_txt) {
  df <- df %>%
    filter(term %in% names(forest_labels)) %>%
    mutate(
      term_label = factor(forest_labels[term], levels = rev(order_top_to_bottom)),
      outcome    = factor(outcome, levels = c("DBP", "SBP")))   # SBP listed last → drawn on top
  ggplot(df, aes(x = OR, y = term_label, color = outcome, shape = sig_005)) +
    geom_vline(xintercept = 1, color = "gray60", linetype = "dashed", linewidth = 0.4) +
    geom_errorbar(aes(xmin = OR_lcl, xmax = OR_ucl), width = 0.4, linewidth = 0.5,
                  orientation = "y", position = position_dodge(width = 0.6)) +
    geom_point(size = 3, position = position_dodge(width = 0.6)) +
    scale_x_log10(breaks = seq(0.9, 1.4, 0.1), expand = expansion(mult = 0.03)) +
    scale_color_manual(values = outcome_cols, name = NULL,
                       breaks = c("SBP", "DBP"),
                       labels = c(SBP = "SBP-Responder", DBP = "DBP-Responder")) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    labs(x = "Adjusted OR (Log Scale, 95% CI)", y = NULL,
         title = title_txt, subtitle = subtitle_txt) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position    = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.ticks         = element_line(color = "black", linewidth = 0.4),
      axis.title.x       = element_text(face = "bold", margin = margin(t = 10)),
      plot.title         = element_text(face = "bold", size = 11))
}

sep   <- read_csv(file.path(results_dir, "8_coefficients_separate_primary.csv"), show_col_types = FALSE)
joint <- read_csv(file.path(results_dir, "8_coefficients_joint_primary.csv"),    show_col_types = FALSE)

p_sep <- build_forest(sep,
  "WHOOP Engagement Predictors of Clinically Meaningful Blood Pressure Reductions",
  "Per-Predictor Adjusted Logistic Regression")
p_joint <- build_forest(joint,
  "WHOOP Engagement Predictors of Clinically Meaningful Blood Pressure Reductions",
  "Mutually-Adjusted Multivariable Logistic Regression")

ggsave(file.path(figures_dir, "8_forest_separate_primary.png"), p_sep,   width = 9.5, height = 6.5, dpi = 300)
ggsave(file.path(figures_dir, "8_forest_joint_secondary.png"),  p_joint, width = 9.5, height = 6.5, dpi = 300)
cat("Wrote:\n  8_forest_separate_primary.png\n  8_forest_joint_secondary.png\n")
