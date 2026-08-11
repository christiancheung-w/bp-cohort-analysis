# 8_figure_5.R
# SAP §8 Figure 5: engagement-metric distributions by BP-responder status
# (5a SBP-responder, 5b DBP-responder).
#
# Standalone figure script (house-style pass): reads the frozen user analytic
# frame (read-only) and writes only the two PNGs — the §8 generator
# (8_engagement_logistic.R) and all §8 result CSVs stay frozen. Mirrors the
# §5.1 / §5.5 / §7.1 / §7.2 figure-script separation.
#
# House style: publication title (no §-prefix), short descriptive subtitle, no
# caption; gridlines off + black axes/ticks; cleaned Title-Case panel labels;
# centered 4-over-3 panel layout via ggh4x::facet_manual; responder colours
# matching §6.2 / §7.1 (SBP #1f78b4, DBP #33a02c; non-responder #999999);
# legend bottom, responder first, capitalised. (n per panel reported in text,
# not on the figure — all panels now ≈ 19–20k users; AVG_SLEEP_NEED_PCT replaced
# the sparse SLEEP_PCT_OF_TARGET, so the sleep panel is no longer the ≈7.5k outlier.)
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
# Output: Analyses/8_engagement/derived/figures/8_figure_5_violin_{sbp,dbp}.png

library(tidyverse)
library(ggh4x)

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
figures_dir <- file.path(project_dir, "Analyses", "8_engagement", "derived", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# --- Cohort (matches the §8 generator: §2.2 + sex restriction; DBP responder derived) ---
cohort <- read_rds(input_path) %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female")) %>%
  mutate(DBP_RESPONDER_6MO = as.integer((MO6_DBP - BASELINE_DBP) <= -5))

# --- Seven engagement metrics (raw 0–100 / log scales) + cleaned panel labels ---
engagement_preds_violin <- c(
  "PCT_DAYS_JOURNALED",
  "LOG_AVG_DAILY_TIME_IN_APP_MIN",
  "PCT_DAYS_ENGAGED_NOTIFICATIONS",
  "PCT_DAYS_USED_COACH",
  "AVG_SLEEP_NEED_PCT",
  "PCT_DAYS_MET_OPTIMAL_STRAIN",
  "LOG_N_BP_READINGS_STUDY"
)
pred_labels <- c(
  PCT_DAYS_JOURNALED             = "Days Journaled (%)",
  LOG_AVG_DAILY_TIME_IN_APP_MIN  = "Daily App Time (log min)",
  PCT_DAYS_ENGAGED_NOTIFICATIONS = "Days w/ Notifications (%)",
  PCT_DAYS_USED_COACH            = "Days Used AI Coach (%)",
  AVG_SLEEP_NEED_PCT             = "Average Sleep Need Achieved (%)",
  PCT_DAYS_MET_OPTIMAL_STRAIN    = "Days Met Optimal Strain (%)",
  LOG_N_BP_READINGS_STUDY        = "BP Readings (log count)"
)

# Centered 4-over-3 layout: top row 4 panels, bottom row 3 centered (each panel
# spans 2 cols; bottom row padded 1 col each side). Letters A–G map to the 7
# panel labels in factor order.
# Metrics definitionally bounded at 0–100% — y-axis fixed at 0–100.
pct_bounded <- c("PCT_DAYS_JOURNALED", "PCT_DAYS_ENGAGED_NOTIFICATIONS",
                 "PCT_DAYS_USED_COACH", "PCT_DAYS_MET_OPTIMAL_STRAIN")

panel_design <- "AABBCCDD\n#EEFFGG#"

make_violin <- function(outcome_tag, out_png) {
  resp_col   <- paste0(outcome_tag, "_RESPONDER_6MO")
  resp_label <- paste0(outcome_tag, "-Responder")
  fig_color  <- if (outcome_tag == "DBP") "#33a02c" else "#1f78b4"

  vd <- cohort %>%
    filter(!is.na(.data[[resp_col]])) %>%
    select(all_of(c(resp_col, engagement_preds_violin))) %>%
    mutate(responder_label = factor(
      ifelse(.data[[resp_col]] == 1L, resp_label, "Non-Responder"),
      levels = c("Non-Responder", resp_label))) %>%
    pivot_longer(all_of(engagement_preds_violin),
                 names_to = "predictor", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(predictor_label = factor(pred_labels[predictor], levels = pred_labels))

  # Fixed axis bounds per panel; no data clipping — full tails included.
  # pct_bounded: 0–100; AVG_SLEEP_NEED_PCT: 50–100 (realistic bound); log metrics: free.
  limits_df <- bind_rows(
    tibble(predictor = rep(pct_bounded, each = 2),
           value     = rep(c(0, 100), length(pct_bounded))),
    tibble(predictor = "AVG_SLEEP_NEED_PCT", value = c(50, 100))
  ) %>%
    mutate(predictor_label = factor(pred_labels[predictor], levels = pred_labels),
           responder_label = levels(vd$responder_label)[[1]])

  group_cols <- setNames(c("#999999", fig_color), c("Non-Responder", resp_label))

  p <- ggplot(vd, aes(x = responder_label, y = value, fill = responder_label)) +
    geom_blank(data = limits_df, aes(x = responder_label, y = value),
               inherit.aes = FALSE) +
    geom_violin(alpha = 0.75, scale = "width", trim = TRUE,
                color = "gray30", linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
                 color = "gray30", linewidth = 0.3) +
    facet_manual(vars(predictor_label), design = panel_design,
                 scales = "free_y", axes = "all") +
    scale_fill_manual(values = group_cols, name = NULL,
                      breaks = c(resp_label, "Non-Responder")) +
    labs(
      x = NULL, y = NULL,
      title    = "Six-Month WHOOP Engagement Metrics by Responder Status",
      subtitle = "Per-user study-window distributions; inner box = median and IQR"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid       = element_blank(),
      panel.spacing    = grid::unit(1.2, "lines"),
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      axis.text.x      = element_text(size = 8.5),
      strip.text       = element_text(face = "bold", size = 10),
      legend.position  = "bottom",
      plot.title       = element_text(face = "bold")
    )
  ggsave(out_png, p, width = 12, height = 7.5, dpi = 300)
  cat("Wrote figure:\n  ", out_png, "\n")
}

make_violin("SBP", file.path(figures_dir, "8_figure_5_violin_sbp.png"))
make_violin("DBP", file.path(figures_dir, "8_figure_5_violin_dbp.png"))
