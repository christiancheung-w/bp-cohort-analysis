# 7_1_figure_3.R
# SAP §7.1 Figure 3: WHOOP-metric trajectories by responder status (3a SBP, 3b DBP).
#
# Standalone figure script: reads the frozen §7.1 result CSVs and regenerates the
# two trajectory figures — the LMMs are NOT refit here (mirrors the §5.1 / §5.5 /
# §7.2 figure-script separation). Reproduces the figures exactly as produced inside
# 7_1_whoop_trajectories.R, PLUS a per-panel omnibus-test annotation (bottom-left):
#
#     Group: p=…
#     Time: p=…
#     Interaction: p=…
#
# from the Type III ANOVA (7_1_anova.csv): responder_group main effect, Time main
# effect, and the Time × responder_group interaction. This is a vetting aid — the
# Time × responder_group interaction is the formal test of trajectory divergence,
# so these three p-values let the within/between/consecutive significance encodings
# be judged against whether an interaction is actually present.
#
# Inputs (Analyses/7_1/derived/results/):
#   7_1_emmeans.csv, 7_1_delta_contrasts.csv, 7_1_consec_contrasts.csv,
#   7_1_between_group_contrasts.csv, 7_1_anova.csv
# Outputs (Analyses/7_1/derived/figures/):
#   7_1_figure_3a_sbp.png, 7_1_figure_3b_dbp.png

library(tidyverse)
library(ggh4x)

project_dir <- "./PRIMARY"
results_dir <- file.path(project_dir, "Analyses", "7_1", "derived", "results")
figures_dir <- file.path(project_dir, "Analyses", "7_1", "derived", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

emmeans_tbl <- read_csv(file.path(results_dir, "7_1_emmeans.csv"),                show_col_types = FALSE)
delta_tbl   <- read_csv(file.path(results_dir, "7_1_delta_contrasts.csv"),         show_col_types = FALSE)
consec_tbl  <- read_csv(file.path(results_dir, "7_1_consec_contrasts.csv"),        show_col_types = FALSE)
between_tbl <- read_csv(file.path(results_dir, "7_1_between_group_contrasts.csv"), show_col_types = FALSE)
anova_tbl   <- read_csv(file.path(results_dir, "7_1_anova.csv"),                   show_col_types = FALSE)

# --- Figure dataframe with CIs + shared cosmetics (matches the generator) ---
fig_df_all <- emmeans_tbl %>%
  mutate(lcl = asymp.LCL, ucl = asymp.UCL) %>%
  mutate(
    Time = factor(Time, levels = c("baseline", "3mo", "6mo")),
    time_num = case_when(Time == "baseline" ~ 0, Time == "3mo" ~ 3, Time == "6mo" ~ 6),
    metric_label = factor(metric,
                          levels = c("RHR", "SleepConsistency", "SleepDuration", "Strain", "HRV"),
                          labels = c("RHR (bpm)", "Sleep Consistency (%)",
                                     "Sleep Duration (hours)", "Physical Activity (a.u.)", "HRV (ms)"))
  )

# Per-panel y scales — matched by metric_label factor-level order.
y_scales_sbp <- list(
  scale_y_continuous(limits = c(62.75, 64.25), breaks = c(63, 64)),                 # RHR
  scale_y_continuous(limits = c(69.0, 70.55), breaks = c(69.5, 70, 70.5)),          # Sleep Consistency (floor lowered to clear bottom-left annotation)
  scale_y_continuous(limits = c(6.96, 7.15),   breaks = c(7.00, 7.05, 7.10, 7.15)), # Sleep Duration (floor lowered slightly to clear bottom-left annotation)
  scale_y_continuous(limits = c(8, 9.5),       breaks = c(8, 8.5, 9, 9.5)),         # Physical Activity
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14)))                      # HRV
)
y_scales_wide <- list(
  scale_y_continuous(limits = c(62.5, 64.4), breaks = c(63, 64)),                   # RHR
  scale_y_continuous(limits = c(68.85, 70.6), breaks = c(69, 69.5, 70, 70.5)),      # Sleep Consistency (floor lowered to clear bottom-left annotation)
  scale_y_continuous(limits = c(6.95, 7.15),  breaks = c(7.00, 7.05, 7.10, 7.15)),  # Sleep Duration (floor lowered to clear bottom-left annotation)
  scale_y_continuous(limits = c(8, 9.6),     breaks = c(8, 8.5, 9, 9.5)),           # Physical Activity
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14)))                      # HRV
)

# Omnibus p-value formatter: truncate below 0.001, else 3 decimals.
fmt_omni <- function(p) ifelse(p < 0.001, "<0.001", sprintf("=%.3f", p))

make_parallel_figure <- function(df_def, responder_def, out_png, y_scales) {
  rdef <- responder_def
  resp_label <- sprintf("%s-Responder", responder_def)
  resp_color <- if (responder_def == "DBP") "#33a02c" else "#1f78b4"  # DBP green matches §6.2
  alpha_sig  <- 0.005
  time_of    <- c("baseline" = 0, "3mo" = 3, "6mo" = 6)

  # Interaction gate: the Time × responder_group interaction is §7.1's formal
  # test of trajectory divergence, so the post-hoc significance marks (point
  # fill, segment style, between-group asterisk) are only shown where it clears
  # α. Metrics whose interaction p ≥ α are rendered in the "fail-to-reject" style
  # — uniform OPEN points + DASHED segments, no asterisk — matching the figure's
  # own vocabulary (open = not rejected vs baseline; dashed = no rejected adjacent
  # change), so a non-rejected divergence null reads as "no claim." Inference is
  # carried by the omnibus p-block. (Avoids the difference-in-significance fallacy
  # where filled-vs-open within-group marks imply a divergence the interaction
  # does not license.)
  gated_metrics <- anova_tbl %>%
    filter(responder_def == rdef, term == "Time:responder_group",
           `Pr(>F)` >= alpha_sig) %>%
    pull(metric)

  # Point fill: significant vs BASELINE (delta / trt.vs.ctrl); baseline = empty (reference).
  # Uses Holm-adjusted p (p.holm) per SAP v0.84 — consistent with 7_1_whoop_trajectories.R.
  pt_sig <- delta_tbl %>% filter(responder_def == rdef) %>%
    transmute(metric, responder_group,
              time_num = unname(time_of[sub(" .*$", "", contrast)]),
              sig_base = p.holm < alpha_sig)
  df_plot <- df_def %>%
    mutate(responder_label = factor(responder_group,
                                    levels = c("non_responder", "responder"),
                                    labels = c("Non-Responder", resp_label))) %>%
    left_join(pt_sig, by = c("metric", "responder_group", "time_num")) %>%
    mutate(pt_filled = dplyr::coalesce(sig_base, FALSE),
           pt_filled = ifelse(metric %in% gated_metrics, FALSE, pt_filled))  # gated → uniform open (fail to reject)
  group_cols <- setNames(c("#999999", resp_color), c("Non-Responder", resp_label))

  # Line segments: solid if significant change vs previous timepoint (consec), else dashed.
  em_lk <- df_plot %>% distinct(metric, metric_label, responder_group, responder_label, time_num, emmean)
  seg <- consec_tbl %>% filter(responder_def == rdef) %>%
    tidyr::separate(contrast, into = c("to", "from"), sep = " - ", remove = TRUE) %>%
    transmute(metric, responder_group,
              x = unname(time_of[from]), xend = unname(time_of[to]),
              seg_sig = p.holm < alpha_sig) %>%
    left_join(em_lk, by = c("metric", "responder_group", "x" = "time_num")) %>% rename(y = emmean) %>%
    left_join(em_lk %>% select(metric, responder_group, time_num, emmean),
              by = c("metric", "responder_group", "xend" = "time_num")) %>% rename(yend = emmean) %>%
    mutate(seg_sig = ifelse(metric %in% gated_metrics, FALSE, seg_sig))  # gated → uniform dashed (fail to reject)

  # Asterisk: significant Responder vs Non-Responder difference at a timepoint.
  ast <- between_tbl %>%
    filter(responder_def == rdef, p.holm < alpha_sig, !(metric %in% gated_metrics)) %>%
    transmute(metric, time_num = unname(time_of[Time])) %>%
    left_join(df_def %>% group_by(metric, metric_label, time_num) %>%
                summarise(y = max(ucl), .groups = "drop"),
              by = c("metric", "time_num"))

  # Omnibus annotation (bottom-left of each panel): Group / Time / Interaction p-values.
  omni_df <- anova_tbl %>%
    filter(responder_def == rdef,
           term %in% c("responder_group", "Time", "Time:responder_group")) %>%
    select(metric, term, p = `Pr(>F)`) %>%
    pivot_wider(names_from = term, values_from = p) %>%
    mutate(
      metric_label = factor(metric,
                            levels = c("RHR", "SleepConsistency", "SleepDuration", "Strain", "HRV"),
                            labels = c("RHR (bpm)", "Sleep Consistency (%)",
                                       "Sleep Duration (hours)", "Physical Activity (a.u.)", "HRV (ms)")),
      label = sprintf("Group: p%s\nTime: p%s\nInteraction: p%s",
                      fmt_omni(responder_group), fmt_omni(Time), fmt_omni(`Time:responder_group`))
    )

  p <- ggplot(df_plot, aes(x = time_num, y = emmean,
                           color = responder_label, fill = responder_label,
                           group = responder_label)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.18, color = NA) +
    geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend,
                                 color = responder_label, linetype = seg_sig),
                 linewidth = 0.9, inherit.aes = FALSE) +
    geom_point(aes(shape = pt_filled), size = 2.5) +
    geom_text(data = ast, aes(x = time_num, y = y, label = "*"),
              vjust = 0, size = 6, fontface = "bold", color = "black", inherit.aes = FALSE) +
    geom_text(data = omni_df, aes(label = label), x = -Inf, y = -Inf,
              hjust = -0.05, vjust = -0.35, size = 2.6, lineheight = 0.95,
              color = "gray20", inherit.aes = FALSE) +
    facet_manual(vars(metric_label),
                 design = "AE#\nAED\nCBD\nCB#",
                 scales = "free_y", strip.position = "left", axes = "x") +
    facetted_pos_scales(y = y_scales) +
    scale_x_continuous(breaks = c(0, 3, 6),
                       labels = c("Baseline", "3 Months", "6 Months")) +
    scale_color_manual(values = group_cols, name = NULL,
                       breaks = c(resp_label, "Non-Responder")) +
    scale_fill_manual(values  = group_cols, name = NULL,
                      breaks = c(resp_label, "Non-Responder")) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed"), guide = "none") +
    labs(
      x = NULL, y = NULL,
      title    = "Six-Month Wearable Derived Biometric Trajectories by Responder Status",
      subtitle = "Adjusted LMM estimated marginal means with 95% CI"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid        = element_blank(),
      panel.spacing     = grid::unit(1.5, "lines"),
      plot.margin       = margin(t = 6, r = 16, b = 6, l = 6),
      axis.line         = element_line(color = "black", linewidth = 0.4),
      axis.ticks        = element_line(color = "black", linewidth = 0.4),
      axis.text.x       = element_text(face = "bold"),
      strip.placement   = "outside",
      strip.background  = element_blank(),
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 11),
      legend.position   = "bottom",
      plot.title        = element_text(face = "bold")
    )
  ggsave(out_png, p, width = 11, height = 7.5, dpi = 300)
  cat("Wrote figure:\n  ", out_png, "\n")
}

make_parallel_figure(
  fig_df_all %>% filter(responder_def == "SBP"),
  "SBP",
  file.path(figures_dir, "7_1_figure_3a_sbp.png"),
  y_scales = y_scales_sbp
)
make_parallel_figure(
  fig_df_all %>% filter(responder_def == "DBP"),
  "DBP",
  file.path(figures_dir, "7_1_figure_3b_dbp.png"),
  y_scales = y_scales_wide
)

# Echo the omnibus table for vetting.
cat("\nOmnibus p-values (Type III ANOVA) by responder_def × metric:\n")
print(anova_tbl %>%
        filter(term %in% c("responder_group", "Time", "Time:responder_group")) %>%
        select(responder_def, metric, term, `F value`, `Pr(>F)`) %>%
        arrange(responder_def, metric, term), n = Inf, width = Inf)
