# 5_5_figure_2.R
# SAP §5.5 Figure 2: BP trajectory by AHA baseline stratum, 4-timepoint
# adjusted LMM (baseline / 3mo / 6mo / 12mo) — the reported lead, matching Fig 1A
# / Table 2 (adjusted model incl. baseline-HTN per CHANGELOG v0.79).
#
# Mirrors the Figure 1A house style (5_1_figure_1a.R) extended to the MO12
# timepoint: AHA-chart stratum palette (navy-grey Elevated), rotated full-name
# left y-axis titles, no gridlines + visible black axes/ticks, 5-mmHg breaks,
# dotted AHA cutoffs (unlabelled). Wider CI ribbons at
# MO12 (smaller analytic subset) are part of the visual story; the MO12 n is
# reported in the manuscript text, not on the figure.
#
# Input:  Analyses/5_5/derived/emmeans_results.rds
# Output: Analyses/5_5/derived/figures/5_5_figure_2.png

library(tidyverse)

project_dir <- "./PRIMARY"
in_path     <- file.path(project_dir, "Analyses", "5_5", "derived", "emmeans_results.rds")
fig_dir     <- file.path(project_dir, "Analyses", "5_5", "derived", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

emmeans_results <- read_rds(in_path)

# CI column detection (handles asymp.LCL/UCL or lower.CL/upper.CL)
get_lcl_ucl <- function(tib) {
  if ("asymp.LCL" %in% names(tib)) list(lcl = tib$asymp.LCL, ucl = tib$asymp.UCL)
  else if ("lower.CL" %in% names(tib)) list(lcl = tib$lower.CL, ucl = tib$upper.CL)
  else stop("No CI columns found in emmeans output")
}

# Categorical Time → numeric months for continuous x.
time_to_month <- c("BASELINE" = 0, "MO3" = 3, "MO6" = 6, "MO12" = 12)

# Full outcome names become the (left) y-axis titles.
sbp_lab <- "Systolic Blood Pressure (mmHg)"
dbp_lab <- "Diastolic Blood Pressure (mmHg)"

prep <- function(tib, outcome) {
  ci <- get_lcl_ucl(tib)
  tibble(
    Time    = as.character(tib$Time),
    stratum = factor(tib$BASELINE_AHA,
                     levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
                     labels = c("Elevated", "Stage 1", "Stage 2")),
    emmean  = tib$emmean,
    lcl     = ci$lcl,
    ucl     = ci$ucl,
    outcome = outcome
  ) %>%
    mutate(month = unname(time_to_month[Time]))
}

plot_df <- bind_rows(
  prep(emmeans_results$sbp_adj$emmeans, sbp_lab),
  prep(emmeans_results$dbp_adj$emmeans, dbp_lab)
) %>%
  mutate(outcome = factor(outcome, levels = c(sbp_lab, dbp_lab)))

# --- Significance indicators (α = 0.005), from the saved adjusted contrasts (Holm) ---
#   Points : filled = timepoint differs from BASELINE (within_stratum); open = not.
#   Segments: solid = significant change between adjacent timepoints (consecutive);
#            dashed = plateau (no significant change).
alpha <- 0.005
relabel_stratum <- function(x) factor(x,
  levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
  labels = c("Elevated", "Stage 1", "Stage 2"))

point_sig <- bind_rows(
  emmeans_results$sbp_adj$within_stratum %>% mutate(outcome = sbp_lab),
  emmeans_results$dbp_adj$within_stratum %>% mutate(outcome = dbp_lab)
) %>%
  transmute(outcome  = factor(outcome, levels = c(sbp_lab, dbp_lab)),
            stratum  = relabel_stratum(BASELINE_AHA),
            month    = unname(time_to_month[sub(" .*$", "", contrast)]),
            sig_base = p.holm < alpha) %>%
  filter(!is.na(month))

plot_df <- plot_df %>%
  left_join(point_sig, by = c("outcome", "stratum", "month")) %>%
  # Baseline is the reference (no change from itself) → open, like any non-significant point.
  mutate(pt_filled = dplyr::coalesce(sig_base, FALSE))

em_lookup <- plot_df %>% distinct(outcome, stratum, month, emmean)
seg_df <- bind_rows(
  emmeans_results$sbp_adj$consecutive %>% mutate(outcome = sbp_lab),
  emmeans_results$dbp_adj$consecutive %>% mutate(outcome = dbp_lab)
) %>%
  tidyr::separate(contrast, into = c("to", "from"), sep = " - ", remove = TRUE) %>%
  transmute(outcome = factor(outcome, levels = c(sbp_lab, dbp_lab)),
            stratum = relabel_stratum(BASELINE_AHA),
            x       = unname(time_to_month[from]),
            xend    = unname(time_to_month[to]),
            seg_sig = p.holm < alpha) %>%
  left_join(em_lookup, by = c("outcome", "stratum", "x" = "month")) %>% rename(y = emmean) %>%
  left_join(em_lookup, by = c("outcome", "stratum", "xend" = "month")) %>% rename(yend = emmean)

# AHA cutoff lines per outcome (facet-aware geom_hline). Lines only — no labels.
cutoff_df <- bind_rows(
  tibble(outcome = sbp_lab, cutoff = c(120, 130, 140)),
  tibble(outcome = dbp_lab, cutoff = c(80, 90))
) %>%
  mutate(outcome = factor(outcome, levels = c(sbp_lab, dbp_lab)))

# Per-panel y-limits (SBP 115–145, DBP 70–95) via invisible anchors —
# facet free_y can't take per-panel limits from a single scale.
limits_df <- bind_rows(
  tibble(outcome = sbp_lab, month = 0, y = c(115, 145)),
  tibble(outcome = dbp_lab, month = 0, y = c(70, 95))
) %>%
  mutate(outcome = factor(outcome, levels = c(sbp_lab, dbp_lab)))

# Stratum palette: Stage 1/2 AHA-chart hues; Elevated navy-grey (#44546A) — matched
# to Fig 1A / the abstract (recoloured from AHA gold per the 2026-06 review).
stratum_cols <- c("Elevated" = "#44546A", "Stage 1" = "#FFB600", "Stage 2" = "#BB3A03")

p <- ggplot(plot_df, aes(x = month, y = emmean, color = stratum, fill = stratum, group = stratum)) +
  geom_hline(data = cutoff_df, aes(yintercept = cutoff),
             color = "gray60", linetype = "dotted", linewidth = 0.7) +
  geom_blank(data = limits_df, aes(x = month, y = y), inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.18, color = NA) +
  geom_segment(data = seg_df, aes(x = x, xend = xend, y = y, yend = yend,
                                  color = stratum, linetype = seg_sig),
               linewidth = 0.9, inherit.aes = FALSE) +
  geom_point(aes(shape = pt_filled), size = 2.4) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.35, linewidth = 0.5) +
  facet_wrap(~ outcome, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_x_continuous(breaks = c(0, 3, 6, 12),
                     labels = c("Baseline", "3 Months", "6 Months", "12 Months")) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(breaks = scales::breaks_width(5), expand = expansion(mult = c(0, 0))) +
  scale_color_manual(values = stratum_cols, name = "AHA Baseline Stratum") +
  scale_fill_manual(values  = stratum_cols, name = "AHA Baseline Stratum") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     name = "vs. Baseline", labels = c(`TRUE` = "p < 0.005", `FALSE` = "p ≥ 0.005")) +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed"),
                        name = "Adjacent Change", labels = c(`TRUE` = "p < 0.005", `FALSE` = "p ≥ 0.005")) +
  guides(fill = "none", shape = "none", linetype = "none",  # significance keys → caption
         color = guide_legend(order = 1)) +
  labs(
    x = NULL, y = NULL,
    title    = "12-Month Blood Pressure Trajectories by Baseline AHA Classification",
    subtitle = "Adjusted LMM estimated marginal means with 95% CI"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid        = element_blank(),
    panel.spacing.y   = grid::unit(1.5, "lines"),
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.4),
    axis.text.x       = element_text(face = "bold"),
    strip.placement   = "outside",
    strip.background  = element_blank(),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 12),
    legend.position   = "bottom",
    plot.title        = element_text(face = "bold", size = 12.5)
  )

out_png <- file.path(fig_dir, "5_5_figure_2.png")
ggsave(out_png, p, width = 7, height = 8, dpi = 300)

cat("Wrote:\n  ", out_png, "\n")
cat("\nCell estimates (adjusted) used in figure:\n")
print(plot_df %>% select(outcome, stratum, Time, month, emmean, lcl, ucl) %>%
        arrange(outcome, stratum, month))
