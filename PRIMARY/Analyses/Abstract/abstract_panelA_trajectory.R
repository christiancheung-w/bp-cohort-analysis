# abstract_panelA_trajectory.R
# Abstract Panel A: abstract-specific copy of SAP §5.1 Figure 1A.
#
# Identical to Analyses/5_1/5_1_figure_1a.R in every respect EXCEPT one visual
# change for the abstract: the "Elevated" stratum is recoloured from the AHA gold
# (#EDCE44) to a navy grey (#44546A), so it reads as visually distinct from the
# warm Stage 1 / Stage 2 colours. All geometry, cutoff lines, significance
# encoding, theme, and figure dimensions are unchanged.
#
# Reads the same frozen §5.1 emmeans output as the locked figure; writes to the
# Abstract figures folder so the locked 5_1_figure_1a.png is NOT overwritten.
#
# Input:  Analyses/5_1/derived/emmeans_results.rds   (frozen; shared with 5_1_figure_1a.R)
# Output: Analyses/Abstract/derived/figures/abstract_panelA_trajectory.png

library(tidyverse)

project_dir <- "./PRIMARY"
in_path     <- file.path(project_dir, "Analyses", "5_1", "derived", "emmeans_results.rds")
fig_dir     <- file.path(project_dir, "Analyses", "Abstract", "derived", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

emmeans_results <- read_rds(in_path)

# Categorical Time → numeric months for continuous x.
time_to_month <- c("BASELINE" = 0, "MO3" = 3, "MO6" = 6)

# Full outcome names become the (left) y-axis titles.
sbp_lab <- "Systolic Blood Pressure (mmHg)"
dbp_lab <- "Diastolic Blood Pressure (mmHg)"

prep <- function(tib, outcome) {
  tib %>%
    mutate(
      outcome = outcome,
      month   = unname(time_to_month[as.character(Time)]),
      stratum = factor(BASELINE_AHA,
                       levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
                       labels = c("Elevated", "Stage 1", "Stage 2"))
    )
}

plot_df <- bind_rows(
  prep(emmeans_results$sbp_unadj$emmeans, sbp_lab),
  prep(emmeans_results$dbp_unadj$emmeans, dbp_lab)
) %>%
  mutate(outcome = factor(outcome, levels = c(sbp_lab, dbp_lab)))

# --- Significance indicators (α = 0.005), from the saved unadjusted contrasts ---
#   Points : filled = timepoint differs from BASELINE (within_stratum / trt.vs.ctrl);
#            open = not. Baseline point is the reference (drawn filled).
#   Segments: solid = significant change between adjacent timepoints (consecutive);
#            dashed = no significant change (plateau).
alpha <- 0.005
relabel_stratum <- function(x) factor(x,
  levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
  labels = c("Elevated", "Stage 1", "Stage 2"))
month_of <- c("BASELINE" = 0, "MO3" = 3, "MO6" = 6, "MO12" = 12)

# Point significance vs baseline
point_sig <- bind_rows(
  emmeans_results$sbp_unadj$within_stratum %>% mutate(outcome = sbp_lab),
  emmeans_results$dbp_unadj$within_stratum %>% mutate(outcome = dbp_lab)
) %>%
  transmute(outcome  = factor(outcome, levels = c(sbp_lab, dbp_lab)),
            stratum  = relabel_stratum(BASELINE_AHA),
            month    = unname(month_of[sub(" .*$", "", contrast)]),  # left token of "MOx - BASELINE"
            sig_base = p.value < alpha) %>%
  filter(!is.na(month))

plot_df <- plot_df %>%
  left_join(point_sig, by = c("outcome", "stratum", "month")) %>%
  # Baseline is the reference (no change from itself) → open, like any non-significant point.
  mutate(pt_filled = dplyr::coalesce(sig_base, FALSE))

# Segment significance between adjacent timepoints
em_lookup <- plot_df %>% distinct(outcome, stratum, month, emmean)
seg_df <- bind_rows(
  emmeans_results$sbp_unadj$consecutive %>% mutate(outcome = sbp_lab),
  emmeans_results$dbp_unadj$consecutive %>% mutate(outcome = dbp_lab)
) %>%
  tidyr::separate(contrast, into = c("to", "from"), sep = " - ", remove = TRUE) %>%
  transmute(outcome = factor(outcome, levels = c(sbp_lab, dbp_lab)),
            stratum = relabel_stratum(BASELINE_AHA),
            x       = unname(month_of[from]),
            xend    = unname(month_of[to]),
            seg_sig = p.value < alpha) %>%
  left_join(em_lookup, by = c("outcome", "stratum", "x" = "month")) %>% rename(y = emmean) %>%
  left_join(em_lookup, by = c("outcome", "stratum", "xend" = "month")) %>% rename(yend = emmean)

# AHA cutoff lines per outcome (facet-aware geom_hline).
cutoff_df <- bind_rows(
  tibble(outcome = sbp_lab, cutoff = c(120, 130, 140),
         label = c("Elevated (≥120)", "Stage 1 (≥130)", "Stage 2 (≥140)")),
  tibble(outcome = dbp_lab, cutoff = c(80, 90),
         label = c("Stage 1 (≥80)", "Stage 2 (≥90)"))
) %>%
  mutate(outcome = factor(outcome, levels = c(sbp_lab, dbp_lab)))

# Per-panel y-limits (SBP 120–145, DBP 70–95) via invisible anchors —
# facet free_y can't take per-panel limits from a single scale.
limits_df <- bind_rows(
  tibble(outcome = sbp_lab, month = 0, y = c(120, 145)),
  tibble(outcome = dbp_lab, month = 0, y = c(70, 95))
) %>%
  mutate(outcome = factor(outcome, levels = c(sbp_lab, dbp_lab)))

# Stratum palette: Stage 1 / Stage 2 sampled from the AHA Blood Pressure Categories
# chart. ABSTRACT-SPECIFIC CHANGE: "Elevated" recoloured from AHA gold (#EDCE44) to
# navy grey (#44546A) so it is visually distinct from the warm Stage 1 / Stage 2 hues.
stratum_cols <- c("Elevated" = "#44546A", "Stage 1" = "#FFB600", "Stage 2" = "#BB3A03")

p <- ggplot(plot_df, aes(x = month, y = emmean, color = stratum, fill = stratum, group = stratum)) +
  geom_hline(data = cutoff_df, aes(yintercept = cutoff),
             color = "gray60", linetype = "dotted", linewidth = 0.7) +
  geom_blank(data = limits_df, aes(x = month, y = y), inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = 0.18, color = NA) +
  geom_segment(data = seg_df, aes(x = x, xend = xend, y = y, yend = yend,
                                  color = stratum, linetype = seg_sig),
               linewidth = 0.9, inherit.aes = FALSE) +
  geom_point(aes(shape = pt_filled), size = 2.4) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.35, linewidth = 0.5) +
  geom_text(data = dplyr::filter(cutoff_df, !grepl("Elevated", label)),
            aes(y = cutoff, label = label), x = Inf,
            hjust = 1, vjust = -0.7, size = 4, fontface = "bold", color = "black", inherit.aes = FALSE) +
  geom_text(data = dplyr::filter(cutoff_df, grepl("Elevated", label)),
            aes(y = cutoff, label = label), x = Inf,
            hjust = 1, vjust = 1.6, size = 4, fontface = "bold", color = "black", inherit.aes = FALSE) +
  facet_wrap(~ outcome, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_x_continuous(breaks = c(0, 3, 6), labels = c("Baseline", "3 Months", "6 Months")) +
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
    title    = "6-Month Blood Pressure Trajectories by Baseline AHA Classification",
    subtitle = "Unadjusted LMM estimated means with 95% CI"
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
    plot.title        = element_text(face = "bold", size = 12.5),
    # Top margin opens a band above the title for the panel tag (drawn below via cowplot).
    plot.margin       = margin(t = 30, r = 5, b = 5, l = 5)
  )

# Panel tag "A": drawn on the full canvas (cowplot draw_label uses reliable 0-1 device
# coordinates), on its own line above the title, horizontally in line with the y-axis
# (panel left edge). y is set so the letter sits ~62 px from the top at 300 dpi, matched
# (absolute) to Panels B/C; x = 0.135 aligns with the panel left edge.
p <- cowplot::ggdraw(p) +
  cowplot::draw_label("A", x = 0.102, y = 1 - 62 / (8 * 300),
                      fontface = "bold", size = 18, hjust = 0.5, vjust = 0.5)

out_png <- file.path(fig_dir, "abstract_panelA_trajectory.png")
ggsave(out_png, p, width = 7, height = 8, dpi = 300)

cat("Wrote:\n  ", out_png, "\n")
