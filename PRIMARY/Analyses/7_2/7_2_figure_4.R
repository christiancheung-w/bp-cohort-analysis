# 7_2_figure_4.R
# SAP §7.2 Figure 4: forest plot of standardized W/B-decomposed daily-signal
# coefficients (β per SD of predictor) predicting each reading's (BP − baseline).
#
# Mirrors the Figure 5 / §6.2 forest grammar (6_2_logistic.R): a SINGLE panel
# with the two outcomes (SBP, DBP) overlaid as vertically-dodged coloured points
# (SBP blue / DBP green), and the between-/within-person decomposition drawn as
# a nested bracket OUTSIDE the y-axis (legendry::guide_axis_nested). The one
# necessary divergence from §6.2: the x-axis is the linear β scale (mmHg per SD),
# symmetric about a dashed 0 reference — NOT the log-OR-about-1.0 scale.
#
# Kept separate from 7_2_wb_lmm.R (which fits the cluster-robust LMMs and writes
# the coefficient CSV) so the figure can be restyled without refitting — the same
# figure-script separation used in §5.1 (5_1_figure_1a.R) and §5.5 (5_5_figure_2.R).
#
# Input:  Analyses/7_2/derived/results/7_2_coefficients.csv
# Output: Analyses/7_2/derived/figures/7_2_figure_4.png

library(tidyverse)
suppressPackageStartupMessages(library(legendry))

project_dir <- "./PRIMARY"
in_path  <- file.path(project_dir, "Analyses", "7_2", "derived", "results", "7_2_coefficients.csv")
fig_dir  <- file.path(project_dir, "Analyses", "7_2", "derived", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

coefs_all <- read_csv(in_path, show_col_types = FALSE)

# The 10 W/B fixed effects (5 metrics × {between (_b), within (_w)}).
wb_terms <- c(outer(c("RHR", "HRV", "Strain", "SleepDuration", "SleepConsistency"),
                    c("_b", "_w"), paste0))

# Display labels matched to §6.2 / §7.1: Strain → "Physical Activity",
# clusters ordered RHR → HRV → Sleep Duration → Sleep Consistency → Physical Activity.
metric_labels <- c(RHR = "RHR", HRV = "HRV",
                   SleepDuration = "Sleep Duration",
                   SleepConsistency = "Sleep Consistency",
                   Strain = "Physical Activity")
metric_order  <- c("RHR", "HRV", "SleepDuration", "SleepConsistency", "Strain")

# Between-Person cluster on top, within-person below, with a non-breaking-space
# spacer row between them (a blank gap; avoids an "NA" tick) — mirrors §6.2.
order_top_to_bottom <- c(paste(metric_labels[metric_order], "Between-Person", sep = "&"),
                         " ",
                         paste(metric_labels[metric_order], "Within-Person",  sep = "&"))
real_levels <- c(paste(metric_labels[metric_order], "Between-Person", sep = "&"),
                 paste(metric_labels[metric_order], "Within-Person",  sep = "&"))

forest_df <- coefs_all %>%
  filter(term %in% wb_terms) %>%
  mutate(
    metric    = sub("_(b|w)$", "", term),
    form      = ifelse(grepl("_b$", term), "Between-Person", "Within-Person"),
    term_label = factor(paste(metric_labels[metric], form, sep = "&"),
                        levels = rev(order_top_to_bottom)),
    # SBP listed last → position_dodge draws it at the higher (top) y position;
    # legend order restored to SBP-first via scale_color_manual(breaks=).
    outcome   = factor(outcome, levels = c("DBP-baseline", "SBP-baseline")),
    sig_label = ifelse(p_value < 0.005, "p < 0.005", "p ≥ 0.005")
  )

outcome_cols <- c("SBP-baseline" = "#00538A", "DBP-baseline" = "#5B9BD5")

p <- ggplot(forest_df,
            aes(x = estimate, y = term_label,
                color = outcome, shape = sig_label)) +
  geom_vline(xintercept = 0, color = "gray60", linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0.4,
                linewidth = 0.5, orientation = "y",
                position = position_dodge(width = 0.6)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_x_continuous(limits = c(-0.7, 0.7), breaks = seq(-0.6, 0.6, by = 0.2),
                     labels = c("-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "0.6"),
                     expand = expansion(0)) +   # plot area capped ±0.7, ticks flush to the box
  scale_y_discrete(drop = FALSE,
                   guide = guide_axis_nested(
                     key = key_range_auto(sep = "&"),
                     regular_key = key_manual(aesthetic = real_levels,
                                              label = sub("&.*$", "", real_levels)),
                     bracket = "square",
                     levels_text = list(
                       element_text(face = "plain", angle = 0),                  # inner metric tick labels
                       element_text(face = "bold",  angle = 90, hjust = 0.5)))) + # centred Between/Within titles
  scale_color_manual(values = outcome_cols, name = NULL,
                     breaks = c("SBP-baseline", "DBP-baseline"),
                     labels = c("SBP-baseline" = "Systolic Blood Pressure",
                                "DBP-baseline" = "Diastolic Blood Pressure")) +
  scale_shape_manual(values = c("p < 0.005" = 16, "p ≥ 0.005" = 1), guide = "none") +
  labs(
    x = "β (mmHg from Baseline per SD of Predictor, 95% CI)",
    y = NULL,
    title    = "Within- and Between-Person Wearable Biometric Associations With Blood Pressure",
    subtitle = "Baseline-Anchored Within/Between-Person Decomposed Linear Mixed Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position    = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.4),
    axis.title.x       = element_text(face = "bold", margin = margin(t = 10)),
    axis.text.y        = element_text(face = "plain"),
    plot.title         = element_text(face = "bold", size = 11)
  )

out_png <- file.path(fig_dir, "7_2_figure_4.png")
suppressWarnings(  # benign key_range_auto split note from the blank spacer row
  ggsave(out_png, p, width = 9.5, height = 7, dpi = 300))
cat("Wrote forest plot:\n  ", out_png, "\n")

cat("\nβ estimates used in figure:\n")
print(forest_df %>%
        select(outcome, form, metric, estimate, conf_low, conf_high, sig_label) %>%
        arrange(outcome, form, metric))
