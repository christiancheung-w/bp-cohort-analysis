# abstract_forests.R
# Targeted 3-panel abstract figure — Panels B and C only.
#   Panel A is the locked Figure 1A (5_1/derived/figures/5_1_figure_1a.png); reused as-is.
#   Panel B = WITHIN-person subset of §7.2 Figure 4: RHR + Sleep Duration, SBP & DBP,
#             β (mmHg from baseline per SD), dashed 0 reference.
#   Panel C = OR forest splicing TWO models. Predictor names match the §8 / §6.2
#             figure style ("Daily App Time", "RHR"); the two models' N is carried
#             in the subtitle so the composite is not read as one model:
#               - Daily App Time (log) from §8 joint primary  (n = 7,576)
#               - RHR, 6-mo study mean, from §6.2             (n = 14,189)
#             Odds ratios on a log scale, dashed 1.0 reference.
#
# Mirrors the house forest grammar of 7_2_figure_4.R / 8_figure_forest.R:
#   SBP blue (#1f78b4) / DBP green (#33a02c), vertically dodged; filled marker =
#   significant at alpha = 0.005; boxed panel, no y-gridlines, bold x-title,
#   bottom legend. Reads only frozen coefficient CSVs; writes only PNGs.
#
# Run with working directory = dID_Analysis/  (same relative-path scheme as the
# other figure scripts: project_dir <- "./PRIMARY").
#
# Inputs:  Analyses/7_2/derived/results/7_2_coefficients.csv
#          Analyses/6_2/derived/results/6_2_coefficients_primary.csv
#          Analyses/8_engagement/derived/results/8_coefficients_joint_primary.csv
# Outputs: Analyses/Abstract/derived/figures/abstract_panelB_within.png
#          Analyses/Abstract/derived/figures/abstract_panelC_or.png

library(tidyverse)

project_dir <- "./PRIMARY"
fig_dir <- file.path(project_dir, "Analyses", "Abstract", "derived", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

outcome_cols <- c("SBP" = "#1f78b4", "DBP" = "#33a02c")

# Shared theme so B and C read as a matched pair.
forest_theme <- theme_minimal(base_size = 11) +
  theme(
    legend.position    = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.4),
    axis.title.x       = element_text(face = "bold", margin = margin(t = 10)),
    plot.title         = element_text(face = "bold", size = 11),
    # Top margin opens a band above the title for the panel tag (drawn below via cowplot).
    plot.margin        = margin(t = 30, r = 5, b = 5, l = 5)
  )

# Panel tag helper: draw the letter on the full canvas with cowplot (reliable 0-1 device
# coordinates), on its own line above the title. y places the letter ~62 px from the top
# at 300 dpi over the 3.6-inch-tall figure — matched (absolute) to Panel A. x is passed
# per panel so the letter sits over the middle of the predictor (y-axis) labels.
add_panel_tag <- function(plot, letter, x) {
  cowplot::ggdraw(plot) +
    cowplot::draw_label(letter, x = x, y = 1 - 62 / (3.6 * 300),
                        fontface = "bold", size = 18, hjust = 0.5, vjust = 0.5)
}

# ----------------------------------------------------------------------------
# Panel B — §7.2 within-person, RHR + Sleep Duration only (β per SD)
# ----------------------------------------------------------------------------
b_in <- read_csv(file.path(project_dir, "Analyses", "7_2", "derived", "results",
                           "7_2_coefficients.csv"), show_col_types = FALSE)

b_labels <- c(RHR_w = "RHR", SleepDuration_w = "Sleep Duration")   # match §7.2/§6.2 metric labels
b_order  <- c("RHR", "Sleep Duration")   # y-axis top -> bottom

panelB_df <- b_in %>%
  filter(term %in% names(b_labels)) %>%
  mutate(
    term_label = factor(b_labels[term], levels = rev(b_order)),
    outcome    = factor(sub("-baseline$", "", outcome), levels = c("DBP", "SBP")), # SBP last -> drawn on top
    sig_label  = ifelse(p_value < 0.005, "p < 0.005", "p >= 0.005")
  )

pB <- ggplot(panelB_df,
             aes(x = estimate, y = term_label, color = outcome, shape = sig_label)) +
  geom_vline(xintercept = 0, color = "gray60", linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0.35, linewidth = 0.5,
                orientation = "y", position = position_dodge(width = 0.6)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_x_continuous(limits = c(-0.4, 0.8), breaks = seq(-0.4, 0.8, 0.2),
                     expand = expansion(0)) +
  scale_color_manual(values = outcome_cols, name = NULL, breaks = c("SBP", "DBP"),
                     labels = c(SBP = "Systolic Blood Pressure",
                                DBP = "Diastolic Blood Pressure")) +
  scale_shape_manual(values = c("p < 0.005" = 16, "p >= 0.005" = 1), guide = "none") +
  labs(x = "β (mmHg from Baseline per SD of Predictor, 95% CI)", y = NULL,
       title    = "Within-Person Wearable Biometrics vs. Same-Day Blood Pressure",
       subtitle = "Baseline-Anchored Within-Person Linear Mixed Model") +
  forest_theme

ggsave(file.path(fig_dir, "abstract_panelB_within.png"), add_panel_tag(pB, "B", 0.063),
       width = 8, height = 3.6, dpi = 300)

# ----------------------------------------------------------------------------
# Panel C — responder ORs: Daily App Time (§8 joint) + RHR 6-mo (§6.2)
# Two different models; predictor labels match §8/§6.2 style, the two N in the subtitle.
# ----------------------------------------------------------------------------
c62 <- read_csv(file.path(project_dir, "Analyses", "6_2", "derived", "results",
                          "6_2_coefficients_primary.csv"), show_col_types = FALSE) %>%
  filter(term == "RHR_STUDY_MEAN") %>%
  transmute(outcome = sub("_responder$", "", outcome),
            predictor = "Δ RHR",
            OR, OR_lcl, OR_ucl, sig_005)

c8 <- read_csv(file.path(project_dir, "Analyses", "8_engagement", "derived", "results",
                         "8_coefficients_joint_primary.csv"), show_col_types = FALSE) %>%
  filter(term == "LOG_AVG_DAILY_TIME_IN_APP_MIN") %>%
  transmute(outcome,
            predictor = "Daily App Time",
            OR, OR_lcl, OR_ucl, sig_005)

c_order <- c("Δ RHR", "Daily App Time")

panelC_df <- bind_rows(c62, c8) %>%
  mutate(
    term_label = factor(predictor, levels = rev(c_order)),
    outcome    = factor(outcome, levels = c("DBP", "SBP")),
    sig_005    = as.logical(sig_005)
  )

pC <- ggplot(panelC_df,
             aes(x = OR, y = term_label, color = outcome, shape = sig_005)) +
  geom_vline(xintercept = 1, color = "gray60", linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(xmin = OR_lcl, xmax = OR_ucl), width = 0.35, linewidth = 0.5,
                orientation = "y", position = position_dodge(width = 0.6)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_x_log10(breaks = seq(0.9, 1.4, 0.1), expand = expansion(mult = 0.03)) +
  scale_color_manual(values = outcome_cols, name = NULL, breaks = c("SBP", "DBP"),
                     labels = c(SBP = "SBP-Responder", DBP = "DBP-Responder")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(x = "Adjusted OR (Log Scale, 95% CI)", y = NULL,
       title    = "Predictors of Clinically Meaningful Blood Pressure Reduction",
       subtitle = "Adjusted Logistic Regression — 6-Month Responder Status") +
  forest_theme

ggsave(file.path(fig_dir, "abstract_panelC_or.png"), add_panel_tag(pC, "C", 0.063),
       width = 8, height = 3.6, dpi = 300)

# --- Echo the estimates used, for verification ---
cat("\nPanel B (beta per SD):\n")
print(panelB_df %>% select(outcome, term_label, estimate, conf_low, conf_high, sig_label) %>%
        arrange(term_label, outcome))
cat("\nPanel C (OR):\n")
print(panelC_df %>% select(outcome, term_label, OR, OR_lcl, OR_ucl, sig_005) %>%
        arrange(term_label, outcome))
cat("\nWrote:\n  ", file.path(fig_dir, "abstract_panelB_within.png"),
    "\n  ", file.path(fig_dir, "abstract_panelC_or.png"), "\n")
