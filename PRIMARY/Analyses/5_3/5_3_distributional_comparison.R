# 5_3_distributional_comparison.R
# SAP §5.3: cohort-overall distributional comparison of BP at baseline vs MO6.
#
# Overlaid pre/post density plot, cohort-overall (no AHA stratification),
# separately for SBP and DBP. Supplementary numeric: cohort-overall mean Δ
# in mmHg + standardized effect size d = mean(Δ) / pooled baseline SD of the
# outcome (baseline-SD standardization; Kazis et al., Med Care 1989).
# Harmonized with §5.1 Table 2 / §5.5 Table S3 (same denominator); numerator
# is the empirical cohort-overall mean change (this section is model-free).
#
# Per SAP §5.3: no formal hypothesis test. Wilson CIs not applicable (continuous
# outcome). 95% CI on mean Δ added for transparency.
#
# Cohort: §2.2 hypertensive (Elevated + Stage 1 + Stage 2) with sex recorded
# (matches §5.1/§5.5). Users without both baseline and MO6 for the given
# outcome are excluded from that outcome's Δ summary (baseline SD uses the full
# cohort baseline, matching §5.1/§5.5).
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/5_3/derived/):
#   results/5_3_summary_stats.csv         — outcome × {n, mean(Δ), SD(Δ), 95% CI, d}
#   results/5_3_paper.md                  — paper-ready text
#   results/5_3_distributional_comparison.log
#   figures/5_3_figure_1c.png             — overlaid baseline/MO6 densities, 2-panel

library(tidyverse)
library(ggtext)   # geom_richtext(): markdown in annotations (italic d)

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "5_3", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "5_3_distributional_comparison.log")
sink(log_path, split = TRUE)
cat("§5.3 distributional comparison —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read frame and filter to §2.2 cohort ---
user_frame <- read_rds(input_path)
cohort <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female"))
cat("§2.2 cohort (sex-restricted, matches §5.1/§5.5):", nrow(cohort), "users\n\n")

# --- Cohort baseline SD: denominator for d (one value per outcome) ---
# Pooled across strata; full-cohort baseline (matches §5.1/§5.5). Kazis (1989).
baseline_sd <- cohort %>%
  summarise(SBP = sd(BASELINE_SBP, na.rm = TRUE),
            DBP = sd(BASELINE_DBP, na.rm = TRUE)) %>%
  unlist()
cat("Cohort baseline SD (d denominator):\n")
cat("  SBP:", round(baseline_sd[["SBP"]], 4), " DBP:", round(baseline_sd[["DBP"]], 4), "\n\n")

# --- Long format for density plot (one row per user × outcome × timepoint) ---
density_long <- cohort %>%
  select(STUDY_ID, BASELINE_SBP, MO6_SBP, BASELINE_DBP, MO6_DBP) %>%
  pivot_longer(-STUDY_ID,
               names_to  = c("timepoint", "outcome"),
               names_sep = "_",
               values_to = "value") %>%
  mutate(timepoint = factor(timepoint, levels = c("BASELINE", "MO6"),
                            labels = c("Baseline", "6 Months")),
         outcome   = factor(outcome,   levels = c("SBP", "DBP"),
                            labels = c("SBP (mmHg)", "DBP (mmHg)"))) %>%
  filter(!is.na(value))

cat("Density-plot inputs per (outcome × timepoint):\n")
print(density_long %>% count(outcome, timepoint))

# --- Per-user paired Δs for summary stats ---
paired <- cohort %>%
  transmute(
    STUDY_ID,
    SBP = MO6_SBP - BASELINE_SBP,
    DBP = MO6_DBP - BASELINE_DBP
  ) %>%
  pivot_longer(-STUDY_ID, names_to = "outcome", values_to = "delta") %>%
  filter(!is.na(delta))

summary_stats <- paired %>%
  group_by(outcome) %>%
  summarise(
    n            = n(),
    mean_delta   = mean(delta),
    sd_delta     = sd(delta),
    se_delta     = sd_delta / sqrt(n),
    ci_lcl       = mean_delta - 1.959964 * se_delta,
    ci_ucl       = mean_delta + 1.959964 * se_delta,
    .groups = "drop"
  ) %>%
  mutate(d_baseline = mean_delta / unname(baseline_sd[outcome]))

write_csv(summary_stats, file.path(results_dir, "5_3_summary_stats.csv"))

cat("\n--- Cohort-overall paired-difference summary ---\n")
print(summary_stats)

# --- Figure 1C: overlaid densities with AHA cutoff lines ---
cutoff_df <- bind_rows(
  tibble(outcome = "SBP (mmHg)", cutoff = c(120, 130, 140), label = c("Elevated", "Stage 1", "Stage 2")),
  tibble(outcome = "DBP (mmHg)", cutoff = c(80, 90),        label = c("Stage 1", "Stage 2"))
) %>%
  mutate(outcome = factor(outcome, levels = c("SBP (mmHg)", "DBP (mmHg)")))

tp_cols <- c("Baseline" = "#9E9E9E", "6 Months" = "#2A9D8F")

# Label tibble for in-panel mean Δ + d annotation
label_df <- summary_stats %>%
  mutate(outcome = factor(ifelse(outcome == "SBP", "SBP (mmHg)", "DBP (mmHg)"),
                          levels = c("SBP (mmHg)", "DBP (mmHg)")),
         lbl = sprintf("Mean Δ = %+.2f mmHg<br>95%% CI [%+.2f, %+.2f]<br>*d* = %.3f<br>n = %s",
                       mean_delta, ci_lcl, ci_ucl, abs(d_baseline),
                       format(n, big.mark = ",")))

# Full outcome names become the per-panel x-axis titles (facet strip at bottom).
outcome_titles <- c("SBP (mmHg)" = "Systolic Blood Pressure (mmHg)",
                    "DBP (mmHg)" = "Diastolic Blood Pressure (mmHg)")

# Extend the DBP panel's x-axis down to 50 so the lower tail isn't curtailed.
x_floor <- tibble(outcome = factor("DBP (mmHg)", levels = c("SBP (mmHg)", "DBP (mmHg)")),
                  value = 50)

p <- ggplot(density_long, aes(x = value, fill = timepoint, color = timepoint)) +
  geom_density(alpha = 0.35, linewidth = 0.6, adjust = 1.0) +
  geom_blank(data = x_floor, aes(x = value), inherit.aes = FALSE) +
  geom_vline(data = cutoff_df, aes(xintercept = cutoff),
             color = "gray30", linetype = "dotted", linewidth = 0.4) +
  geom_text(data = cutoff_df, aes(x = cutoff, label = label), y = Inf,
            hjust = -0.1, vjust = 1.4, size = 2.8, color = "gray25", inherit.aes = FALSE) +
  facet_wrap(~ outcome, ncol = 1, scales = "free", strip.position = "bottom",
             labeller = as_labeller(outcome_titles)) +
  geom_richtext(data = label_df,
            aes(label = lbl, x = Inf, y = Inf),
            hjust = 1, vjust = 1, inherit.aes = FALSE,
            size = 3.4, lineheight = 1.1, fontface = "bold", color = "gray20",
            fill = NA, label.color = NA, label.padding = grid::unit(c(4, 5, 4, 5), "pt")) +
  scale_fill_manual(values = tp_cols, name = "Timepoint") +
  scale_color_manual(values = tp_cols, name = "Timepoint") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = "Density (per mmHg)",
    title = "Baseline-to-6-Month Blood Pressure Distribution Shift"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    panel.spacing.y  = grid::unit(1.4, "lines"),
    axis.line        = element_line(color = "black", linewidth = 0.4),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    axis.title.y     = element_text(face = "bold", size = 12),
    strip.placement  = "outside",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 12),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold", size = 12),
    legend.text      = element_text(size = 11),
    plot.title       = element_text(face = "bold", size = 12.5, margin = margin(b = 12))
  )

out_png <- file.path(figures_dir, "5_3_figure_1c.png")
ggsave(out_png, p, width = 7, height = 8, dpi = 300)
cat("\nWrote figure:\n  ", out_png, "\n")

# --- Paper-ready markdown ---
fmt_signed <- function(x, digits = 2) sprintf("%+.*f", digits, x)

md_lines <- c(
  "# §5.3 — Cohort-overall distributional comparison",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "## Supplementary text (paper-quote-ready)",
  ""
)
for (otc in c("SBP", "DBP")) {
  s <- summary_stats %>% filter(outcome == otc)
  md_lines <- c(md_lines, sprintf(
    "- **%s** (n = %s): mean Δ = %s mmHg (95%% CI [%s, %s]), SD(Δ) = %.2f mmHg, *d* = %s.",
    otc, format(s$n, big.mark = ","),
    fmt_signed(s$mean_delta), fmt_signed(s$ci_lcl), fmt_signed(s$ci_ucl),
    s$sd_delta, fmt_signed(s$d_baseline, 3)
  ))
}
md_lines <- c(md_lines,
  "",
  "## Notes",
  "",
  "- Cohort = §2.2 hypertensive (Elevated + Stage 1 + Stage 2). Cohort-overall means *no* stratification by AHA — these statistics aggregate across all three strata weighted by their cohort proportions.",
  "- Per SAP §5.3, no formal hypothesis test. Effect-size frame is mean Δ in mmHg + standardized *d* (supplementary).",
  "- The Stage 2 stratum is overrepresented in the absolute Δ contribution relative to its cohort share, because its per-user Δ magnitude is largest. Per-stratum effect sizes are reported in §5.1 Table 2; this section presents the cohort-level mean for distributional-shift reporting only.",
  "- *d* here = (cohort-overall mean of paired MO6 − Baseline differences) / (pooled baseline SD of the outcome, all strata combined; Kazis et al., Med Care 1989). Same baseline-SD standardizer as §5.1 Table 2 / §5.5 Table S3; the numerator is the empirical cohort-overall mean change (this section is model-free), whereas §5.1/§5.5 use the adjusted LMM Δ per stratum. Signed (negative = BP reduction).",
  "- Pairs with: §5.1 trajectory analysis (provides per-stratum Δ), Figure 1C (this artifact)."
)
writeLines(md_lines, file.path(results_dir, "5_3_paper.md"))
cat("\nWrote markdown:\n  ", file.path(results_dir, "5_3_paper.md"), "\n")

cat("\n=== §5.3 distributional comparison done ===\n")
sink()
