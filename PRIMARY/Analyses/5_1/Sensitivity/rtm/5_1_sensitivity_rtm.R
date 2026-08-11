# 5_1_sensitivity_rtm.R
# SAP §5.1 sensitivity: regression-to-the-mean (RTM) adjustment
# (Oldham 1962 / Blomqvist 1977 framing).
#
# Approach: fit lm(Δ ~ baseline_BP) on the §2.2 cohort. The slope β captures
# the baseline-Δ correlation; the magnitude of β × (X̄_stratum − X̄_cohort)
# represents the share of each stratum's apparent mean Δ that is attributable
# to RTM. Report per-stratum decomposition:
#
#   Δ_observed_s     = mean(Δ) within stratum s
#   Δ_RTM_attrib_s   = β × (X̄_s − X̄_cohort)
#   Δ_adjusted_s     = Δ_observed_s − Δ_RTM_attrib_s
#   pct_RTM_attrib_s = Δ_RTM_attrib_s / Δ_observed_s
#
# Interpretation note (load-bearing for the manuscript):
# The Oldham/Blomqvist approach attributes ALL baseline-Δ correlation to RTM.
# If a real dose-response of intervention exists (e.g., wearable engagement
# works proportionally better at higher baseline BP), that real signal will
# be falsely absorbed into the RTM estimate. The "RTM-adjusted" Δs reported
# here are therefore CONSERVATIVE LOWER BOUNDS on the true effect — they
# over-correct in the presence of any real dose-response. Frame in the
# manuscript as "RTM-conservative bounds," not "RTM-corrected."
#
# Cohort: §2.2 hypertensive (Elevated + Stage 1 + Stage 2), matching §5.1.
# Each user contributes one Δ per outcome (MO6 − BASELINE).
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/5_1/Sensitivity/rtm/derived/):
#   results/regression_coefficients.csv     — β, α, R², ρ_test_retest, n
#   results/per_stratum_decomposition.csv   — per-stratum: observed, RTM-attrib, adjusted, pct
#   results/cross_stratum_decomposition.csv — RTM-adjusted cross-stratum Δ-difference
#   results/5_1_sensitivity_rtm_paper.md
#   results/5_1_sensitivity_rtm.log
#   figures/5_1_sensitivity_rtm_scatter.png — Δ vs baseline_BP, stratum-colored, with regression line

library(tidyverse)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")  # emit_sizes() (no LMM/contrast families here)

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "5_1", "Sensitivity", "rtm", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "5_1_sensitivity_rtm.log")
sink(log_path, split = TRUE)
cat("§5.1 sensitivity — Oldham/Blomqvist RTM adjustment —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

user_frame <- read_rds(input_path)
cohort <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female")) %>%
  transmute(
    STUDY_ID,
    stratum     = factor(BASELINE_AHA,
                         levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
                         labels = c("Elevated", "Stage 1", "Stage 2")),
    BASELINE_SBP,
    BASELINE_DBP,
    DELTA_SBP_MO6 = MO6_SBP - BASELINE_SBP,
    DELTA_DBP_MO6 = MO6_DBP - BASELINE_DBP,
    MO6_SBP,
    MO6_DBP
  )
cat("§2.2 cohort:", nrow(cohort), "users\n\n")

# Capture sizes now (the `cohort` symbol is reused downstream) for Table A.2.
.sz_sbp <- sum(!is.na(cohort$BASELINE_SBP) & !is.na(cohort$MO6_SBP))
.sz_dbp <- sum(!is.na(cohort$BASELINE_DBP) & !is.na(cohort$MO6_DBP))
.sz_ceiling <- nrow(cohort)

# --- Fit Oldham regression for each outcome ---
fit_rtm <- function(outcome_label, baseline_col, delta_col, mo6_col, data) {
  d <- data %>% select(stratum,
                       baseline = all_of(baseline_col),
                       delta    = all_of(delta_col),
                       mo6      = all_of(mo6_col)) %>%
    filter(!is.na(baseline), !is.na(delta))

  n  <- nrow(d)
  x_bar <- mean(d$baseline)
  rho   <- cor(d$baseline, d$mo6)  # test-retest correlation

  fit <- lm(delta ~ baseline, data = d)
  beta  <- coef(fit)["baseline"]
  alpha <- coef(fit)["(Intercept)"]
  r2    <- summary(fit)$r.squared
  beta_se <- summary(fit)$coefficients["baseline", "Std. Error"]

  # Per-stratum decomposition
  per_stratum <- d %>%
    group_by(stratum) %>%
    summarise(
      n_users          = n(),
      mean_baseline    = mean(baseline),
      mean_delta_obs   = mean(delta),
      sd_delta         = sd(delta),
      .groups = "drop"
    ) %>%
    mutate(
      x_centered       = mean_baseline - x_bar,
      delta_rtm_attrib = beta * x_centered,
      delta_adjusted   = mean_delta_obs - delta_rtm_attrib,
      pct_rtm_attrib   = ifelse(abs(mean_delta_obs) > 0.5,
                                100 * delta_rtm_attrib / mean_delta_obs,
                                NA_real_),
      outcome          = outcome_label
    )

  # Cross-stratum decomposition: Δ-difference between strata, RTM-adjusted
  pairs_df <- expand_grid(s_a = per_stratum$stratum, s_b = per_stratum$stratum) %>%
    filter(as.integer(s_a) < as.integer(s_b)) %>%
    rowwise() %>%
    mutate(
      mean_delta_a  = per_stratum$mean_delta_obs[per_stratum$stratum == s_a],
      mean_delta_b  = per_stratum$mean_delta_obs[per_stratum$stratum == s_b],
      diff_obs      = mean_delta_a - mean_delta_b,
      adj_a         = per_stratum$delta_adjusted[per_stratum$stratum == s_a],
      adj_b         = per_stratum$delta_adjusted[per_stratum$stratum == s_b],
      diff_adjusted = adj_a - adj_b,
      diff_rtm      = diff_obs - diff_adjusted,
      pct_rtm       = ifelse(abs(diff_obs) > 0.5, 100 * diff_rtm / diff_obs, NA_real_),
      outcome       = outcome_label
    ) %>%
    ungroup() %>%
    transmute(outcome, s_a, s_b, diff_obs, diff_rtm, diff_adjusted, pct_rtm)

  list(
    outcome     = outcome_label,
    n           = n,
    beta        = beta, beta_se = beta_se,
    alpha       = alpha,
    r2          = r2,
    rho         = rho,
    x_bar       = x_bar,
    per_stratum = per_stratum,
    cross       = pairs_df,
    fit         = fit,
    data        = d
  )
}

sbp_rtm <- fit_rtm("SBP", "BASELINE_SBP", "DELTA_SBP_MO6", "MO6_SBP", cohort)
dbp_rtm <- fit_rtm("DBP", "BASELINE_DBP", "DELTA_DBP_MO6", "MO6_DBP", cohort)

# --- Regression coefficients table ---
reg_table <- bind_rows(
  tibble(outcome = "SBP", n = sbp_rtm$n,
         beta = sbp_rtm$beta, beta_se = sbp_rtm$beta_se,
         alpha = sbp_rtm$alpha, r2 = sbp_rtm$r2,
         rho_test_retest = sbp_rtm$rho,
         cohort_mean_baseline = sbp_rtm$x_bar),
  tibble(outcome = "DBP", n = dbp_rtm$n,
         beta = dbp_rtm$beta, beta_se = dbp_rtm$beta_se,
         alpha = dbp_rtm$alpha, r2 = dbp_rtm$r2,
         rho_test_retest = dbp_rtm$rho,
         cohort_mean_baseline = dbp_rtm$x_bar)
)
write_csv(reg_table, file.path(results_dir, "regression_coefficients.csv"))

cat("--- Regression coefficients (Δ ~ baseline) ---\n")
print(reg_table)

# --- Per-stratum decomposition table ---
per_stratum_table <- bind_rows(sbp_rtm$per_stratum, dbp_rtm$per_stratum) %>%
  select(outcome, stratum, n_users, mean_baseline, x_centered,
         mean_delta_obs, delta_rtm_attrib, delta_adjusted, pct_rtm_attrib, sd_delta)
write_csv(per_stratum_table, file.path(results_dir, "per_stratum_decomposition.csv"))

cat("\n--- Per-stratum RTM decomposition ---\n")
print(per_stratum_table)

# --- Cross-stratum decomposition ---
cross_stratum_table <- bind_rows(sbp_rtm$cross, dbp_rtm$cross)
write_csv(cross_stratum_table, file.path(results_dir, "cross_stratum_decomposition.csv"))

cat("\n--- Cross-stratum RTM decomposition (s_a − s_b) ---\n")
print(cross_stratum_table)

# --- Scatter plot: Δ vs baseline_BP, stratum-colored, with regression line ---
stratum_cols <- c("Elevated" = "#fdae61", "Stage 1" = "#f46d43", "Stage 2" = "#a50026")

plot_df <- bind_rows(
  sbp_rtm$data %>% mutate(outcome = "SBP (mmHg)"),
  dbp_rtm$data %>% mutate(outcome = "DBP (mmHg)")
) %>%
  mutate(outcome = factor(outcome, levels = c("SBP (mmHg)", "DBP (mmHg)")))

# Regression-line data
line_df <- bind_rows(
  tibble(outcome = "SBP (mmHg)",
         baseline = seq(min(sbp_rtm$data$baseline), max(sbp_rtm$data$baseline), length.out = 200)) %>%
    mutate(delta_pred = sbp_rtm$alpha + sbp_rtm$beta * baseline),
  tibble(outcome = "DBP (mmHg)",
         baseline = seq(min(dbp_rtm$data$baseline), max(dbp_rtm$data$baseline), length.out = 200)) %>%
    mutate(delta_pred = dbp_rtm$alpha + dbp_rtm$beta * baseline)
) %>%
  mutate(outcome = factor(outcome, levels = c("SBP (mmHg)", "DBP (mmHg)")))

# Annotation labels for in-panel regression equation
ann_df <- tibble(
  outcome  = factor(c("SBP (mmHg)", "DBP (mmHg)"),
                    levels = c("SBP (mmHg)", "DBP (mmHg)")),
  label    = c(
    sprintf("β = %.3f mmHg/mmHg  (SE %.4f)\nα at X̄ (%.1f) = %.2f mmHg\nR² = %.3f, ρ_test-retest = %.3f, n = %s",
            sbp_rtm$beta, sbp_rtm$beta_se, sbp_rtm$x_bar, sbp_rtm$alpha + sbp_rtm$beta * sbp_rtm$x_bar,
            sbp_rtm$r2, sbp_rtm$rho, format(sbp_rtm$n, big.mark = ",")),
    sprintf("β = %.3f mmHg/mmHg  (SE %.4f)\nα at X̄ (%.1f) = %.2f mmHg\nR² = %.3f, ρ_test-retest = %.3f, n = %s",
            dbp_rtm$beta, dbp_rtm$beta_se, dbp_rtm$x_bar, dbp_rtm$alpha + dbp_rtm$beta * dbp_rtm$x_bar,
            dbp_rtm$r2, dbp_rtm$rho, format(dbp_rtm$n, big.mark = ","))
  )
)

p <- ggplot(plot_df, aes(x = baseline, y = delta)) +
  geom_point(aes(color = stratum), alpha = 0.18, size = 0.5) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.4) +
  geom_line(data = line_df, aes(x = baseline, y = delta_pred),
            color = "black", linewidth = 0.9) +
  geom_text(data = ann_df, aes(label = label, x = -Inf, y = -Inf),
            hjust = -0.05, vjust = -0.3, inherit.aes = FALSE,
            size = 3.2, lineheight = 1.0, color = "gray20") +
  facet_wrap(~ outcome, ncol = 1, scales = "free") +
  scale_color_manual(values = stratum_cols, name = "AHA baseline stratum",
                     guide = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(
    x = "Baseline BP (mmHg)",
    y = "Δ BP (MO6 − Baseline, mmHg)",
    title    = "Oldham regression of Δ on baseline BP — §2.2 cohort",
    subtitle = "Slope β quantifies RTM-attributable share of per-stratum Δ",
    caption  = "RTM-adjusted per-stratum Δs are CONSERVATIVE lower bounds — Oldham approach over-corrects in presence of real dose-response."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 12),
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold"),
    plot.caption     = element_text(size = 9, color = "gray40", hjust = 0)
  )

out_png <- file.path(figures_dir, "5_1_sensitivity_rtm_scatter.png")
ggsave(out_png, p, width = 8, height = 9, dpi = 300)
cat("\nWrote figure:\n  ", out_png, "\n")

# --- Paper-ready markdown ---
fmt_signed <- function(x, d = 2) sprintf("%+.*f", d, x)

md_lines <- c(
  "# §5.1 sensitivity — Oldham/Blomqvist RTM adjustment",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "## Regression of Δ on baseline (Oldham framing)",
  "",
  sprintf("- **SBP** (n = %s): β = %s mmHg/mmHg (SE %.4f); test-retest ρ(baseline, MO6) = %.3f; R² = %.3f. Cohort-mean baseline = %.1f mmHg.",
          format(sbp_rtm$n, big.mark = ","), fmt_signed(sbp_rtm$beta, 3),
          sbp_rtm$beta_se, sbp_rtm$rho, sbp_rtm$r2, sbp_rtm$x_bar),
  sprintf("- **DBP** (n = %s): β = %s mmHg/mmHg (SE %.4f); test-retest ρ(baseline, MO6) = %.3f; R² = %.3f. Cohort-mean baseline = %.1f mmHg.",
          format(dbp_rtm$n, big.mark = ","), fmt_signed(dbp_rtm$beta, 3),
          dbp_rtm$beta_se, dbp_rtm$rho, dbp_rtm$r2, dbp_rtm$x_bar),
  "",
  "## Per-stratum RTM decomposition",
  ""
)
for (otc in c("SBP", "DBP")) {
  md_lines <- c(md_lines, paste0("### ", otc))
  rows <- per_stratum_table %>% filter(outcome == otc)
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    md_lines <- c(md_lines, sprintf(
      "- **%s** (n = %s, mean baseline = %.1f mmHg): observed Δ = %s, RTM-attributable = %s, RTM-adjusted Δ = %s (%.0f%% of observed Δ is RTM-attributable).",
      r$stratum, format(r$n_users, big.mark = ","), r$mean_baseline,
      fmt_signed(r$mean_delta_obs), fmt_signed(r$delta_rtm_attrib),
      fmt_signed(r$delta_adjusted),
      ifelse(is.na(r$pct_rtm_attrib), NA, r$pct_rtm_attrib)
    ))
  }
  md_lines <- c(md_lines, "")
}
md_lines <- c(md_lines,
  "## Notes",
  "",
  "- **Conservative-bound framing**: the Oldham approach attributes ALL baseline-Δ correlation to RTM. If a real dose-response of intervention is present (e.g., wearable engagement works proportionally better at higher baseline), that real signal is wrongly absorbed into the RTM share. The 'RTM-adjusted' Δs reported here are therefore **conservative lower bounds on the true effect** — over-correcting in the presence of any real dose-response. Phrase in the manuscript as 'RTM-conservative bounds,' not 'RTM-corrected estimates.'",
  "- The β coefficient is the slope of Δ on baseline; |β| close to 1 implies strong baseline-Δ correlation (large RTM contribution); |β| close to 0 implies negligible correlation (RTM is small relative to true change).",
  "- ρ_test-retest (correlation between baseline and MO6 BP within the §2.2 cohort) is reported as a sanity check; under classical bivariate-normal assumptions, β = ρ × σ_MO6/σ_baseline − 1.",
  "- The cohort here is §2.2 hypertensive (n = 20,021), matching §5.1 primary. The full-cohort sensitivity (4 strata, n = 28,259) already showed Normal-stratum Δ = +2.7 SBP / +1.3 DBP — direct empirical evidence of RTM. This script provides the *quantitative decomposition* of the §2.2 hypertensive cohort's Δ into RTM-attributable and RTM-adjusted shares.",
  "- The intercept α at the cohort mean baseline X̄ represents the population-mean Δ after removing baseline-dependence — i.e., what the average user would experience if their baseline were at the cohort mean. Read as: 'a typical cohort member's RTM-adjusted Δ.'",
  "- Pairs with: §5.1 primary (unadjusted per-stratum Δs), §5.1 2xWeekly sensitivity (independent RTM bound via measurement-noise reduction), §5.1 full-cohort sensitivity (Normal-stratum direct RTM evidence)."
)
writeLines(md_lines, file.path(results_dir, "5_1_sensitivity_rtm_paper.md"))

cat("\nWrote markdown:\n  ", file.path(results_dir, "5_1_sensitivity_rtm_paper.md"), "\n")
# --- Emit standardized sizes for Table A.2 (paired baseline→MO6 n per outcome) ---
emit_sizes(
  tibble(
    section = "5.1", role = "sensitivity",
    analysis = "RTM decomposition (Oldham/Blomqvist)", cohort = "hypertensive",
    outcome = c("SBP", "DBP"), predictor_set = "RTM-adjusted Δ",
    n_users = c(.sz_sbp, .sz_dbp), n_obs = c(.sz_sbp, .sz_dbp),
    n_events = NA_integer_, ceiling_n = .sz_ceiling
  ),
  "5_1_rtm", results_dir
)

cat("\n=== §5.1 RTM sensitivity done ===\n")
sink()
