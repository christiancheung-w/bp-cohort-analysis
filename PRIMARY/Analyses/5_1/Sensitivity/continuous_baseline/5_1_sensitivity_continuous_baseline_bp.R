# 5_1_sensitivity_continuous_baseline_bp.R
# SAP §5.1 sensitivity: replace AHA categorical stratum with continuous
# baseline BP + Time:baseline_BP interaction. Tests whether the §5.1 stratum-
# level findings are robust to the AHA cutpoint choice — does the same
# gradient emerge from continuous baseline, or does the categorical model
# capture something continuous misses?
#
# Implementation: fit on follow-up timepoints only (MO3, MO6) with
# baseline_BP_centered (at cohort mean) as a static user-level covariate.
# BASELINE rows are dropped because Y = baseline_BP at baseline is tautological
# and creates degeneracy when baseline_BP is included as a covariate.
#
# Model structure:
#   Unadjusted: Y_t ~ Time + baseline_BP_centered + Time:baseline_BP_centered + (1|STUDY_ID)
#   Adjusted:   + AGE_AT_START + SEX + BASELINE_BMI + cyclic time_of_day + SEASON
#               + Time:AGE_AT_START + Time:SEX + Time:BASELINE_BMI
#
# Cohort: §2.2 hypertensive (Elevated + Stage 1 + Stage 2), matching §5.1.
# Outputs include predicted Δ from baseline at three baseline_BP values
# corresponding to AHA stratum means, for direct comparison to §5.1.
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/5_1/Sensitivity/continuous_baseline/derived/):
#   fitted_models.rds
#   results/regression_coefficients.csv
#   results/predicted_delta_at_stratum_means.csv  — for direct comparison to §5.1
#   results/5_1_sensitivity_continuous_baseline_paper.md
#   results/5_1_sensitivity_continuous_baseline.log
#   figures/5_1_continuous_baseline_curve.png      — predicted Δ vs baseline_BP

library(tidyverse)
library(lme4)
library(lmerTest)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")  # htn_baseline_on(), emit_sizes()

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "5_1", "Sensitivity",
                         "continuous_baseline", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(out_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "5_1_sensitivity_continuous_baseline.log")
sink(log_path, split = TRUE)
cat("§5.1 sensitivity — continuous baseline BP —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

user_frame <- read_rds(input_path)

# §2.2 cohort + sex restriction
analysis_users <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female"))
cat("§2.2 cohort + sex restriction:", nrow(analysis_users), "users\n")

# Long format: follow-up timepoints only
long_outcomes <- analysis_users %>%
  select(STUDY_ID, BASELINE_AHA, AGE_AT_START, SEX, BASELINE_BMI, HTN_MED_BASELINE,
         BASELINE_SBP, MO3_SBP, MO6_SBP,
         BASELINE_DBP, MO3_DBP, MO6_DBP,
         MO3_TIME_OF_DAY, MO6_TIME_OF_DAY,
         MO3_SEASON,      MO6_SEASON) %>%
  pivot_longer(
    cols = c(starts_with("MO3_"), starts_with("MO6_")),
    names_to = c("Time", ".value"),
    names_pattern = "(MO3|MO6)_(.+)"
  ) %>%
  mutate(
    Time         = factor(Time, levels = c("MO3", "MO6")),
    SEX          = factor(SEX),
    SEASON       = factor(SEASON, levels = c("winter", "spring", "summer", "fall")),
    BASELINE_AHA = factor(BASELINE_AHA,
                          levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")),
    HTN_BASELINE_ON = htn_baseline_on(HTN_MED_BASELINE)  # baseline-only HTN trait (v0.79)
  )

# Cohort-mean baseline for centering (per outcome)
x_bar_sbp <- mean(analysis_users$BASELINE_SBP, na.rm = TRUE)
x_bar_dbp <- mean(analysis_users$BASELINE_DBP, na.rm = TRUE)
cat("Cohort-mean baseline: SBP =", round(x_bar_sbp, 2),
    " | DBP =", round(x_bar_dbp, 2), "\n")

long_outcomes <- long_outcomes %>%
  mutate(
    BASELINE_SBP_c = BASELINE_SBP - x_bar_sbp,
    BASELINE_DBP_c = BASELINE_DBP - x_bar_dbp
  )

# Stratum-mean baselines (for predicted-Δ comparison)
stratum_means <- analysis_users %>%
  group_by(BASELINE_AHA) %>%
  summarise(
    n            = n(),
    mean_BL_SBP  = mean(BASELINE_SBP),
    mean_BL_DBP  = mean(BASELINE_DBP),
    .groups = "drop"
  ) %>%
  mutate(
    stratum     = factor(BASELINE_AHA,
                         levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
                         labels = c("Elevated", "Stage 1", "Stage 2")),
    SBP_c       = mean_BL_SBP - x_bar_sbp,
    DBP_c       = mean_BL_DBP - x_bar_dbp
  )
cat("\nStratum-mean baselines (for prediction comparison):\n")
print(stratum_means)

# --- Helper: fit LMM ---
fit_lmm_cont <- function(outcome_col, baseline_centered_col, adjusted, data) {
  if (adjusted) {
    f <- as.formula(paste0(
      outcome_col, " ~ Time + ", baseline_centered_col,
      " + Time:", baseline_centered_col,
      " + AGE_AT_START + SEX + BASELINE_BMI + HTN_BASELINE_ON",
      " + sin(2 * pi * TIME_OF_DAY / 24) + cos(2 * pi * TIME_OF_DAY / 24)",
      " + SEASON",
      " + Time:AGE_AT_START + Time:SEX + Time:BASELINE_BMI + Time:HTN_BASELINE_ON",
      " + (1 | STUDY_ID)"
    ))
  } else {
    f <- as.formula(paste0(
      outcome_col, " ~ Time + ", baseline_centered_col,
      " + Time:", baseline_centered_col,
      " + (1 | STUDY_ID)"
    ))
  }
  d <- data %>% filter(!is.na(.data[[outcome_col]]),
                       !is.na(.data[[baseline_centered_col]]))
  if (adjusted) {
    d <- d %>% filter(!is.na(AGE_AT_START), !is.na(SEX),
                      !is.na(BASELINE_BMI), !is.na(HTN_BASELINE_ON),
                      !is.na(TIME_OF_DAY), !is.na(SEASON))
  }
  cat(sprintf("  Fitting on %d rows, %d users…\n",
              nrow(d), n_distinct(d$STUDY_ID)))
  lmer(f, data = d, REML = TRUE)
}

# --- Fit 4 models ---
cat("\n\n========================================\nSBP UNADJUSTED — continuous baseline\n========================================\n")
fit_sbp_unadj <- fit_lmm_cont("SBP", "BASELINE_SBP_c", FALSE, long_outcomes)
cat("\n--- Type III ANOVA ---\n")
print(suppressMessages(anova(fit_sbp_unadj, type = 3)))

cat("\n\n========================================\nSBP ADJUSTED — continuous baseline\n========================================\n")
fit_sbp_adj <- fit_lmm_cont("SBP", "BASELINE_SBP_c", TRUE, long_outcomes)
cat("\n--- Type III ANOVA ---\n")
print(suppressMessages(anova(fit_sbp_adj, type = 3)))

cat("\n\n========================================\nDBP UNADJUSTED — continuous baseline\n========================================\n")
fit_dbp_unadj <- fit_lmm_cont("DBP", "BASELINE_DBP_c", FALSE, long_outcomes)
cat("\n--- Type III ANOVA ---\n")
print(suppressMessages(anova(fit_dbp_unadj, type = 3)))

cat("\n\n========================================\nDBP ADJUSTED — continuous baseline\n========================================\n")
fit_dbp_adj <- fit_lmm_cont("DBP", "BASELINE_DBP_c", TRUE, long_outcomes)
cat("\n--- Type III ANOVA ---\n")
print(suppressMessages(anova(fit_dbp_adj, type = 3)))

# --- Regression coefficients ---
extract_coefs <- function(fit, outcome, adjusted_label) {
  s <- summary(fit)
  cf <- s$coefficients
  as_tibble(cf, rownames = "term") %>%
    rename(estimate = Estimate, SE = `Std. Error`, df = df, t_value = `t value`, p_value = `Pr(>|t|)`) %>%
    mutate(outcome = outcome, model = adjusted_label) %>%
    select(outcome, model, term, estimate, SE, df, t_value, p_value)
}

reg_table <- bind_rows(
  extract_coefs(fit_sbp_unadj, "SBP", "unadj"),
  extract_coefs(fit_sbp_adj,   "SBP", "adj"),
  extract_coefs(fit_dbp_unadj, "DBP", "unadj"),
  extract_coefs(fit_dbp_adj,   "DBP", "adj")
)
write_csv(reg_table, file.path(results_dir, "regression_coefficients.csv"))

cat("\n\n--- Key coefficients (Time + baseline_BP_centered + Time:baseline_BP_centered) ---\n")
print(reg_table %>%
        filter(term %in% c("(Intercept)", "TimeMO6",
                           "BASELINE_SBP_c", "BASELINE_DBP_c",
                           "TimeMO6:BASELINE_SBP_c", "TimeMO6:BASELINE_DBP_c")) %>%
        select(outcome, model, term, estimate, SE, p_value))

# --- Predicted Δ at AHA stratum mean baselines ---
predict_delta <- function(fit, outcome_label, model_label,
                          baseline_centered_col, stratum_centered_col) {
  # newdata: one row per (Time × stratum), with covariates at reference levels
  d <- expand_grid(
    Time = factor(c("MO3", "MO6"), levels = c("MO3", "MO6")),
    stratum = stratum_means$stratum
  ) %>%
    left_join(stratum_means, by = "stratum") %>%
    mutate(
      AGE_AT_START = mean(analysis_users$AGE_AT_START, na.rm = TRUE),
      SEX          = factor("male", levels = levels(long_outcomes$SEX)),
      BASELINE_BMI = mean(analysis_users$BASELINE_BMI, na.rm = TRUE),
      TIME_OF_DAY  = mean(long_outcomes$TIME_OF_DAY, na.rm = TRUE),
      SEASON       = factor("winter", levels = levels(long_outcomes$SEASON)),
      HTN_BASELINE_ON = factor("not_on_at_baseline", levels = levels(long_outcomes$HTN_BASELINE_ON))
    )
  d[[baseline_centered_col]] <- d[[stratum_centered_col]]

  # Predicted Y at the follow-up timepoint
  d$predicted_Y <- predict(fit, newdata = d, re.form = NA)

  # Baseline value at this stratum
  d$baseline_value <- if (baseline_centered_col == "BASELINE_SBP_c") d$mean_BL_SBP else d$mean_BL_DBP

  # Predicted Δ from baseline
  d %>%
    mutate(predicted_delta = predicted_Y - baseline_value,
           outcome         = outcome_label,
           model           = model_label) %>%
    select(outcome, model, Time, stratum, baseline_value, predicted_Y, predicted_delta)
}

pred_table <- bind_rows(
  predict_delta(fit_sbp_unadj, "SBP", "unadj", "BASELINE_SBP_c", "SBP_c"),
  predict_delta(fit_sbp_adj,   "SBP", "adj",   "BASELINE_SBP_c", "SBP_c"),
  predict_delta(fit_dbp_unadj, "DBP", "unadj", "BASELINE_DBP_c", "DBP_c"),
  predict_delta(fit_dbp_adj,   "DBP", "adj",   "BASELINE_DBP_c", "DBP_c")
) %>%
  arrange(outcome, model, stratum, Time)

write_csv(pred_table, file.path(results_dir, "predicted_delta_at_stratum_means.csv"))

cat("\n\n--- Predicted Δ at AHA stratum mean baselines (continuous model) ---\n")
print(pred_table, n = Inf, width = Inf)

# --- Save fits ---
fitted_models <- list(
  sbp_unadj = fit_sbp_unadj, sbp_adj = fit_sbp_adj,
  dbp_unadj = fit_dbp_unadj, dbp_adj = fit_dbp_adj
)
write_rds(fitted_models, file.path(out_dir, "fitted_models.rds"))

# --- Emit standardized model sizes for Table A.2 (from the fits) ---
# NB: no emmeans contrast families here (raw regression coefficients), so the
# universal-Holm convention does not apply to this sensitivity.
emit_sizes(
  imap_dfr(fitted_models, function(fit, key) tibble(
    section = "5.1", role = "sensitivity",
    analysis = "Continuous-baseline-BP LMM (cutpoint-free)", cohort = "hypertensive",
    outcome = if (startsWith(key, "sbp")) "SBP" else "DBP",
    predictor_set = if (grepl("unadj", key)) "unadjusted" else "adjusted",
    n_users = as.integer(lme4::ngrps(fit)), n_obs = as.integer(nobs(fit)),
    n_events = NA_integer_, ceiling_n = nrow(analysis_users)
  )),
  "5_1_continuous_baseline", results_dir
)

# --- Figure: predicted Δ at MO6 (adjusted) as a function of baseline_BP ---
make_curve <- function(fit, outcome_label, baseline_centered_col, x_bar) {
  grid <- tibble(
    baseline_value = if (outcome_label == "SBP") seq(115, 160, by = 1) else seq(70, 105, by = 1),
    Time           = factor("MO6", levels = c("MO3", "MO6")),
    AGE_AT_START   = mean(analysis_users$AGE_AT_START, na.rm = TRUE),
    SEX            = factor("male", levels = levels(long_outcomes$SEX)),
    BASELINE_BMI   = mean(analysis_users$BASELINE_BMI, na.rm = TRUE),
    TIME_OF_DAY    = mean(long_outcomes$TIME_OF_DAY, na.rm = TRUE),
    SEASON         = factor("winter", levels = levels(long_outcomes$SEASON)),
    HTN_BASELINE_ON = factor("not_on_at_baseline", levels = levels(long_outcomes$HTN_BASELINE_ON))
  )
  grid[[baseline_centered_col]] <- grid$baseline_value - x_bar
  grid$predicted_Y <- predict(fit, newdata = grid, re.form = NA)
  grid %>%
    mutate(predicted_delta = predicted_Y - baseline_value,
           outcome         = outcome_label) %>%
    select(outcome, baseline_value, predicted_delta)
}

curve_df <- bind_rows(
  make_curve(fit_sbp_adj, "SBP", "BASELINE_SBP_c", x_bar_sbp),
  make_curve(fit_dbp_adj, "DBP", "BASELINE_DBP_c", x_bar_dbp)
) %>%
  mutate(outcome = factor(outcome, levels = c("SBP", "DBP")))

# Stratum marker points (using the adjusted predictions at MO6)
stratum_markers <- pred_table %>%
  filter(Time == "MO6", model == "adj") %>%
  mutate(outcome = factor(outcome, levels = c("SBP", "DBP")))

stratum_cols <- c("Elevated" = "#fdae61", "Stage 1" = "#f46d43", "Stage 2" = "#a50026")

p <- ggplot(curve_df, aes(x = baseline_value, y = predicted_delta)) +
  geom_hline(yintercept = 0, color = "gray60", linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 1, color = "gray20") +
  geom_point(data = stratum_markers, aes(color = stratum), size = 3.5) +
  geom_text(data = stratum_markers,
            aes(label = sprintf("%s\nbaseline=%.0f\nΔ=%.2f",
                                stratum, baseline_value, predicted_delta),
                color = stratum),
            hjust = -0.15, vjust = 0.5, size = 3, lineheight = 0.9) +
  facet_wrap(~ outcome, ncol = 1, scales = "free") +
  scale_color_manual(values = stratum_cols, name = NULL, guide = "none") +
  labs(
    x = "Baseline BP (mmHg)",
    y = "Predicted Δ at MO6 (mmHg)",
    title    = "Continuous baseline BP model — predicted MO6 Δ as function of baseline",
    subtitle = "Adjusted LMM on §2.2 cohort; AHA stratum-mean baselines marked",
    caption  = "Curve shows continuous-baseline model prediction. Marker points are predicted Δ at each AHA stratum's mean baseline — directly comparable to §5.1 per-stratum estimates."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 12),
    plot.title       = element_text(face = "bold"),
    plot.caption     = element_text(size = 9, color = "gray40", hjust = 0)
  )

out_png <- file.path(figures_dir, "5_1_continuous_baseline_curve.png")
ggsave(out_png, p, width = 9, height = 9, dpi = 300)
cat("\nWrote figure:\n  ", out_png, "\n")

# --- Paper-ready markdown ---
fmt_signed <- function(x, d = 2) sprintf("%+.*f", d, x)

# Extract key Time:baseline interaction p-values to test cutpoint sensitivity
key_terms <- reg_table %>%
  filter(term %in% c("TimeMO6:BASELINE_SBP_c", "TimeMO6:BASELINE_DBP_c"))

md_lines <- c(
  "# §5.1 sensitivity — continuous baseline BP (alt LMM structure)",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "## Predicted Δ at AHA stratum mean baselines",
  "",
  "Continuous-baseline LMM, MO6 timepoint, **adjusted** model. Compare to §5.1 per-stratum estimates.",
  ""
)
mo6_adj <- pred_table %>% filter(Time == "MO6", model == "adj")
for (otc in c("SBP", "DBP")) {
  md_lines <- c(md_lines, paste0("### ", otc))
  rows <- mo6_adj %>% filter(outcome == otc)
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    md_lines <- c(md_lines, sprintf(
      "- **%s** (baseline = %.1f mmHg): predicted MO6 Δ = %s mmHg.",
      r$stratum, r$baseline_value, fmt_signed(r$predicted_delta)
    ))
  }
  md_lines <- c(md_lines, "")
}

md_lines <- c(md_lines,
  "## Key regression coefficients (Time:baseline_BP interaction)",
  ""
)
for (i in seq_len(nrow(key_terms))) {
  r <- key_terms[i, ]
  md_lines <- c(md_lines, sprintf(
    "- **%s**, term `%s`: estimate = %s mmHg per mmHg (SE %.4f, p = %.2e).",
    paste0(r$outcome, " ", r$model),
    r$term, fmt_signed(r$estimate, 4), r$SE, r$p_value
  ))
}

md_lines <- c(md_lines,
  "",
  "## Notes",
  "",
  "- Continuous baseline model fitted on follow-up timepoints (MO3, MO6) with baseline_BP centered at cohort mean; BASELINE rows excluded as tautological (Y = baseline_BP at baseline creates degeneracy when baseline_BP is included as a covariate).",
  "- The `Time:baseline_BP` interaction term tests whether the slope of follow-up BP on baseline BP differs between MO3 and MO6. Significant interaction indicates the regression-to-mean dynamic operates differently at the two timepoints.",
  "- Predicted Δ values at AHA stratum mean baselines are directly comparable to §5.1 per-stratum estimates. If the continuous and categorical formulations agree closely, the AHA cutpoints add nothing beyond what continuous baseline already captures. If they disagree substantially, the cutpoints capture some nonlinearity that continuous misses.",
  "- A reviewer asking 'are your results driven by the AHA cutpoint choice' is answered by this sensitivity: the continuous model is cutpoint-free.",
  "- The curve shape (Δ as function of baseline_BP) directly visualizes the dose-response — and the per-stratum dose-response is just three sample points along this curve.",
  "- Pairs with: §5.1 primary (categorical AHA stratum), §5.1 RTM sensitivity (Oldham/Blomqvist quantifies share of curve slope attributable to RTM)."
)
writeLines(md_lines, file.path(results_dir, "5_1_sensitivity_continuous_baseline_paper.md"))

cat("\nWrote markdown:\n  ", file.path(results_dir, "5_1_sensitivity_continuous_baseline_paper.md"), "\n")
cat("\n=== §5.1 continuous-baseline sensitivity done ===\n")
sink()
