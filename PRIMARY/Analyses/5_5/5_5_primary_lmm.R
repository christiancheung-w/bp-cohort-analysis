# 5_5_primary_lmm.R
# SAP §5.5: secondary BP-trajectory LMM with Time extended to 4 levels
# (baseline, 3mo, 6mo, 12mo). Tests whether the §5.1 6-month per-stratum
# decreases are preserved at 12 months.
#
# Identical specification to §5.1 except:
#   - Time factor adds MO12 as a 4th level
#   - Users without MO12 contribute only baseline / 3mo / 6mo rows (REML handles
#     unbalanced data; no users are filtered)
#
# Fixed effects, covariates, AHA reference level, estimand extraction: all
# identical to §5.1.
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/5_5/derived/):
#   fitted_models.rds, emmeans_results.rds
#   results/{anova_tables, r2, emmeans, delta_contrasts, cross_stratum_diffs}.csv
#   results/5_5_primary_lmm.log

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")  # htn_baseline_on(), tidy_contrast(), emit_sizes()

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "5_5", "derived")
results_dir <- file.path(out_dir, "results")
dir.create(out_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "5_5_primary_lmm.log")
sink(log_path, split = TRUE)
cat("§5.5 primary LMM (4-timepoint) —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read user-level analytic frame ---
user_frame <- read_rds(input_path)
cat("User frame:", nrow(user_frame), "rows ×", ncol(user_frame), "cols\n")

# --- §2.2 cohort + sex restriction ---
analysis_users <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female"))
cat("§2.2 cohort + sex restriction:", nrow(analysis_users), "users\n")

# Quick check: how many have an MO12 reading?
cat("Users with MO12_SBP:", sum(!is.na(analysis_users$MO12_SBP)),
    "of", nrow(analysis_users), "(",
    sprintf("%.1f%%", 100 * mean(!is.na(analysis_users$MO12_SBP))), ")\n")
cat("Users with MO12_DBP:", sum(!is.na(analysis_users$MO12_DBP)), "\n\n")

# --- Pivot to long: now 4 timepoints (BASELINE, MO3, MO6, MO12) ---
long_outcomes <- analysis_users %>%
  select(STUDY_ID, BASELINE_AHA, AGE_AT_START, SEX, BASELINE_BMI, HTN_MED_BASELINE,
         BASELINE_SBP,           MO3_SBP,           MO6_SBP,         MO12_SBP,
         BASELINE_DBP,           MO3_DBP,           MO6_DBP,         MO12_DBP,
         BASELINE_TIME_OF_DAY,   MO3_TIME_OF_DAY,   MO6_TIME_OF_DAY, MO12_TIME_OF_DAY,
         BASELINE_SEASON,        MO3_SEASON,        MO6_SEASON,      MO12_SEASON) %>%
  pivot_longer(
    cols = c(starts_with("BASELINE_S"), starts_with("MO3_S"),  starts_with("MO6_S"),  starts_with("MO12_S"),
             starts_with("BASELINE_D"), starts_with("MO3_D"),  starts_with("MO6_D"),  starts_with("MO12_D"),
             starts_with("BASELINE_T"), starts_with("MO3_T"),  starts_with("MO6_T"),  starts_with("MO12_T")),
    names_to = c("Time", ".value"),
    names_pattern = "(BASELINE|MO3|MO6|MO12)_(.+)"
  )

long_outcomes <- long_outcomes %>%
  mutate(
    Time         = factor(Time, levels = c("BASELINE", "MO3", "MO6", "MO12")),
    BASELINE_AHA = factor(BASELINE_AHA,
                          levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")),
    SEX          = factor(SEX),
    SEASON       = factor(SEASON, levels = c("winter", "spring", "summer", "fall")),
    HTN_BASELINE_ON = htn_baseline_on(HTN_MED_BASELINE)  # baseline-only HTN trait (v0.79)
  )

cat("\nLong-format frame:", nrow(long_outcomes), "rows total\n")
cat("Rows per (Time × stratum) cell (SBP non-NA only):\n")
print(long_outcomes %>%
        filter(!is.na(SBP)) %>%
        count(Time, BASELINE_AHA) %>%
        pivot_wider(names_from = BASELINE_AHA, values_from = n))

# --- Helper: fit + summarize one LMM ---
fit_lmm <- function(outcome_col, adjusted = FALSE, data) {
  if (adjusted) {
    f <- as.formula(paste0(
      outcome_col, " ~ Time * BASELINE_AHA",
      " + AGE_AT_START + SEX + BASELINE_BMI + HTN_BASELINE_ON",
      " + sin(2 * pi * TIME_OF_DAY / 24) + cos(2 * pi * TIME_OF_DAY / 24)",
      " + SEASON",
      " + Time:AGE_AT_START + Time:SEX + Time:BASELINE_BMI + Time:HTN_BASELINE_ON",
      " + (1 | STUDY_ID)"
    ))
  } else {
    f <- as.formula(paste0(
      outcome_col, " ~ Time * BASELINE_AHA + (1 | STUDY_ID)"
    ))
  }
  d <- data %>% filter(!is.na(.data[[outcome_col]]), !is.na(BASELINE_AHA))
  if (adjusted) {
    d <- d %>% filter(!is.na(AGE_AT_START), !is.na(SEX),
                      !is.na(BASELINE_BMI), !is.na(HTN_BASELINE_ON),
                      !is.na(TIME_OF_DAY), !is.na(SEASON))
  }
  cat(sprintf("  Fitting on %d rows, %d users…\n",
              nrow(d), n_distinct(d$STUDY_ID)))
  lmer(f, data = d, REML = TRUE)
}

print_model_summary <- function(fit, label) {
  cat("\n\n=== ", label, " ===\n", sep = "")
  cat("Formula:", deparse(formula(fit)), "\n")
  cat("n_obs:", nobs(fit), "; n_users:", lme4::ngrps(fit), "\n\n")
  cat("--- Type III ANOVA (Satterthwaite df) ---\n")
  print(suppressMessages(anova(fit, type = 3)))
  cat("\n--- Nakagawa R² ---\n")
  print(performance::r2_nakagawa(fit))
}

# --- Fit 4 models ---
cat("\n\n========================================\n")
cat("Fitting Unadjusted LMM (SBP) — 4-timepoint\n")
cat("========================================\n")
fit_sbp_unadj <- fit_lmm("SBP", adjusted = FALSE, data = long_outcomes)
print_model_summary(fit_sbp_unadj, "SBP — UNADJUSTED (4-timepoint)")

cat("\n\n========================================\n")
cat("Fitting Adjusted LMM (SBP) — 4-timepoint\n")
cat("========================================\n")
fit_sbp_adj <- fit_lmm("SBP", adjusted = TRUE, data = long_outcomes)
print_model_summary(fit_sbp_adj, "SBP — ADJUSTED (4-timepoint)")

cat("\n\n========================================\n")
cat("Fitting Unadjusted LMM (DBP) — 4-timepoint\n")
cat("========================================\n")
fit_dbp_unadj <- fit_lmm("DBP", adjusted = FALSE, data = long_outcomes)
print_model_summary(fit_dbp_unadj, "DBP — UNADJUSTED (4-timepoint)")

cat("\n\n========================================\n")
cat("Fitting Adjusted LMM (DBP) — 4-timepoint\n")
cat("========================================\n")
fit_dbp_adj <- fit_lmm("DBP", adjusted = TRUE, data = long_outcomes)
print_model_summary(fit_dbp_adj, "DBP — ADJUSTED (4-timepoint)")

# --- Estimand extraction: emmeans, within-stratum Δ, cross-stratum Δ-diff ---
extract_emmeans <- function(fit, label) {
  cat("\n\n=== Estimand extraction:", label, " ===\n")
  emm <- suppressMessages(emmeans(fit, ~ Time | BASELINE_AHA))
  cat("\n--- Estimated marginal means (Time × Stratum) ---\n")
  print(emm)
  cat("\n--- Within-stratum Δ from BASELINE (3mo, 6mo, 12mo) ---\n")
  delta_emm <- suppressMessages(contrast(emm, method = "trt.vs.ctrl", adjust = "none"))
  print(delta_emm)
  cat("\n--- Cross-stratum Δ-difference at each follow-up timepoint ---\n")
  diff_emm <- suppressMessages(contrast(delta_emm, "pairwise", by = "contrast", adjust = "none"))
  print(diff_emm)
  cat("\n--- Within-stratum CONSECUTIVE Δ (adjacent timepoints: MO3-BASELINE, MO6-MO3, MO12-MO6) ---\n")
  consec_emm <- suppressMessages(contrast(emm, method = "consec", adjust = "none"))
  print(consec_emm)
  list(
    emmeans        = as_tibble(as.data.frame(emm)),
    within_stratum = tidy_contrast(delta_emm),    # 95% CI + within-family Holm (v0.79)
    cross_stratum  = tidy_contrast(diff_emm),
    consecutive    = tidy_contrast(consec_emm)
  )
}

emm_sbp_unadj <- extract_emmeans(fit_sbp_unadj, "SBP — UNADJUSTED")
emm_sbp_adj   <- extract_emmeans(fit_sbp_adj,   "SBP — ADJUSTED")
emm_dbp_unadj <- extract_emmeans(fit_dbp_unadj, "DBP — UNADJUSTED")
emm_dbp_adj   <- extract_emmeans(fit_dbp_adj,   "DBP — ADJUSTED")

# --- Save rds ---
fitted_models <- list(
  sbp_unadj = fit_sbp_unadj, sbp_adj = fit_sbp_adj,
  dbp_unadj = fit_dbp_unadj, dbp_adj = fit_dbp_adj
)
emmeans_results <- list(
  sbp_unadj = emm_sbp_unadj, sbp_adj = emm_sbp_adj,
  dbp_unadj = emm_dbp_unadj, dbp_adj = emm_dbp_adj
)
write_rds(fitted_models,   file.path(out_dir, "fitted_models.rds"))
write_rds(emmeans_results, file.path(out_dir, "emmeans_results.rds"))

# --- Emit standardized model sizes for Table A (from the fits) ---
# §5.5 is a secondary (not sensitivity) analysis → role = "primary".
emit_sizes(
  imap_dfr(fitted_models, function(fit, key) tibble(
    section = "5.5", role = "primary",
    analysis = "12-month BP-trajectory LMM (4-timepoint)", cohort = "hypertensive",
    outcome = if (startsWith(key, "sbp")) "SBP" else "DBP",
    predictor_set = if (grepl("unadj", key)) "unadjusted" else "adjusted",
    n_users = as.integer(lme4::ngrps(fit)), n_obs = as.integer(nobs(fit)),
    n_events = NA_integer_, ceiling_n = nrow(analysis_users)
  )),
  "5_5", results_dir
)

# --- Structured CSV exports ---
anova_table <- imap_dfr(fitted_models, function(fit, name) {
  a <- suppressMessages(anova(fit, type = 3))
  as_tibble(a, rownames = "term") %>% mutate(model = name) %>% select(model, everything())
})
write_csv(anova_table, file.path(results_dir, "5_5_anova_tables.csv"))

r2_table <- imap_dfr(fitted_models, function(fit, name) {
  r2 <- suppressMessages(r2_nakagawa(fit))
  tibble(model = name,
         marginal_r2    = as.numeric(r2$R2_marginal),
         conditional_r2 = as.numeric(r2$R2_conditional))
})
write_csv(r2_table, file.path(results_dir, "5_5_r2.csv"))

emmeans_table <- imap_dfr(emmeans_results, function(r, name) {
  r$emmeans %>% mutate(model = name) %>% select(model, everything())
})
write_csv(emmeans_table, file.path(results_dir, "5_5_emmeans.csv"))

delta_table <- imap_dfr(emmeans_results, function(r, name) {
  r$within_stratum %>% mutate(model = name) %>% select(model, everything())
})
write_csv(delta_table, file.path(results_dir, "5_5_delta_contrasts.csv"))

cross_table <- imap_dfr(emmeans_results, function(r, name) {
  r$cross_stratum %>% mutate(model = name) %>% select(model, everything())
})
write_csv(cross_table, file.path(results_dir, "5_5_cross_stratum_diffs.csv"))

# Within-stratum CONSECUTIVE Δ contrasts (adjacent timepoints: MO3-BASELINE, MO6-MO3, MO12-MO6)
# Added for figure significance indicators + later verification of plateau/continued-change.
consec_table <- imap_dfr(emmeans_results, function(r, name) {
  r$consecutive %>% mutate(model = name) %>% select(model, everything())
})
write_csv(consec_table, file.path(results_dir, "5_5_consec_contrasts.csv"))

cat("\n\n=== §5.5 primary LMM done ===\n")
cat("Saved rds:\n  ", file.path(out_dir, "fitted_models.rds"),
    "\n  ", file.path(out_dir, "emmeans_results.rds"), "\n")
cat("Saved CSVs to:\n  ", results_dir, "\n")
for (f in list.files(results_dir, pattern = "\\.csv$")) cat("    ", f, "\n")
cat("Saved log:\n  ", log_path, "\n")
cat("\nNext: 5_5_diagnostics.R, then 5_5_figure_2.R, then companions for §5.2/§5.4 at 12mo.\n")

sink()
