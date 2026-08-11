# 5_1_sensitivity_htn_naive.R
# SAP §5.1 sensitivity bullet 7: refit the primary LMM on the
# HTN-medication-naive subset (HTN_MED_BASELINE == "never_logged").
#
# Tests whether the §5.1 BP-reduction findings are driven by users on or
# initiating HTN medications during the study window. If the finding survives
# restriction to users with no HTN-med entries in the 180-day lookback OR
# the study window, it cannot be attributed to pharmacological treatment.
#
# Cohort: §2.2 hypertensive (Elevated + Stage 1 + Stage 2) AND
# HTN_MED_BASELINE == "never_logged" — ~16,336 users (81.6% of §2.2).
#
# Otherwise structurally identical to §5.1 primary LMM.
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/5_1/Sensitivity/htn_naive/derived/):
#   fitted_models.rds, emmeans_results.rds
#   results/{anova_tables, r2, emmeans, delta_contrasts, cross_stratum_diffs}.csv
#   results/5_1_sensitivity_htn_naive.log

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")  # tidy_contrast(), emit_sizes() — NB: no baseline-HTN covariate here (constant: all never_logged)

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "5_1", "Sensitivity", "htn_naive", "derived")
results_dir <- file.path(out_dir, "results")
dir.create(out_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "5_1_sensitivity_htn_naive.log")
sink(log_path, split = TRUE)
cat("§5.1 sensitivity — HTN-medication-naive subset —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

user_frame <- read_rds(input_path)
cat("User frame:", nrow(user_frame), "rows ×", ncol(user_frame), "cols\n")

# --- §2.2 cohort + HTN-naive filter + sex restriction ---
analysis_users <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         HTN_MED_BASELINE == "never_logged",
         SEX %in% c("male", "female"))
cat("§2.2 cohort:", sum(user_frame$BASELINE_AHA %in%
                          c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")), "\n")
cat("After HTN-naive filter:", sum(user_frame$BASELINE_AHA %in%
                          c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension") &
                          user_frame$HTN_MED_BASELINE == "never_logged"), "\n")
cat("After sex restriction:", nrow(analysis_users), "\n")
cat("Per-stratum HTN-naive counts:\n")
print(analysis_users %>% count(BASELINE_AHA))

# --- Long format ---
long_outcomes <- analysis_users %>%
  select(STUDY_ID, BASELINE_AHA, AGE_AT_START, SEX, BASELINE_BMI,
         BASELINE_SBP,           MO3_SBP,           MO6_SBP,
         BASELINE_DBP,           MO3_DBP,           MO6_DBP,
         BASELINE_TIME_OF_DAY,   MO3_TIME_OF_DAY,   MO6_TIME_OF_DAY,
         BASELINE_SEASON,        MO3_SEASON,        MO6_SEASON) %>%
  pivot_longer(
    cols = c(starts_with("BASELINE_S"), starts_with("MO3_S"), starts_with("MO6_S"),
             starts_with("BASELINE_D"), starts_with("MO3_D"), starts_with("MO6_D"),
             starts_with("BASELINE_T"), starts_with("MO3_T"), starts_with("MO6_T")),
    names_to = c("Time", ".value"),
    names_pattern = "(BASELINE|MO3|MO6)_(.+)"
  ) %>%
  mutate(
    Time         = factor(Time, levels = c("BASELINE", "MO3", "MO6")),
    BASELINE_AHA = factor(BASELINE_AHA,
                          levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")),
    SEX          = factor(SEX),
    SEASON       = factor(SEASON, levels = c("winter", "spring", "summer", "fall"))
  )

cat("\nLong-format frame:", nrow(long_outcomes), "rows total\n")
cat("Rows per (Time × stratum) cell (SBP non-NA):\n")
print(long_outcomes %>%
        filter(!is.na(SBP)) %>%
        count(Time, BASELINE_AHA) %>%
        pivot_wider(names_from = BASELINE_AHA, values_from = n))

# --- Helpers (identical to primary §5.1) ---
fit_lmm <- function(outcome_col, adjusted = FALSE, data) {
  if (adjusted) {
    f <- as.formula(paste0(
      outcome_col, " ~ Time * BASELINE_AHA",
      " + AGE_AT_START + SEX + BASELINE_BMI",
      " + sin(2 * pi * TIME_OF_DAY / 24) + cos(2 * pi * TIME_OF_DAY / 24)",
      " + SEASON",
      " + Time:AGE_AT_START + Time:SEX + Time:BASELINE_BMI",
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
                      !is.na(BASELINE_BMI), !is.na(TIME_OF_DAY), !is.na(SEASON))
  }
  cat(sprintf("  Fitting on %d rows, %d users…\n",
              nrow(d), n_distinct(d$STUDY_ID)))
  lmer(f, data = d, REML = TRUE)
}

print_model_summary <- function(fit, label) {
  cat("\n\n=== ", label, " ===\n", sep = "")
  cat("Formula:", deparse(formula(fit)), "\n")
  cat("n_obs:", nobs(fit), "; n_users:", lme4::ngrps(fit), "\n\n")
  cat("--- Type III ANOVA ---\n")
  print(suppressMessages(anova(fit, type = 3)))
  cat("\n--- Nakagawa R² ---\n")
  print(performance::r2_nakagawa(fit))
}

# --- Fit 4 models ---
cat("\n\n========================================\nSBP UNADJUSTED — HTN-naive\n========================================\n")
fit_sbp_unadj <- fit_lmm("SBP", adjusted = FALSE, data = long_outcomes)
print_model_summary(fit_sbp_unadj, "SBP — UNADJUSTED (HTN-naive)")

cat("\n\n========================================\nSBP ADJUSTED — HTN-naive\n========================================\n")
fit_sbp_adj <- fit_lmm("SBP", adjusted = TRUE, data = long_outcomes)
print_model_summary(fit_sbp_adj, "SBP — ADJUSTED (HTN-naive)")

cat("\n\n========================================\nDBP UNADJUSTED — HTN-naive\n========================================\n")
fit_dbp_unadj <- fit_lmm("DBP", adjusted = FALSE, data = long_outcomes)
print_model_summary(fit_dbp_unadj, "DBP — UNADJUSTED (HTN-naive)")

cat("\n\n========================================\nDBP ADJUSTED — HTN-naive\n========================================\n")
fit_dbp_adj <- fit_lmm("DBP", adjusted = TRUE, data = long_outcomes)
print_model_summary(fit_dbp_adj, "DBP — ADJUSTED (HTN-naive)")

# --- Estimand extraction ---
extract_emmeans <- function(fit, label) {
  cat("\n\n=== Estimand extraction:", label, " ===\n")
  emm <- suppressMessages(emmeans(fit, ~ Time | BASELINE_AHA))
  cat("\n--- emmeans (Time × Stratum) ---\n"); print(emm)
  delta_emm <- suppressMessages(contrast(emm, method = "trt.vs.ctrl", adjust = "none"))
  cat("\n--- Within-stratum Δ from BASELINE ---\n"); print(delta_emm)
  diff_emm  <- suppressMessages(contrast(delta_emm, "pairwise", by = "contrast", adjust = "none"))
  cat("\n--- Cross-stratum Δ-difference ---\n"); print(diff_emm)
  list(
    emmeans        = as_tibble(as.data.frame(emm)),
    within_stratum = tidy_contrast(delta_emm),   # 95% CI + within-family Holm (v0.79)
    cross_stratum  = tidy_contrast(diff_emm)
  )
}

emm_sbp_unadj <- extract_emmeans(fit_sbp_unadj, "SBP — UNADJUSTED (HTN-naive)")
emm_sbp_adj   <- extract_emmeans(fit_sbp_adj,   "SBP — ADJUSTED (HTN-naive)")
emm_dbp_unadj <- extract_emmeans(fit_dbp_unadj, "DBP — UNADJUSTED (HTN-naive)")
emm_dbp_adj   <- extract_emmeans(fit_dbp_adj,   "DBP — ADJUSTED (HTN-naive)")

# --- Save ---
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

# --- Emit standardized model sizes for Table A.2 (from the fits) ---
emit_sizes(
  imap_dfr(fitted_models, function(fit, key) tibble(
    section = "5.1", role = "sensitivity",
    analysis = "BP-trajectory LMM (HTN-naive subset)", cohort = "htn_naive",
    outcome = if (startsWith(key, "sbp")) "SBP" else "DBP",
    predictor_set = if (grepl("unadj", key)) "unadjusted" else "adjusted",
    n_users = as.integer(lme4::ngrps(fit)), n_obs = as.integer(nobs(fit)),
    n_events = NA_integer_, ceiling_n = nrow(analysis_users)
  )),
  "5_1_htn_naive", results_dir
)

anova_table <- imap_dfr(fitted_models, function(fit, name) {
  a <- suppressMessages(anova(fit, type = 3))
  as_tibble(a, rownames = "term") %>% mutate(model = name) %>% select(model, everything())
})
write_csv(anova_table, file.path(results_dir, "anova_tables.csv"))

r2_table <- imap_dfr(fitted_models, function(fit, name) {
  r2 <- suppressMessages(r2_nakagawa(fit))
  tibble(model = name,
         marginal_r2    = as.numeric(r2$R2_marginal),
         conditional_r2 = as.numeric(r2$R2_conditional))
})
write_csv(r2_table, file.path(results_dir, "r2.csv"))

write_csv(imap_dfr(emmeans_results, function(r, name) r$emmeans %>% mutate(model = name) %>% select(model, everything())),
          file.path(results_dir, "emmeans.csv"))
write_csv(imap_dfr(emmeans_results, function(r, name) r$within_stratum %>% mutate(model = name) %>% select(model, everything())),
          file.path(results_dir, "delta_contrasts.csv"))
write_csv(imap_dfr(emmeans_results, function(r, name) r$cross_stratum %>% mutate(model = name) %>% select(model, everything())),
          file.path(results_dir, "cross_stratum_diffs.csv"))

cat("\n\n=== §5.1 HTN-naive sensitivity done ===\n")
cat("Saved rds + CSVs to:\n  ", out_dir, "\n")
sink()
