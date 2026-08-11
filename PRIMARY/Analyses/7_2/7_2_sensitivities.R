# 7_2_sensitivities.R
# SAP §7.2 sensitivity analyses, three SAP-prespecified:
#
#   (1) Raw-BP sensitivity. Refit the §7.2 primary LMM with raw SBP / DBP as
#       outcome and baseline_BP re-introduced as a between-person covariate.
#       Addresses RTM: convergence between anchored and raw-BP coefficients in
#       direction and rough magnitude = empirical evidence RTM is not driving
#       findings; divergence flags RTM as candidate explanation.
#
#   (2) Full-cohort regression (including Normal stratum). Same model structure,
#       fit on all four AHA strata with the ≥2-non-baseline-readings filter.
#       Tests robustness to the hypertensive-only inclusion criterion.
#
#   (3) HTN-medication-naive subset replication. Refit the §7.2 LMM on the
#       subset with never_logged HTN-medication status. Tests whether daily-
#       signal associations are driven by medication-treated users.
#
# Each sensitivity refits two LMMs (SBP, DBP outcomes) = 6 LMMs total.
# Cluster-robust SEs via clubSandwich::vcovCR (CR1; CR2 OOM at this scale).
#
# Input:  Build/derived/qc/bp_reading_analytic_frame_1X.rds
#
# Outputs (Analyses/7_2/derived/):
#   results/7_2_sensitivity_coefficients.csv     — long-format: sensitivity × outcome × term → β + CI + p
#   results/7_2_sensitivity_comparison.csv       — side-by-side vs §7.2 primary
#   results/7_2_sensitivity_direction_check.csv  — favorable-direction concordance per cell
#   results/7_2_sensitivities_paper.md           — paper-quote-ready markdown
#   results/7_2_sensitivities.log                — full sink log

library(tidyverse)
library(lme4)
library(lmerTest)
library(clubSandwich)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc",
                         "bp_reading_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "7_2", "derived")
results_dir <- file.path(out_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "7_2_sensitivities.log")
sink(log_path, split = TRUE)
cat("§7.2 sensitivities — raw-BP + full-cohort + HTN-naive —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read frame ---
bp_frame <- read_rds(input_path)
cat("BP-reading frame:", nrow(bp_frame), "rows ×", ncol(bp_frame), "cols\n\n")

# Ceiling counts (users before ≥2-readings filter, per cohort definition)
ceiling_n_primary <- n_distinct(
  bp_frame %>%
    filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
           SEX %in% c("male", "female"),
           !IS_BASELINE_READING, PERIOD == "DURING") %>%
    pull(STUDY_ID)
)
ceiling_n_full <- n_distinct(
  bp_frame %>%
    filter(BASELINE_AHA %in% c("Normal", "Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
           SEX %in% c("male", "female"),
           !IS_BASELINE_READING, PERIOD == "DURING") %>%
    pull(STUDY_ID)
)
ceiling_n_htn_naive <- n_distinct(
  bp_frame %>%
    filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
           SEX %in% c("male", "female"),
           HTN_MED_BASELINE == "never_logged",
           !IS_BASELINE_READING, PERIOD == "DURING") %>%
    pull(STUDY_ID)
)

# --- Predictor definitions (match §7.2 primary) ---
daily_predictors <- c(
  RHR              = "RHR_LOOKBACK_MEAN",
  HRV              = "HRV_LOOKBACK_MEAN",
  Strain           = "STRAIN_LOOKBACK_MEAN",
  SleepDuration    = "SLEEP_DURATION_LOOKBACK_MEAN",
  SleepConsistency = "SLEEP_CONSISTENCY_LOOKBACK_MEAN"
)

predictor_favorable_beta <- c(
  RHR              = "positive",  # higher RHR → less BP drop (positive β favorable)
  HRV              = "negative",  # higher HRV → more BP drop (negative β favorable)
  Strain           = "negative",
  SleepDuration    = "negative",
  SleepConsistency = "negative"
)

# --- Cohort builders ---
# §2.2 hypertensive + non-baseline + PERIOD == "DURING" + ≥2 readings + sex restriction (primary)
build_primary_cohort <- function(d) {
  d %>%
    filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
           SEX %in% c("male", "female"),
           !IS_BASELINE_READING,
           PERIOD == "DURING") %>%
    group_by(STUDY_ID) %>%
    filter(n() >= 2) %>%
    ungroup()
}

# Full-cohort (4 strata)
build_full_cohort <- function(d) {
  d %>%
    filter(BASELINE_AHA %in% c("Normal", "Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
           SEX %in% c("male", "female"),
           !IS_BASELINE_READING,
           PERIOD == "DURING") %>%
    group_by(STUDY_ID) %>%
    filter(n() >= 2) %>%
    ungroup()
}

# HTN-naive
build_htn_naive_cohort <- function(d) {
  build_primary_cohort(d) %>%
    filter(HTN_MED_BASELINE == "never_logged")
}

# --- Helper: build W/B decomposition (per cohort) ---
build_wb <- function(d, pred_col, prefix) {
  d %>%
    group_by(STUDY_ID) %>%
    mutate(!!paste0(prefix, "_b") := mean(.data[[pred_col]], na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(!!paste0(prefix, "_w") := .data[[pred_col]] - .data[[paste0(prefix, "_b")]])
}

# --- Helper: prep cohort (W/B decompose, complete-case, standardize) ---
prep_cohort <- function(d, include_baseline_bp = FALSE) {
  for (m in names(daily_predictors)) {
    d <- build_wb(d, daily_predictors[[m]], m)
  }
  wb_cols <- c(unlist(lapply(names(daily_predictors),
                              function(m) c(paste0(m, "_b"), paste0(m, "_w")))))
  required <- c(wb_cols, "AGE_AT_START", "SEX", "BASELINE_BMI",
                "SBP_MINUS_BASELINE", "DBP_MINUS_BASELINE",
                "SBP", "DBP", "BASELINE_SBP", "BASELINE_DBP")
  d <- d %>%
    filter(if_all(all_of(required), ~ !is.na(.))) %>%
    mutate(SEX = factor(SEX),
           HTN_BASELINE_ON = htn_baseline_on(HTN_MED_BASELINE))
  # Standardize: _b on user-level distribution (one row per user); _w row-level.
  # Per audit finding 7.2.4 + SAP_CHANGELOG v0.63.
  for (m in names(daily_predictors)) {
    b_col <- paste0(m, "_b")
    w_col <- paste0(m, "_w")
    user_b <- distinct(d[, c("STUDY_ID", b_col)])
    b_mean <- mean(user_b[[b_col]], na.rm = TRUE)
    b_sd   <- sd(user_b[[b_col]],   na.rm = TRUE)
    d[[b_col]] <- (d[[b_col]] - b_mean) / b_sd
    d[[w_col]] <- as.numeric(scale(d[[w_col]]))
  }
  d
}

# --- Helper: fit LMM + extract CR1 coefficients ---
fit_and_extract <- function(formula_str, data, label) {
  cat(sprintf("\n--- %s ---\n", label))
  cat("n_obs:", nrow(data), "; n_users:", n_distinct(data$STUDY_ID), "\n")
  f <- as.formula(formula_str)
  t0 <- Sys.time()
  fit <- lmer(f, data = data, REML = TRUE)
  cat(sprintf("  LMM fit: %.2f min\n", as.numeric(Sys.time() - t0, units = "mins")))
  t0 <- Sys.time()
  vcov_cr <- clubSandwich::vcovCR(fit, cluster = data$STUDY_ID, type = "CR1")
  cat(sprintf("  vcovCR1: %.2f min\n", as.numeric(Sys.time() - t0, units = "mins")))
  ct <- clubSandwich::coef_test(fit, vcov = vcov_cr,
                                test = "z", coefs = "All")
  ci <- clubSandwich::conf_int(fit, vcov = vcov_cr,
                               test = "z", coefs = "All", level = 0.95)
  tibble(
    label      = label,
    term       = rownames(ct),
    estimate   = ct$beta,
    SE_robust  = ct$SE,
    z_stat     = ct$tstat,
    p_value    = ct$p_z,
    conf_low   = ci$CI_L,
    conf_high  = ci$CI_U
  ) %>%
    mutate(sig_005 = p_value < 0.005,
           n_obs = nrow(data),
           n_users = n_distinct(data$STUDY_ID))
}

# --- Build sensitivity cohorts (each W/B-decomposed within its own cohort) ---
cat("Preparing cohorts...\n")
cohort_primary  <- prep_cohort(build_primary_cohort(bp_frame))
cohort_full     <- prep_cohort(build_full_cohort(bp_frame))
cohort_htnnaive <- prep_cohort(build_htn_naive_cohort(bp_frame))
cat(sprintf("  Primary (§2.2): %d readings, %d users\n",
            nrow(cohort_primary), n_distinct(cohort_primary$STUDY_ID)))
cat(sprintf("  Full cohort (4 strata): %d readings, %d users\n",
            nrow(cohort_full), n_distinct(cohort_full$STUDY_ID)))
cat(sprintf("  HTN-naive: %d readings, %d users\n",
            nrow(cohort_htnnaive), n_distinct(cohort_htnnaive$STUDY_ID)))

# --- Build W/B term names (consistent across cohorts) ---
wb_terms <- c()
for (m in names(daily_predictors)) {
  for (form in c("_b", "_w")) {
    wb_terms <- c(wb_terms, paste0(m, form))
  }
}
covariates_primary  <- c("AGE_AT_START", "SEX", "BASELINE_BMI", "HTN_BASELINE_ON")
covariates_htn_naive <- c("AGE_AT_START", "SEX", "BASELINE_BMI")

# --- Fit all sensitivity LMMs (6 total: 3 sensitivities × 2 outcomes) ---
all_coefs <- list()

cat("\n\n========================================\n")
cat("SENSITIVITY 1: Raw-BP outcome (baseline_BP as covariate)\n")
cat("========================================\n")
# Outcome: raw SBP / DBP. Covariate: baseline_BP as between-person.
formula_sbp_raw <- paste0(
  "SBP ~ BASELINE_SBP + ",
  paste(c(covariates_primary, wb_terms), collapse = " + "),
  " + (1 | STUDY_ID)"
)
formula_dbp_raw <- paste0(
  "DBP ~ BASELINE_DBP + ",
  paste(c(covariates_primary, wb_terms), collapse = " + "),
  " + (1 | STUDY_ID)"
)
all_coefs[["RawBP__SBP"]] <- fit_and_extract(formula_sbp_raw, cohort_primary, "RawBP__SBP") %>%
  mutate(sensitivity = "RawBP", outcome = "SBP")
all_coefs[["RawBP__DBP"]] <- fit_and_extract(formula_dbp_raw, cohort_primary, "RawBP__DBP") %>%
  mutate(sensitivity = "RawBP", outcome = "DBP")

cat("\n\n========================================\n")
cat("SENSITIVITY 2: Full cohort (4 strata including Normal)\n")
cat("========================================\n")
formula_sbp_full <- paste0(
  "SBP_MINUS_BASELINE ~ ",
  paste(c(covariates_primary, wb_terms), collapse = " + "),
  " + (1 | STUDY_ID)"
)
formula_dbp_full <- paste0(
  "DBP_MINUS_BASELINE ~ ",
  paste(c(covariates_primary, wb_terms), collapse = " + "),
  " + (1 | STUDY_ID)"
)
all_coefs[["FullCohort__SBP"]] <- fit_and_extract(formula_sbp_full, cohort_full, "FullCohort__SBP") %>%
  mutate(sensitivity = "FullCohort", outcome = "SBP")
all_coefs[["FullCohort__DBP"]] <- fit_and_extract(formula_dbp_full, cohort_full, "FullCohort__DBP") %>%
  mutate(sensitivity = "FullCohort", outcome = "DBP")

cat("\n\n========================================\n")
cat("SENSITIVITY 3: HTN-medication-naive subset\n")
cat("========================================\n")
formula_sbp_htn <- paste0(
  "SBP_MINUS_BASELINE ~ ",
  paste(c(covariates_htn_naive, wb_terms), collapse = " + "),
  " + (1 | STUDY_ID)"
)
formula_dbp_htn <- paste0(
  "DBP_MINUS_BASELINE ~ ",
  paste(c(covariates_htn_naive, wb_terms), collapse = " + "),
  " + (1 | STUDY_ID)"
)
all_coefs[["HTN_naive__SBP"]] <- fit_and_extract(formula_sbp_htn, cohort_htnnaive, "HTN_naive__SBP") %>%
  mutate(sensitivity = "HTN_naive", outcome = "SBP")
all_coefs[["HTN_naive__DBP"]] <- fit_and_extract(formula_dbp_htn, cohort_htnnaive, "HTN_naive__DBP") %>%
  mutate(sensitivity = "HTN_naive", outcome = "DBP")

# --- Combine + write CSVs ---
coefs_df <- bind_rows(all_coefs) %>%
  select(sensitivity, outcome, term, estimate, SE_robust, conf_low, conf_high,
         z_stat, p_value, sig_005, n_obs, n_users)
write_csv(coefs_df, file.path(results_dir, "7_2_sensitivity_coefficients.csv"))

# --- emit_sizes: §7.2 sensitivities ---
sens_meta <- list(
  RawBP      = list(analysis = "W/B daily LMM (raw-BP)",      cohort = "hypertensive", ceiling = ceiling_n_primary),
  FullCohort = list(analysis = "W/B daily LMM (full-cohort)", cohort = "full",          ceiling = ceiling_n_full),
  HTN_naive  = list(analysis = "W/B daily LMM (HTN-naive)",   cohort = "htn_naive",     ceiling = ceiling_n_htn_naive)
)
sizes_72s <- imap_dfr(sens_meta, function(meta, skey) {
  rows <- coefs_df %>%
    filter(sensitivity == skey) %>%
    distinct(outcome, n_obs, n_users)
  tibble(
    section       = "7.2",
    role          = "sensitivity",
    analysis      = meta$analysis,
    cohort        = meta$cohort,
    outcome       = rows$outcome,
    predictor_set = "anchored",
    n_users       = as.integer(rows$n_users),
    n_obs         = as.integer(rows$n_obs),
    n_events      = NA_integer_,
    ceiling_n     = as.integer(meta$ceiling)
  )
})
emit_sizes(sizes_72s, "7_2_sensitivities", results_dir)

cat("\n\n--- Sensitivity coefficients (W/B daily predictors only) ---\n")
print(coefs_df %>%
        filter(term %in% wb_terms) %>%
        select(sensitivity, outcome, term, estimate, conf_low, conf_high,
               p_value, sig_005),
      n = Inf, width = Inf)

# --- Compare to §7.2 primary ---
primary <- read_csv(file.path(results_dir, "7_2_coefficients.csv"),
                    show_col_types = FALSE) %>%
  filter(term %in% wb_terms) %>%
  transmute(outcome = sub("-baseline", "", outcome),
            term, primary_estimate = estimate, primary_p = p_value,
            primary_sig = sig_005)

comparison <- coefs_df %>%
  filter(term %in% wb_terms) %>%
  left_join(primary, by = c("outcome", "term")) %>%
  mutate(
    diff_vs_primary = estimate - primary_estimate,
    sig_concordant  = sig_005 == primary_sig
  ) %>%
  arrange(sensitivity, outcome, term)

write_csv(comparison, file.path(results_dir, "7_2_sensitivity_comparison.csv"))

cat("\n--- Sensitivity vs primary comparison (significant predictors in primary only) ---\n")
print(comparison %>%
        filter(primary_sig) %>%
        select(sensitivity, outcome, term,
               sens_β = estimate, primary_β = primary_estimate,
               diff_vs_primary, sens_p = p_value, sens_sig = sig_005),
      n = Inf, width = Inf)

# --- Convergence summary ---
convergence <- comparison %>%
  filter(primary_sig) %>%
  group_by(sensitivity) %>%
  summarise(
    n_cells              = n(),
    n_sig_in_sensitivity = sum(sig_005),
    n_direction_concord  = sum(sign(estimate) == sign(primary_estimate)),
    .groups = "drop"
  )

cat("\n--- Convergence summary (primary-significant predictors only) ---\n")
print(convergence)

# --- Direction-favorability check ---
direction_check <- coefs_df %>%
  filter(term %in% wb_terms) %>%
  mutate(
    metric = sub("_(b|w)$", "", term),
    form   = ifelse(grepl("_b$", term), "between", "within"),
    favorable_direction = predictor_favorable_beta[metric],
    favorable_concordant = case_when(
      favorable_direction == "negative" ~ estimate < 0,
      favorable_direction == "positive" ~ estimate > 0,
      TRUE                              ~ NA
    )
  ) %>%
  select(sensitivity, outcome, form, metric, estimate, p_value, sig_005,
         favorable_direction, favorable_concordant)
write_csv(direction_check, file.path(results_dir, "7_2_sensitivity_direction_check.csv"))

direction_summary <- direction_check %>%
  group_by(sensitivity, form) %>%
  summarise(n_metrics = n(),
            n_favorable = sum(favorable_concordant, na.rm = TRUE),
            n_sig = sum(sig_005),
            .groups = "drop")
cat("\n--- Direction-favorability summary per sensitivity × form ---\n")
print(direction_summary)

# --- Paper-ready markdown ---
fmt_b <- function(est, lcl, ucl, p) {
  sig <- ifelse(p < 0.005, "**", "")
  sprintf("%s%+.3f (%.3f, %.3f)%s [p = %s]",
          sig, est, lcl, ucl, sig,
          ifelse(p < 1e-10, formatC(p, format = "e", digits = 1),
                 ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                        sprintf("%.3f", p))))
}

md_lines <- c(
  "# §7.2 sensitivity analyses — raw-BP + full-cohort + HTN-naive",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "Three SAP-prespecified sensitivities applied to the §7.2 W/B-decomposed daily LMM. All use CR1 cluster-robust SEs on STUDY_ID (CR2 OOM at this scale).",
  "",
  "## Convergence summary (primary-significant predictors only)",
  ""
)
for (i in seq_len(nrow(convergence))) {
  r <- convergence[i, ]
  md_lines <- c(md_lines, sprintf(
    "- **%s**: %d/%d primary-sig cells direction-concordant; %d/%d also sig at α = 0.005.",
    r$sensitivity,
    r$n_direction_concord, r$n_cells,
    r$n_sig_in_sensitivity, r$n_cells
  ))
}

md_lines <- c(md_lines, "",
  "## Notes",
  "",
  "- **Raw-BP sensitivity** addresses the RTM concern flagged in §7.2 primary (between-person DBP RHR_b coefficient in unfavorable direction, likely due to baseline-severity confound). With `baseline_BP` re-introduced as a between-person covariate, the between-person RHR_b coefficient should change toward zero or reverse direction if RTM was driving it.",
  "- **Full-cohort sensitivity** includes Normal-stratum users (~8k additional users) to test whether the W/B-decomposed findings depend on the §2.2 hypertensive-only inclusion criterion.",
  "- **HTN-naive sensitivity** restricts to users with `never_logged` HTN-medication status. Convergence indicates daily-signal associations are not driven by medication-treated users.",
  "- Cluster-robust SE: CR1 type (CR2 OOM at this scale per §7.2 primary). z-test inference at large n_clusters is standard.",
  "- Pairs with: §7.2 primary (`7_2_wb_lmm.R`)."
)

writeLines(md_lines, file.path(results_dir, "7_2_sensitivities_paper.md"))
cat("\nWrote markdown:\n  ", file.path(results_dir, "7_2_sensitivities_paper.md"), "\n")

cat("\n=== §7.2 sensitivities done ===\n")
sink()
