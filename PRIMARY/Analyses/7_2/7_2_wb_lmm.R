# 7_2_wb_lmm.R
# SAP §7.2: Daily WHOOP signals predicting same-day BP.
#
# Two primary LMMs (one per outcome):
#   SBP − baseline_SBP  ~ W/B-decomposed daily predictors + covariates + (1|STUDY_ID)
#   DBP − baseline_DBP  ~ W/B-decomposed daily predictors + covariates + (1|STUDY_ID)
#
# Predictors (5 daily signals, each W/B-decomposed = 10 fixed effects):
#   - RHR_LOOKBACK_MEAN
#   - HRV_LOOKBACK_MEAN
#   - STRAIN_LOOKBACK_MEAN
#   - SLEEP_DURATION_LOOKBACK_MEAN
#   - SLEEP_CONSISTENCY_LOOKBACK_MEAN
# (All are 7-day lookback means computed in the BP-reading-level build.)
#
# W/B decomposition per predictor X:
#   X_b (between) = user-mean of X across the user's non-baseline readings
#   X_w (within)  = X − X_b (deviation from the user's typical level at each reading)
# Both X_b and X_w are standardized (z-scored on the analytic cohort) so that
# coefficients are interpretable as "mmHg per SD of predictor" and directly
# comparable across the 5 daily predictors.
#
# Covariates: age, sex, baseline_BMI.
# Random effects: random intercept per STUDY_ID.
# Standard errors: cluster-robust on STUDY_ID (via clubSandwich::vcovCR, CR1 type).
# (CR1 chosen over CR2 because CR2's per-cluster bias-correction matrices OOM
# at this scale — n_clusters ≈ 17k. At large cluster counts the CR1 vs CR2
# inference difference is negligible; CR2's small-sample correction matters
# most when n_clusters < ~50. SAP §7.2 specifies "cluster-robust on STUDY_ID";
# both CR1 and CR2 satisfy this — the choice is implementation detail.)
#
# Cohort: §2.2 hypertensive + ≥2 non-baseline manual-cuff BP readings during
# the study window.
#
# Per SAP §9: α = 0.005 formal-test threshold; 95% CIs retained as the effect-
# size convention.
#
# Input:  Build/derived/qc/bp_reading_analytic_frame_1X.rds
#
# Outputs (Analyses/7_2/derived/):
#   fitted_models.rds                  — 2 lmerMod objects (SBP, DBP)
#   results/7_2_coefficients.csv       — β (per SD) + 95% CI + p (cluster-robust)
#   results/7_2_correlation_matrix.csv — pairwise correlation among 10 W/B predictors
#   results/7_2_model_summaries.csv    — n rows, n users, R², AIC
#   results/7_2_paper.md               — paper-quote-ready text
#   results/7_2.log                    — full sink log
# (figures/7_2_figure_4.png is generated downstream by the standalone
#  7_2_figure_4.R, which reads results/7_2_coefficients.csv.)

library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)
library(clubSandwich)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc",
                         "bp_reading_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "7_2", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(out_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "7_2.log")
sink(log_path, split = TRUE)
cat("§7.2 W/B-decomposed daily LMM —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read frame ---
bp_frame <- read_rds(input_path)
cat("BP-reading frame:", nrow(bp_frame), "rows ×", ncol(bp_frame), "cols\n\n")

# --- Predictor definitions ---
daily_predictors <- c(
  RHR              = "RHR_LOOKBACK_MEAN",
  HRV              = "HRV_LOOKBACK_MEAN",
  Strain           = "STRAIN_LOOKBACK_MEAN",
  SleepDuration    = "SLEEP_DURATION_LOOKBACK_MEAN",
  SleepConsistency = "SLEEP_CONSISTENCY_LOOKBACK_MEAN"
)

# Favorable direction for the responder-like interpretation (mmHg drop = favorable):
#   "lower" predictor → favorable: more of this predictor associates with greater BP drop (β negative)
#   "higher" predictor → favorable: more of this predictor associates with greater BP drop (β negative)
# Actually for ALL daily signals where "more is better" cardiologically, we expect β < 0.
# For RHR where "less is better", we expect β > 0 (higher RHR → less BP drop).
predictor_favorable_beta <- c(
  RHR              = "positive",  # higher RHR → less BP drop (less negative Δ)
  HRV              = "negative",  # higher HRV → more BP drop (more negative Δ)
  Strain           = "negative",
  SleepDuration    = "negative",
  SleepConsistency = "negative"
)

# --- Cohort filter ---
# §2.2 + non-baseline readings + PERIOD == "DURING" (study window) + ≥2 readings per user
cohort_pre <- bp_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female"),
         !IS_BASELINE_READING,
         PERIOD == "DURING")
ceiling_n_72 <- n_distinct(cohort_pre$STUDY_ID)
cohort <- cohort_pre %>%
  group_by(STUDY_ID) %>%
  filter(n() >= 2) %>%
  ungroup()
cat("Cohort after §2.2 + non-baseline + DURING + ≥2-readings filter:",
    nrow(cohort), "rows from", n_distinct(cohort$STUDY_ID), "users\n\n")

# --- Build W/B decomposition per predictor (computed on the analytic cohort) ---
build_wb <- function(d, pred_col, prefix) {
  d %>%
    group_by(STUDY_ID) %>%
    mutate(!!paste0(prefix, "_b") := mean(.data[[pred_col]], na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(!!paste0(prefix, "_w") := .data[[pred_col]] - .data[[paste0(prefix, "_b")]])
}

cohort_wb <- cohort
for (m in names(daily_predictors)) {
  cohort_wb <- build_wb(cohort_wb, daily_predictors[[m]], m)
}
cohort_wb <- cohort_wb %>%
  mutate(HTN_BASELINE_ON = htn_baseline_on(HTN_MED_BASELINE))

# Complete-case on all predictors + covariates + outcomes
required_cols <- c(
  unlist(lapply(names(daily_predictors), function(m) c(paste0(m, "_b"), paste0(m, "_w")))),
  "AGE_AT_START", "SEX", "BASELINE_BMI", "HTN_BASELINE_ON",
  "SBP_MINUS_BASELINE", "DBP_MINUS_BASELINE"
)
cohort_wb <- cohort_wb %>%
  filter(if_all(all_of(required_cols), ~ !is.na(.))) %>%
  mutate(SEX = factor(SEX))

cat("After complete-case on all predictors + covariates + outcomes:",
    nrow(cohort_wb), "rows from", n_distinct(cohort_wb$STUDY_ID), "users\n\n")

# --- Standardize the 10 W/B predictors ---
# _w (within-person deviation from user mean): row-level z-score is appropriate
#     because _w varies within user; row-level SD reflects the spread of within-
#     person deviations across all reading-rows.
# _b (between-person user mean): user-level z-score (one row per user) because
#     _b is a user-level constant. Row-level z-scoring would weight users with
#     more readings, biasing the SD toward high-engagement users.
#     Per audit finding 7.2.4 + SAP_CHANGELOG v0.63.
wb_terms <- c()
for (m in names(daily_predictors)) {
  b_col <- paste0(m, "_b")
  w_col <- paste0(m, "_w")
  # _b: user-level z-score (one row per user, then map back to rows)
  user_b <- distinct(cohort_wb[, c("STUDY_ID", b_col)])
  b_mean <- mean(user_b[[b_col]], na.rm = TRUE)
  b_sd   <- sd(user_b[[b_col]],   na.rm = TRUE)
  cohort_wb[[b_col]] <- (cohort_wb[[b_col]] - b_mean) / b_sd
  # _w: row-level z-score (deviation from user mean varies within user)
  cohort_wb[[w_col]] <- as.numeric(scale(cohort_wb[[w_col]]))
  wb_terms <- c(wb_terms, b_col, w_col)
}
cat("Standardized W/B predictors (_b: user-level SD; _w: row-level SD):\n  ",
    paste(wb_terms, collapse = ", "), "\n\n")

# --- Pairwise correlation matrix among the 10 W/B predictors ---
cor_mat <- cor(cohort_wb[, wb_terms], use = "pairwise.complete.obs")
cor_long <- as_tibble(cor_mat, rownames = "predictor1") %>%
  pivot_longer(-predictor1, names_to = "predictor2", values_to = "r")
write_csv(cor_long, file.path(results_dir, "7_2_correlation_matrix.csv"))

cat("--- High pairwise correlations among W/B predictors (|r| > 0.4) ---\n")
print(cor_long %>%
        filter(predictor1 != predictor2, abs(r) > 0.4) %>%
        arrange(desc(abs(r))) %>%
        head(20))

# --- Fit primary LMMs (SBP-baseline, DBP-baseline) ---
fit_wb_lmm <- function(outcome_col, data) {
  rhs <- paste(c("AGE_AT_START", "SEX", "BASELINE_BMI", "HTN_BASELINE_ON", wb_terms),
               collapse = " + ")
  f <- as.formula(paste0(outcome_col, " ~ ", rhs, " + (1 | STUDY_ID)"))
  cat(sprintf("Fitting LMM: %s ~ ...\n", outcome_col))
  t0 <- Sys.time()
  fit <- lmer(f, data = data, REML = TRUE)
  cat(sprintf("  Elapsed: %.1f min\n",
              as.numeric(Sys.time() - t0, units = "mins")))
  fit
}

fit_sbp <- fit_wb_lmm("SBP_MINUS_BASELINE", cohort_wb)
fit_dbp <- fit_wb_lmm("DBP_MINUS_BASELINE", cohort_wb)

# --- Checkpoint: save LMM fits + cohort_wb BEFORE the slow clubSandwich step ---
# If the clubSandwich CR2 extraction below is too slow and gets interrupted,
# fitted_models.rds and cohort_wb_checkpoint.rds let a recovery script pick
# up from here without refitting the LMMs.
write_rds(list(SBP = fit_sbp, DBP = fit_dbp), file.path(out_dir, "fitted_models.rds"))
write_rds(cohort_wb, file.path(out_dir, "cohort_wb_checkpoint.rds"))
cat("\n[Checkpoint] Saved fitted_models.rds + cohort_wb_checkpoint.rds before clubSandwich step.\n")

# --- Extract coefficients with cluster-robust SEs (CR1) ---
# CR1 chosen over CR2 because CR2 OOM'd at this scale (~17k clusters); the
# small-sample bias correction in CR2 is most meaningful when n_clusters < ~50.
# At large cluster counts, CR1 vs CR2 inference is functionally identical.
# z-test inference (standard for large-n_clusters cluster-robust).
extract_robust_coefs <- function(fit, outcome_label, data) {
  cat(sprintf("\nExtracting cluster-robust coefficients (CR1) for %s...\n", outcome_label))
  t0 <- Sys.time()
  vcov_cr <- clubSandwich::vcovCR(fit, cluster = data$STUDY_ID, type = "CR1")
  cat(sprintf("  vcovCR1 elapsed: %.2f min\n",
              as.numeric(Sys.time() - t0, units = "mins")))
  ct <- clubSandwich::coef_test(fit, vcov = vcov_cr,
                                test = "z", coefs = "All")
  ci <- clubSandwich::conf_int(fit, vcov = vcov_cr,
                               test = "z", coefs = "All", level = 0.95)
  tibble(
    outcome    = outcome_label,
    term       = rownames(ct),
    estimate   = ct$beta,
    SE_robust  = ct$SE,
    z_stat     = ct$tstat,
    p_value    = ct$p_z,
    conf_low   = ci$CI_L,
    conf_high  = ci$CI_U
  ) %>%
    mutate(sig_005 = p_value < 0.005)
}

coefs_sbp <- extract_robust_coefs(fit_sbp, "SBP-baseline", cohort_wb)
coefs_dbp <- extract_robust_coefs(fit_dbp, "DBP-baseline", cohort_wb)
coefs_all <- bind_rows(coefs_sbp, coefs_dbp)
write_csv(coefs_all, file.path(results_dir, "7_2_coefficients.csv"))

cat("\n--- §7.2 coefficients (W/B daily predictors, standardized β per SD) ---\n")
print(coefs_all %>%
        filter(term %in% wb_terms) %>%
        select(outcome, term, estimate, SE_robust, conf_low, conf_high,
               p_value, sig_005),
      n = Inf, width = Inf)

# --- Model summaries ---
model_summ <- bind_rows(
  tibble(model = "SBP-baseline", n_obs = nobs(fit_sbp),
         n_users = lme4::ngrps(fit_sbp),
         aic = AIC(fit_sbp),
         marginal_r2    = as.numeric(performance::r2_nakagawa(fit_sbp)$R2_marginal),
         conditional_r2 = as.numeric(performance::r2_nakagawa(fit_sbp)$R2_conditional)),
  tibble(model = "DBP-baseline", n_obs = nobs(fit_dbp),
         n_users = lme4::ngrps(fit_dbp),
         aic = AIC(fit_dbp),
         marginal_r2    = as.numeric(performance::r2_nakagawa(fit_dbp)$R2_marginal),
         conditional_r2 = as.numeric(performance::r2_nakagawa(fit_dbp)$R2_conditional))
)
write_csv(model_summ, file.path(results_dir, "7_2_model_summaries.csv"))

cat("\n--- Model summaries ---\n")
print(model_summ)

# --- emit_sizes: §7.2 primary ---
sizes_72 <- tibble(
  section       = "7.2",
  role          = "primary",
  analysis      = "W/B daily LMM (anchored BP − baseline_BP)",
  cohort        = "hypertensive",
  outcome       = c("SBP", "DBP"),
  predictor_set = "anchored",
  n_users       = as.integer(n_distinct(cohort_wb$STUDY_ID)),
  n_obs         = as.integer(nrow(cohort_wb)),
  n_events      = NA_integer_,
  ceiling_n     = as.integer(ceiling_n_72)
)
emit_sizes(sizes_72, "7_2", results_dir)

# --- Figure 4: forest plot of standardized β ---
# Generated by the standalone Analyses/7_2/7_2_figure_4.R, which reads the
# coefficient CSV written just above (results/7_2_coefficients.csv) and writes
# figures/7_2_figure_4.png. Kept separate so the figure can be restyled without
# refitting the cluster-robust LMM — the same figure-script separation used in
# §5.1 (5_1_figure_1a.R) and §5.5 (5_5_figure_2.R). After a full re-run of this
# script, regenerate the figure with:  Rscript PRIMARY/Analyses/7_2/7_2_figure_4.R

# --- Direction-favorability check ---
# For each (form × outcome × metric):
#   "negative" favorable → β should be < 0 for "more predictor → more BP drop" (i.e., reduction)
#   "positive" favorable → β should be > 0 for the predictor where "less is better" (RHR)
direction_check <- coefs_all %>%
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
  select(outcome, form, metric, estimate, p_value, sig_005,
         favorable_direction, favorable_concordant)
write_csv(direction_check, file.path(results_dir, "7_2_direction_check.csv"))

cat("\n--- Direction-favorability check ---\n")
print(direction_check, n = Inf, width = Inf)

cat(sprintf("\n  → %d/%d coefficients direction-favorable for BP reduction.\n",
            sum(direction_check$favorable_concordant), nrow(direction_check)))

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
  "# §7.2 Table — Daily WHOOP signals predicting same-day BP",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  paste0("§2.2 cohort, n = ", format(nrow(cohort_wb), big.mark = ","),
         " readings from ", format(n_distinct(cohort_wb$STUDY_ID), big.mark = ","),
         " users."),
  "",
  "Standardized β = mmHg of BP change from baseline per SD of predictor. ",
  "Bold = significant at α = 0.005 (cluster-robust CR1, z-test)."
)
for (outcome_label in c("SBP-baseline", "DBP-baseline")) {
  md_lines <- c(md_lines, "", paste0("## ", outcome_label))
  for (form in c("Between-person", "Within-person")) {
    md_lines <- c(md_lines, "", paste0("### ", form))
    form_suffix <- ifelse(form == "Between-person", "_b", "_w")
    rows <- coefs_all %>%
      filter(outcome == outcome_label, grepl(paste0(form_suffix, "$"), term))
    for (m in names(daily_predictors)) {
      r <- rows %>% filter(term == paste0(m, form_suffix))
      if (nrow(r) == 0) next
      md_lines <- c(md_lines, sprintf(
        "- **%s**: %s",
        m, fmt_b(r$estimate, r$conf_low, r$conf_high, r$p_value)
      ))
    }
  }
}

md_lines <- c(md_lines,
  "",
  "## Notes",
  "",
  "- **Outcome**: `SBP − baseline_SBP` and `DBP − baseline_DBP` at each non-baseline reading. One row per non-baseline BP reading per user.",
  "- **W/B decomposition** per predictor: between-person = user-mean across the user's non-baseline readings; within-person = reading-level deviation from user's own mean. Both standardized (z-scored on the analytic cohort).",
  "- **Cluster-robust SEs**: clubSandwich::vcovCR, CR1 type, STUDY_ID clusters; z-test inference. (CR2 was originally specified but OOM'd at this scale; at n_clusters ≈ 17k the CR1 vs CR2 inference difference is negligible — see SAP_CHANGELOG if applicable.)",
  "- **Random effects**: random intercept per STUDY_ID (level 1 = readings; level 2 = users).",
  "- **Cohort**: §2.2 hypertensive + ≥2 non-baseline readings during study window.",
  "- **Coefficient interpretation**: β = expected change in (BP − baseline_BP) per 1-SD increase in the predictor. Negative β = more of this predictor → more BP reduction. For RHR specifically, positive β = higher RHR → less BP reduction (less favorable).",
  "- **Within-person coefficients are the §7.2 distinctive contribution**: each user serves as their own control for stable confounders.",
  "- Pairs with: §7.1 trajectory analyses (group-level WHOOP trajectories); §6.1/§6.2 (baseline/study-mean predictors of BP outcomes)."
)
writeLines(md_lines, file.path(results_dir, "7_2_paper.md"))
cat("\nWrote markdown:\n  ", file.path(results_dir, "7_2_paper.md"), "\n")

cat("\n=== §7.2 W/B LMM done ===\n")
sink()
