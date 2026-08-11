# 5_1_primary_lmm.R
# SAP §5.1: primary BP-trajectory LMM (co-primary outcome).
#
# Per SAP §5.1 spec:
#   Outcome:           raw SBP and DBP at baseline, 3mo (optional), 6mo
#   Cohort:            §2.2 hypertensive (Elevated + Stage 1 + Stage 2); descriptive
#                      n = 20,021, sex-restricted modeled n = 20,003 (LMM excludes
#                      non-binary / prefer-not-to-say, per §2.3)
#   Fixed effects:     Time + AHA_baseline_stratum + Time × AHA_baseline_stratum
#   Adjusted variant:  + age + sex + baseline_BMI + baseline_HTN_med_status + time_of_day
#                        (cyclic sin/cos pair per SAP v0.49) + season + Time:age + Time:sex
#                        + Time:baseline_BMI + Time:baseline_HTN_med_status
#                        (baseline HTN status = on-at-baseline trait only, with a Time
#                         interaction so it adjusts the trajectory; mid-study initiation
#                         is deliberately excluded — it is a post-baseline mediator)
#   Random effects:    (1 | STUDY_ID)
#   Time coding:       BASELINE (ref), MO3, MO6
#   AHA reference:     Elevated
#   Estimand extract:  emmeans per (Time × Stratum); within-stratum Δ from baseline at
#                      3mo and 6mo with 95% CIs; cross-stratum Δ-difference contrasts.
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds (Phase 1 output)
#
# Output (printed to console):
#   - 4 model summaries (SBP × {unadjusted, adjusted}; DBP × {unadjusted, adjusted})
#   - Type III ANOVA — focus on Time × AHA interaction
#   - emmeans + Δ contrasts + cross-stratum Δ-difference
#   - Nakagawa R² (marginal + conditional)
#
# Saved objects (rds, for downstream Table 2 / Figure 1A scripts):
#   Analyses/5_1/derived/fitted_models.rds  — list of 4 fitted lmerMod objects
#   Analyses/5_1/derived/emmeans_results.rds — list of emmeans / contrast tibbles

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")  # emit_sizes() + Table A contract

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "5_1", "derived")
results_dir <- file.path(out_dir, "results")
dir.create(out_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Capture full console output to a log file (split = TRUE → also show on console)
log_path <- file.path(results_dir, "5_1_primary_lmm.log")
sink(log_path, split = TRUE)
cat("§5.1 primary LMM run —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read user-level analytic frame ---
user_frame <- read_rds(input_path)
cat("User frame:", nrow(user_frame), "rows ×", ncol(user_frame), "cols\n")

# --- Filter to §2.2 cohort + restrict sex to {male, female} for LMM ---
analysis_users <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female"))
cat("§2.2 cohort filter + sex restriction: ", nrow(analysis_users), "users\n")

# --- Pivot to long format: one row per (user × timepoint) ---
# Wide → long for outcomes (SBP / DBP) and per-timepoint covariates (TIME_OF_DAY, SEASON)
long_outcomes <- analysis_users %>%
  select(STUDY_ID, BASELINE_AHA, AGE_AT_START, SEX, BASELINE_BMI, HTN_MED_BASELINE,
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
  rename(SBP = SBP, DBP = DBP, TIME_OF_DAY = TIME_OF_DAY, SEASON = SEASON)

long_outcomes <- long_outcomes %>%
  mutate(
    Time         = factor(Time, levels = c("BASELINE", "MO3", "MO6")),
    BASELINE_AHA = factor(BASELINE_AHA,
                          levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")),
    SEX          = factor(SEX),
    SEASON       = factor(SEASON, levels = c("winter", "spring", "summer", "fall")),
    # Baseline HTN-med status (strictly the pre-baseline trait): on at baseline =
    # stable_on_meds ∪ recent_initiator; mid-study initiators are NOT on at baseline.
    # Added as an adjusted-model covariate WITH a Time interaction — in a random-
    # intercept model a main-effect-only baseline covariate is absorbed by
    # (1|STUDY_ID) and leaves the within-stratum Δ unchanged; the Time interaction
    # is what adjusts the trajectory for baseline-med composition.
    HTN_BASELINE_ON = factor(
      if_else(HTN_MED_BASELINE %in% c("stable_on_meds", "recent_initiator"),
              "on_at_baseline", "not_on_at_baseline"),
      levels = c("not_on_at_baseline", "on_at_baseline"))
  )

cat("\nLong-format frame: ", nrow(long_outcomes), "rows total\n")
cat("Rows per (Time × stratum) cell:\n")
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
  d <- data %>% filter(!is.na(.data[[outcome_col]]),
                       !is.na(BASELINE_AHA))
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

  cat("--- Type III ANOVA (Satterthwaite df; Time:AHA is the formal §5.1 test) ---\n")
  print(suppressMessages(anova(fit, type = 3)))

  cat("\n--- Nakagawa R² (marginal + conditional) ---\n")
  print(performance::r2_nakagawa(fit))
}

# --- Fit 4 models ---
cat("\n\n========================================\n")
cat("Fitting Unadjusted LMM (SBP)…\n")
cat("========================================\n")
fit_sbp_unadj <- fit_lmm("SBP", adjusted = FALSE, data = long_outcomes)
print_model_summary(fit_sbp_unadj, "SBP — UNADJUSTED")

cat("\n\n========================================\n")
cat("Fitting Adjusted LMM (SBP)…\n")
cat("========================================\n")
fit_sbp_adj <- fit_lmm("SBP", adjusted = TRUE, data = long_outcomes)
print_model_summary(fit_sbp_adj, "SBP — ADJUSTED")

cat("\n\n========================================\n")
cat("Fitting Unadjusted LMM (DBP)…\n")
cat("========================================\n")
fit_dbp_unadj <- fit_lmm("DBP", adjusted = FALSE, data = long_outcomes)
print_model_summary(fit_dbp_unadj, "DBP — UNADJUSTED")

cat("\n\n========================================\n")
cat("Fitting Adjusted LMM (DBP)…\n")
cat("========================================\n")
fit_dbp_adj <- fit_lmm("DBP", adjusted = TRUE, data = long_outcomes)
print_model_summary(fit_dbp_adj, "DBP — ADJUSTED")

# --- Estimand extraction per fit: emmeans, Δ contrasts, cross-stratum Δ-diffs ---
extract_emmeans <- function(fit, label) {
  cat("\n\n=== Estimand extraction:", label, " ===\n")

  emm <- suppressMessages(emmeans(fit, ~ Time | BASELINE_AHA))

  cat("\n--- Estimated marginal means (Time × Stratum) ---\n")
  print(emm)

  cat("\n--- Within-stratum Δ from BASELINE (3mo, 6mo) ---\n")
  delta_emm <- suppressMessages(contrast(emm, method = "trt.vs.ctrl", adjust = "none"))
  print(delta_emm)

  cat("\n--- Cross-stratum Δ-difference (each follow-up timepoint × each stratum pair) ---\n")
  diff_emm <- suppressMessages(contrast(delta_emm, "pairwise", by = "contrast", adjust = "none"))
  print(diff_emm)

  cat("\n--- Within-stratum CONSECUTIVE Δ (adjacent timepoints: MO3-BASELINE, MO6-MO3) ---\n")
  consec_emm <- suppressMessages(contrast(emm, method = "consec", adjust = "none"))
  print(consec_emm)

  # Tidy a contrast set with 95% CIs (infer) + a within-family Holm-adjusted p.
  # Reporting convention (§5.1): ALL contrast families are Holm-adjusted WITHIN
  # their own family — within-stratum Δ, cross-stratum Δ-difference, and
  # consecutive Δ are each corrected over their own ~6 contrasts; families are
  # NOT pooled (they test conceptually distinct questions). The primary
  # within-stratum Δ is led by its 95% CI (estimation-first); p.holm accompanies.
  to_tbl <- function(cobj) {
    df <- as_tibble(summary(cobj, infer = c(TRUE, TRUE), adjust = "none"))
    df %>% mutate(p.holm = p.adjust(p.value, method = "holm"))
  }
  list(
    emmeans          = as_tibble(as.data.frame(emm)),
    within_stratum   = to_tbl(delta_emm),
    cross_stratum    = to_tbl(diff_emm),
    consecutive      = to_tbl(consec_emm)
  )
}

emm_sbp_unadj <- extract_emmeans(fit_sbp_unadj, "SBP — UNADJUSTED")
emm_sbp_adj   <- extract_emmeans(fit_sbp_adj,   "SBP — ADJUSTED")
emm_dbp_unadj <- extract_emmeans(fit_dbp_unadj, "DBP — UNADJUSTED")
emm_dbp_adj   <- extract_emmeans(fit_dbp_adj,   "DBP — ADJUSTED")

# --- Save objects for downstream scripts ---
fitted_models <- list(
  sbp_unadj = fit_sbp_unadj,
  sbp_adj   = fit_sbp_adj,
  dbp_unadj = fit_dbp_unadj,
  dbp_adj   = fit_dbp_adj
)
emmeans_results <- list(
  sbp_unadj = emm_sbp_unadj,
  sbp_adj   = emm_sbp_adj,
  dbp_unadj = emm_dbp_unadj,
  dbp_adj   = emm_dbp_adj
)
write_rds(fitted_models,  file.path(out_dir, "fitted_models.rds"))
write_rds(emmeans_results, file.path(out_dir, "emmeans_results.rds"))

# --- Emit standardized model sizes for Table A (numbers read from the fits) ---
model_sizes <- imap_dfr(fitted_models, function(fit, key) {
  tibble(
    section       = "5.1",
    role          = "primary",
    analysis      = "Primary BP-trajectory LMM",
    cohort        = "hypertensive",
    outcome       = if (startsWith(key, "sbp")) "SBP" else "DBP",
    predictor_set = if (grepl("unadj", key)) "unadjusted" else "adjusted",
    n_users       = as.integer(lme4::ngrps(fit)),
    n_obs         = as.integer(nobs(fit)),
    n_events      = NA_integer_,
    ceiling_n     = nrow(analysis_users)
  )
})
emit_sizes(model_sizes, "5_1", results_dir)

# --- Structured CSV exports for downstream reference (no re-run needed) ---

# 1. ANOVA tables
anova_table <- imap_dfr(fitted_models, function(fit, name) {
  a <- suppressMessages(anova(fit, type = 3))
  as_tibble(a, rownames = "term") %>% mutate(model = name) %>% select(model, everything())
})
write_csv(anova_table, file.path(results_dir, "5_1_anova_tables.csv"))

# 2. Nakagawa R²
r2_table <- imap_dfr(fitted_models, function(fit, name) {
  r2 <- suppressMessages(r2_nakagawa(fit))
  tibble(model = name,
         marginal_r2    = as.numeric(r2$R2_marginal),
         conditional_r2 = as.numeric(r2$R2_conditional))
})
write_csv(r2_table, file.path(results_dir, "5_1_r2.csv"))

# 3. emmeans (Time × Stratum cell means)
emmeans_table <- imap_dfr(emmeans_results, function(r, name) {
  r$emmeans %>% mutate(model = name) %>% select(model, everything())
})
write_csv(emmeans_table, file.path(results_dir, "5_1_emmeans.csv"))

# 4. Within-stratum Δ contrasts (3mo - BASELINE, 6mo - BASELINE)
delta_table <- imap_dfr(emmeans_results, function(r, name) {
  r$within_stratum %>% mutate(model = name) %>% select(model, everything())
})
write_csv(delta_table, file.path(results_dir, "5_1_delta_contrasts.csv"))

# 5. Cross-stratum Δ-difference contrasts
cross_table <- imap_dfr(emmeans_results, function(r, name) {
  r$cross_stratum %>% mutate(model = name) %>% select(model, everything())
})
write_csv(cross_table, file.path(results_dir, "5_1_cross_stratum_diffs.csv"))

# 6. Within-stratum CONSECUTIVE Δ contrasts (adjacent timepoints: MO3-BASELINE, MO6-MO3)
#    Added for figure significance indicators + later verification of plateau/continued-change.
consec_table <- imap_dfr(emmeans_results, function(r, name) {
  r$consecutive %>% mutate(model = name) %>% select(model, everything())
})
write_csv(consec_table, file.path(results_dir, "5_1_consec_contrasts.csv"))

# 7. Multiplicity diagnostic: across ALL three contrast families (within-stratum,
#    cross-stratum, consecutive), how many contrasts does within-family Holm push
#    from significant→non-sig at alpha = 0.005, vs uncorrected? Lists any
#    "nudged out".
ALPHA <- 0.005
family_rows <- function(df) {
  df <- as.data.frame(df)
  names(df) <- make.unique(names(df))
  idx <- !vapply(df, is.numeric, logical(1))
  lab <- if (any(idx)) do.call(paste, c(df[idx], sep = " | ")) else as.character(seq_len(nrow(df)))
  tibble(contrast = lab, p.value = df$p.value, p.holm = df$p.holm)
}
holm_nudge <- imap_dfr(emmeans_results, function(r, model) {
  bind_rows(
    family_rows(r$within_stratum) %>% mutate(family = "within_stratum"),
    family_rows(r$cross_stratum)  %>% mutate(family = "cross_stratum"),
    family_rows(r$consecutive)    %>% mutate(family = "consecutive")
  ) %>% mutate(model = model, .before = 1)
}) %>%
  mutate(sig_uncorrected = p.value < ALPHA,
         sig_holm        = p.holm  < ALPHA,
         nudged_out      = sig_uncorrected & !sig_holm)
write_csv(holm_nudge, file.path(results_dir, "5_1_holm_nudge.csv"))

cat("\n\n--- Holm multiplicity diagnostic (all families, within-family; alpha =", ALPHA, ") ---\n")
cat(sprintf("  contrasts total:             %d\n", nrow(holm_nudge)))
cat(sprintf("  significant uncorrected:     %d\n", sum(holm_nudge$sig_uncorrected)))
cat(sprintf("  significant after Holm:      %d\n", sum(holm_nudge$sig_holm)))
cat(sprintf("  NUDGED out by Holm:          %d\n", sum(holm_nudge$nudged_out)))
if (any(holm_nudge$nudged_out)) {
  cat("  Nudged contrasts:\n")
  print(holm_nudge %>% filter(nudged_out) %>% select(model, family, contrast, p.value, p.holm))
} else {
  cat("  (none — no contrast crosses alpha under Holm)\n")
}

cat("\n\n=== §5.1 primary LMMs done ===\n")
cat("Saved rds:\n")
cat("  ", file.path(out_dir, "fitted_models.rds"), "\n")
cat("  ", file.path(out_dir, "emmeans_results.rds"), "\n")
cat("Saved CSVs to:\n  ", results_dir, "\n")
for (f in list.files(results_dir, pattern = "\\.csv$")) cat("    ", f, "\n")
cat("Saved log:\n  ", log_path, "\n")
cat("\nNext: 5_1_diagnostics.R (Q-Q, fitted/residual, BLUP Q-Q plots);\n")
cat("      5_1_table_2.R (publication-format Δ + Cohen's dz);\n")
cat("      5_1_figure_1a.R (trajectory plot with 95% CI ribbons).\n")

# Close sink → flush log
sink()
