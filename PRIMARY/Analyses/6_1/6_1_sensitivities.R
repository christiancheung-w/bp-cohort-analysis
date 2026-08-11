# 6_1_sensitivities.R
# SAP §6.1 sensitivity analyses: full-cohort (4 strata including Normal) and
# HTN-medication-naive subset replication. Each applied to all four primary
# frames per SAP v0.53:
#   - Continuous ANCOVA: SBP, DBP
#   - Logistic-responder: SBP-responder, DBP-responder
# Across all six SAP §6.1 variants (Primary, +Δweight, +ΔHDL, +ΔNon-HDL,
# +ΔTG, +ΔHbA1c). 4 outcomes × 6 variants × 2 sensitivities = 48 fits total.
#
# Convergence in direction and rough magnitude with the primary §6.1 (§2.2
# cohort) result demonstrates robustness to:
#   - Inclusion of the Normal stratum (full-cohort)
#   - Medication-treated users (HTN-naive)
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/6_1/derived/):
#   results/6_1_sensitivity_coefficients.csv     — long-format: sensitivity × outcome × variant × term → β/OR + CI + p
#   results/6_1_sensitivity_cohort_sizes.csv     — sensitivity × variant → n
#   results/6_1_sensitivity_comparison.csv       — side-by-side vs primary for significant LE8 predictors
#   results/6_1_sensitivities_paper.md           — paper-quote-ready summary
#   results/6_1_sensitivities.log                — full sink log

library(tidyverse)
library(broom)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")  # htn_baseline_on(), emit_sizes()

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "6_1", "derived")
results_dir <- file.path(out_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "6_1_sensitivities.log")
sink(log_path, split = TRUE)
cat("§6.1 sensitivities — full-cohort + HTN-naive —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read frame, construct DBP_RESPONDER_6MO ---
user_frame <- read_rds(input_path) %>%
  mutate(DBP_RESPONDER_6MO = as.integer((MO6_DBP - BASELINE_DBP) <= -5),
         HTN_BASELINE_ON   = htn_baseline_on(HTN_MED_BASELINE))  # baseline-only HTN trait (v0.82)

# --- Predictor + covariate definitions (match primary script) ---
le8_primary  <- c("SLEEP_DURATION_STUDY_MEAN", "STRAIN_STUDY_MEAN",
                  "LOG_TOBACCO_NICOTINE_USE_DAYS", "TOBACCO_NICOTINE_EVER_LOGGED")
le8_dweight  <- c(le8_primary, "DELTA_WEIGHT_KG")
covariates   <- c("AGE_AT_START", "SEX", "BASELINE_BMI",
                  "SLEEP_DURATION_BASELINE_MEAN", "STRAIN_BASELINE_MEAN")  # baseline sleep/PA T0 controls per v0.85

variants <- list(
  Primary     = le8_primary,
  PlusDweight = le8_dweight,
  PlusDHDL    = c(le8_primary, "DELTA_HDL_CHOLESTEROL"),
  PlusDNonHDL = c(le8_primary, "DELTA_NON_HDL_CHOLESTEROL"),
  PlusDTG     = c(le8_primary, "DELTA_TRIGLYCERIDES"),
  PlusDHbA1c  = c(le8_primary, "DELTA_HBA1C")
)
le8_all_terms <- unique(unlist(variants))

# --- Sensitivity cohort builders ---
build_full_cohort <- function(d) {
  d %>%
    filter(BASELINE_AHA %in% c("Normal", "Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
           SEX %in% c("male", "female"))
}
build_htn_naive_cohort <- function(d) {
  d %>%
    filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
           SEX %in% c("male", "female"),
           HTN_MED_BASELINE == "never_logged")
}

sensitivity_cohorts <- list(
  FullCohort = build_full_cohort(user_frame),
  HTN_naive  = build_htn_naive_cohort(user_frame)
)

cat("Sensitivity cohorts (before variant complete-case filter):\n")
for (s in names(sensitivity_cohorts)) {
  cat(sprintf("  %s: %d users\n", s, nrow(sensitivity_cohorts[[s]])))
}

# --- Helpers (identical to primary) ---
# covs passed per-cohort: FullCohort uses baseline-HTN; HTN_naive omits it (constant).
fit_continuous <- function(outcome_col, baseline_col, predictors, data, covs) {
  vars_required <- c(predictors, covs, baseline_col, outcome_col)
  d <- data %>% filter(if_all(all_of(vars_required), ~ !is.na(.)))
  rhs <- paste(c(baseline_col, predictors, covs), collapse = " + ")
  f <- as.formula(paste0(outcome_col, " ~ ", rhs))
  fit <- lm(f, data = d)
  list(fit = fit, n = nrow(d))
}

fit_logistic <- function(outcome_col, baseline_col, predictors, data, covs) {
  vars_required <- c(predictors, covs, baseline_col, outcome_col)
  d <- data %>% filter(if_all(all_of(vars_required), ~ !is.na(.)))
  rhs <- paste(c(baseline_col, predictors, covs), collapse = " + ")
  f <- as.formula(paste0(outcome_col, " ~ ", rhs))
  fit <- glm(f, family = binomial(), data = d)
  list(fit = fit, n = nrow(d))
}

# --- Fit all sensitivities ---
all_coefs <- list()
all_sizes <- list()

for (sname in names(sensitivity_cohorts)) {
  scoh <- sensitivity_cohorts[[sname]]
  # FullCohort adopts baseline-HTN (v0.82); HTN_naive cannot (constant — all never_logged).
  cov_set <- if (sname == "HTN_naive") covariates else c(covariates, "HTN_BASELINE_ON")
  cat("\n\n========================================\n")
  cat("Sensitivity:", sname, " | cohort n =", nrow(scoh),
      " | covs:", paste(cov_set, collapse = ", "), "\n")
  cat("========================================\n")

  for (v in names(variants)) {
    preds <- variants[[v]]
    # Per-variant VIF override: PlusDHbA1c drops SLEEP_DURATION_BASELINE_MEAN
    # (VIF > 5 at n~269; same rule as primary script per SAP v0.86).
    v_cov_set <- if (v == "PlusDHbA1c") setdiff(cov_set, "SLEEP_DURATION_BASELINE_MEAN") else cov_set
    cat(sprintf("\n--- %s × %s ---\n", sname, v))

    # Continuous SBP + DBP
    for (otc in c("SBP", "DBP")) {
      outcome_col  <- paste0("MO6_", otc)
      baseline_col <- paste0("BASELINE_", otc)
      f <- fit_continuous(outcome_col, baseline_col, preds, scoh, v_cov_set)
      cat(sprintf("  %s ANCOVA: n = %d\n", otc, f$n))
      coef_tbl <- broom::tidy(f$fit, conf.int = TRUE, conf.level = 0.95) %>%
        mutate(sensitivity = sname, frame = "continuous", outcome = otc,
               variant = v, n = f$n, sig_005 = p.value < 0.005) %>%
        rename(point = estimate, lcl = conf.low, ucl = conf.high)
      all_coefs[[paste0(sname, "__continuous__", otc, "__", v)]] <- coef_tbl
    }

    # Logistic SBP-responder + DBP-responder
    for (resp in c("SBP", "DBP")) {
      outcome_col  <- paste0(resp, "_RESPONDER_6MO")
      baseline_col <- paste0("BASELINE_", resp)
      f <- fit_logistic(outcome_col, baseline_col, preds, scoh, v_cov_set)
      cat(sprintf("  %s-responder logistic: n = %d\n", resp, f$n))
      coef_tbl <- broom::tidy(f$fit, conf.int = TRUE, conf.level = 0.95,
                              exponentiate = TRUE) %>%
        mutate(sensitivity = sname, frame = "logistic", outcome = paste0(resp, "_responder"),
               variant = v, n = f$n, sig_005 = p.value < 0.005) %>%
        rename(point = estimate, lcl = conf.low, ucl = conf.high)
      all_coefs[[paste0(sname, "__logistic__", resp, "__", v)]] <- coef_tbl
    }

    # Track per-variant n (use the SBP ANCOVA n as canonical)
    n_canonical <- all_coefs[[paste0(sname, "__continuous__SBP__", v)]]$n[1]
    all_sizes[[paste0(sname, "__", v)]] <- tibble(
      sensitivity = sname, variant = v, n = n_canonical
    )
  }
}

# --- Combine + write CSVs ---
coefs_df <- bind_rows(all_coefs) %>%
  select(sensitivity, frame, outcome, variant, term, point, std.error,
         lcl, ucl, statistic, p.value, sig_005, n)
write_csv(coefs_df, file.path(results_dir, "6_1_sensitivity_coefficients.csv"))

sizes_df <- bind_rows(all_sizes)
write_csv(sizes_df, file.path(results_dir, "6_1_sensitivity_cohort_sizes.csv"))

# --- Emit standardized sizes for Table A.2 (per sensitivity × variant × outcome) ---
# n_users = complete-case n (continuous frame; logistic shares it). Responder
# events omitted for the sensitivity robustness checks (carried in the primary).
.sens_ceiling <- vapply(sensitivity_cohorts, nrow, integer(1))
le8_sens_sizes <- coefs_df %>%
  filter(frame == "continuous") %>%
  distinct(sensitivity, variant, outcome, n) %>%
  transmute(
    section = "6.1", role = "sensitivity",
    analysis = paste0("LE8 regression sensitivity (", sensitivity, ")"),
    cohort = sensitivity, outcome = outcome, predictor_set = variant,
    n_users = as.integer(n), n_obs = as.integer(n), n_events = NA_integer_,
    ceiling_n = as.integer(.sens_ceiling[sensitivity])
  )
emit_sizes(le8_sens_sizes, "6_1_sensitivities", results_dir)

cat("\n\n--- Sensitivity cohort sizes per variant ---\n")
print(sizes_df, n = Inf)

# --- Compare to primary §6.1 ---
primary_continuous <- read_csv(file.path(results_dir, "6_1_coefficients_continuous.csv"),
                                show_col_types = FALSE) %>%
  filter(term %in% le8_all_terms) %>%
  transmute(frame = "continuous", outcome = outcome, variant, term,
            primary_point = estimate, primary_p = p.value, primary_sig = sig_005)
primary_logistic <- read_csv(file.path(results_dir, "6_1_coefficients_logistic.csv"),
                              show_col_types = FALSE) %>%
  filter(term %in% le8_all_terms) %>%
  transmute(frame = "logistic", outcome = responder_outcome, variant, term,
            primary_point = OR, primary_p = p.value, primary_sig = sig_005)
primary_all <- bind_rows(primary_continuous, primary_logistic)

comparison <- coefs_df %>%
  filter(term %in% le8_all_terms) %>%
  left_join(primary_all, by = c("frame", "outcome", "variant", "term")) %>%
  mutate(diff_vs_primary = point - primary_point) %>%
  select(sensitivity, frame, outcome, variant, term,
         sens_point = point, sens_lcl = lcl, sens_ucl = ucl, sens_p = p.value,
         primary_point, primary_p, diff_vs_primary,
         sens_sig_005 = sig_005, primary_sig)

write_csv(comparison, file.path(results_dir, "6_1_sensitivity_comparison.csv"))

# --- Filter to LE8 predictors that were significant in PRIMARY (for direct comparison) ---
sig_in_primary <- primary_all %>% filter(primary_sig)
cat("\n--- LE8 predictors significant in PRIMARY (α = 0.005); convergence in sensitivities? ---\n")

comparison_sig <- comparison %>%
  semi_join(sig_in_primary, by = c("frame", "outcome", "variant", "term")) %>%
  mutate(
    direction_concordant = sign(sens_point - ifelse(frame == "logistic", 1, 0)) ==
                           sign(primary_point - ifelse(frame == "logistic", 1, 0)),
    magnitude_ratio = ifelse(frame == "logistic",
                             log(sens_point) / log(primary_point),
                             sens_point / primary_point)
  ) %>%
  arrange(sensitivity, frame, outcome, variant, term)

print(comparison_sig %>%
        select(sensitivity, frame, outcome, variant, term,
               sens_point, primary_point, sens_p, primary_p,
               direction_concordant, sens_sig_005),
      n = Inf, width = Inf)

cat("\n\n--- Summary: convergence with primary across sensitivities ---\n")
convergence <- comparison_sig %>%
  group_by(sensitivity, frame) %>%
  summarise(
    n_cells_tested      = n(),
    n_direction_concord = sum(direction_concordant, na.rm = TRUE),
    n_sig_in_sens       = sum(sens_sig_005),
    .groups = "drop"
  )
print(convergence)

# --- Paper-ready markdown ---
fmt_b <- function(est, lcl, ucl, p, alpha = 0.005) {
  sig <- ifelse(p < alpha, "**", "")
  sprintf("%s%+.3f (%.3f, %.3f)%s [p = %s]",
          sig, est, lcl, ucl, sig,
          ifelse(p < 1e-10, formatC(p, format = "e", digits = 1),
                 ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                        sprintf("%.3f", p))))
}
fmt_or <- function(or, lcl, ucl, p, alpha = 0.005) {
  sig <- ifelse(p < alpha, "**", "")
  sprintf("%s%.3f (%.3f, %.3f)%s [p = %s]",
          sig, or, lcl, ucl, sig,
          ifelse(p < 1e-10, formatC(p, format = "e", digits = 1),
                 ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                        sprintf("%.3f", p))))
}

md_lines <- c(
  "# §6.1 sensitivity analyses — full-cohort + HTN-naive",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "Per SAP v0.53, both sensitivities apply to all four primary frames (continuous SBP, continuous DBP, SBP-responder, DBP-responder). Bold = significant at α = 0.005.",
  "",
  "## Convergence summary"
)
md_lines <- c(md_lines, "",
  paste0("| sensitivity | frame | n primary-sig cells tested | direction-concordant | also sig in sensitivity |"),
  "|---|---|---|---|---|"
)
for (i in seq_len(nrow(convergence))) {
  r <- convergence[i, ]
  md_lines <- c(md_lines, sprintf(
    "| %s | %s | %d | %d/%d | %d/%d |",
    r$sensitivity, r$frame,
    r$n_cells_tested,
    r$n_direction_concord, r$n_cells_tested,
    r$n_sig_in_sens, r$n_cells_tested
  ))
}

md_lines <- c(md_lines, "",
  "## Cell-by-cell comparison for primary-significant predictors",
  ""
)
for (sname in c("FullCohort", "HTN_naive")) {
  md_lines <- c(md_lines, paste0("### Sensitivity: ", sname))
  for (otc in c("SBP", "DBP", "SBP_responder", "DBP_responder")) {
    rows <- comparison_sig %>% filter(sensitivity == sname, outcome == otc)
    if (nrow(rows) == 0) next
    md_lines <- c(md_lines, paste0("#### ", otc))
    for (i in seq_len(nrow(rows))) {
      r <- rows[i, ]
      if (r$frame == "continuous") {
        sens_str <- fmt_b(r$sens_point, r$sens_lcl, r$sens_ucl, r$sens_p)
      } else {
        sens_str <- fmt_or(r$sens_point, r$sens_lcl, r$sens_ucl, r$sens_p)
      }
      md_lines <- c(md_lines, sprintf(
        "- `%s` (%s, variant=%s): primary %s = %s, sensitivity = %s, diff = %+.3f.",
        r$term, r$frame, r$variant,
        ifelse(r$frame == "continuous", "β", "OR"),
        sprintf("%.3f", r$primary_point),
        sens_str, r$diff_vs_primary
      ))
    }
    md_lines <- c(md_lines, "")
  }
}

md_lines <- c(md_lines,
  "## Notes",
  "",
  "- Both sensitivities apply to all four primary frames per SAP v0.53 (continuous SBP, continuous DBP, SBP-responder logistic, DBP-responder logistic).",
  "- **Full-cohort sensitivity** adds Normal-stratum users (~8k) to test robustness to the §2.2 hypertensive-only inclusion criterion.",
  "- **HTN-naive sensitivity** restricts to users with no HTN-medication entries in the 180-day lookback or study window — tests whether predictor-outcome associations are driven by medication-treated users.",
  "- Convergence in direction and rough magnitude across sensitivities = robustness; divergence = compositional dependence.",
  "- Direction-concordance for continuous frame = same sign (β); for logistic frame = OR on same side of 1.",
  "- Pairs with: §6.1 primary (`6_1_le8_regression.R`), SAP v0.53."
)

writeLines(md_lines, file.path(results_dir, "6_1_sensitivities_paper.md"))
cat("\nWrote markdown:\n  ", file.path(results_dir, "6_1_sensitivities_paper.md"), "\n")

cat("\n=== §6.1 sensitivities done ===\n")
sink()
