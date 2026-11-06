# 6_2_logistic.R
# SAP §6.2 (v0.71+): logistic regression of BP-responder on 10 WHOOP-signal predictors,
# fit as DUAL-PRIMARY across SBP-responder AND DBP-responder outcomes.
#
# Per v0.71, DBP-responder is promoted from sensitivity to co-primary (parallel to
# §6.1 v0.53, §7.1 v0.64+v0.65, §8 v0.54). All §6.2 fits run for both outcomes:
# 2 outcomes × 3 cohorts (Primary, FullCohort, HTN-naive) = 6 logistic fits.
#
# Outcomes: binary SBP-responder (ΔSBP ≤ −5 mmHg) and binary DBP-responder
#           (ΔDBP ≤ −5 mmHg), both at 6 months.
# Cohort:   §2.2 hypertensive (Elevated + Stage 1 + Stage 2) + sex restriction +
#           ≥ 8 of 30 days of WHOOP coverage in the baseline window (per v0.62).
# Predictors (10): baseline-window mean of {sleep duration, sleep consistency,
#                  strain, HRV, RHR}; same five metrics as study-window mean.
# Covariates: age, sex, baseline_BMI, baseline_HTN_status (main-effect covariate
#             added per SAP v0.83), plus outcome-matched baseline_BP:
#             baseline_SBP (SBP-responder fits) or baseline_DBP (DBP-responder
#             fits). HTN_BASELINE_ON omitted in the HTN-naive sensitivity
#             (constant — all never_logged).
#
# Per SAP §9: α = 0.005 formal-test threshold; 95% CIs reported.
# Per SAP §6.2 diagnostics (v0.58+): pairwise correlation matrix among the 10
# predictors; calibration plot; AUC. (Cook's distance removed per v0.58 — the
# 4/n threshold is uninformative at n > 7,500.)
#
# Sensitivities (each applied to BOTH primary outcomes per v0.71):
#   - Full-cohort (4 strata, including Normal)
#   - HTN-medication-naive subset (never_logged)
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/6_2/derived/):
#   results/6_2_coefficients_primary.csv          — primary fits (both outcomes) with cohort + outcome columns
#   results/6_2_sensitivity_full_cohort.csv       — full-cohort sensitivity × both outcomes
#   results/6_2_sensitivity_htn_naive.csv         — HTN-naive sensitivity × both outcomes
#   results/6_2_correlation_matrix.csv            — 10x10 predictor correlations (Primary cohort)
#   results/6_2_model_summaries.csv               — AUC, deviance, n, AIC across all 6 fits
#   results/6_2_calibration_data.csv              — calibration deciles (both primary outcomes)
#   results/6_2_vif_primary.csv                   — VIF on the primary SBP-responder fit
#   results/6_2_paper_table5.md                   — paper-ready Table 5 (dual-outcome layout)
#   results/6_2_logistic.log                      — sink log (Caveat: when sourced
#                                                   via `source(..., echo = TRUE)` the
#                                                   on.exit + sink interaction captures
#                                                   only the initial echo line; analytic
#                                                   outputs above are produced correctly
#                                                   regardless. To capture the full
#                                                   transcript, run via `Rscript` or
#                                                   capture R console output separately.
#                                                   Same caveat as §5.1 empirical-median.)
#   figures/6_2_forest_primary.png                — forest plot of 10 ORs, both outcomes overlaid
#   figures/6_2_calibration_primary.png           — calibration plots (both outcomes)
#
# Retired per v0.71:
#   results/6_2_sensitivity_dbp_responder.csv     — DBP-responder is no longer a sensitivity;
#                                                   it is co-primary, so its coefficients now
#                                                   live in 6_2_coefficients_primary.csv. The
#                                                   old single-outcome sensitivity CSV is moved
#                                                   to ~/.Trash/ on first run after v0.71.

library(tidyverse)
library(broom)
library(car)        # for vif()
library(pROC)       # for AUC
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")  # htn_baseline_on(), emit_sizes()

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "6_2", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "6_2_logistic.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)
cat("§6.2 logistic — dual-primary (SBP-responder + DBP-responder) ~ 10 WHOOP signals —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read frame, derive DBP-responder indicator inline (parallel to §6.1 / §7.1) ---
user_frame <- read_rds(input_path) %>%
  mutate(DBP_RESPONDER_6MO = as.integer((MO6_DBP - BASELINE_DBP) <= -5),
         HTN_BASELINE_ON   = htn_baseline_on(HTN_MED_BASELINE))
cat("User frame:", nrow(user_frame), "rows ×", ncol(user_frame), "cols\n")

# --- Predictor + covariate definitions ---
baseline_preds <- c(
  "SLEEP_DURATION_BASELINE_MEAN",
  "SLEEP_CONSISTENCY_BASELINE_MEAN",
  "STRAIN_BASELINE_MEAN",
  "HRV_BASELINE_MEAN",
  "RHR_BASELINE_MEAN"
)
study_preds <- c(
  "SLEEP_DURATION_STUDY_MEAN",
  "SLEEP_CONSISTENCY_STUDY_MEAN",
  "STRAIN_STUDY_MEAN",
  "HRV_STUDY_MEAN",
  "RHR_STUDY_MEAN"
)
all_preds <- c(baseline_preds, study_preds)
covariates_primary   <- c("AGE_AT_START", "SEX", "BASELINE_BMI", "HTN_BASELINE_ON")  # baseline-HTN added per v0.83
covariates_htn_naive <- c("AGE_AT_START", "SEX", "BASELINE_BMI")  # HTN_BASELINE_ON constant in never_logged subset

# --- Cohort definitions ---
build_primary_cohort <- function(d) {
  d %>%
    filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
           SEX %in% c("male", "female"),
           BASELINE_N_DAYS >= 8)
}
build_full_cohort <- function(d) {
  d %>%
    filter(BASELINE_AHA %in% c("Normal", "Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
           SEX %in% c("male", "female"),
           BASELINE_N_DAYS >= 8)
}
build_htn_naive_cohort <- function(d) {
  build_primary_cohort(d) %>%
    filter(HTN_MED_BASELINE == "never_logged")
}

primary_cohort   <- build_primary_cohort(user_frame)
full_cohort      <- build_full_cohort(user_frame)
htn_naive_cohort <- build_htn_naive_cohort(user_frame)

cat("\nCohort sizes (before complete-case filter on predictors):\n")
cat("  Primary (§2.2 + ≥8 baseline days):", nrow(primary_cohort), "\n")
cat("  Full cohort (4 strata):           ", nrow(full_cohort), "\n")
cat("  HTN-naive (primary + never_logged):", nrow(htn_naive_cohort), "\n")

# --- Helper: fit logistic, return fit + AUC + n + data ---
fit_logistic <- function(outcome_col, baseline_covar, data, covariates = covariates_primary, label = "") {
  vars_required <- c(all_preds, covariates, baseline_covar, outcome_col)
  d <- data %>% filter(if_all(all_of(vars_required), ~ !is.na(.)))
  cat(sprintf("\n=== %s ===\n", label))
  cat("Complete-case n:", nrow(d), "(of", nrow(data), "in cohort)\n")
  rhs <- paste(c(baseline_covar, all_preds, covariates), collapse = " + ")
  f <- as.formula(paste0(outcome_col, " ~ ", rhs))
  fit <- glm(f, family = binomial(), data = d)
  d$pred_prob <- predict(fit, type = "response")
  auc_obj <- pROC::roc(d[[outcome_col]], d$pred_prob, quiet = TRUE)
  auc_val <- as.numeric(auc_obj$auc)
  cat("AUC:", round(auc_val, 4), "\n")
  list(fit = fit, n = nrow(d), auc = auc_val, data = d)
}

# =============================================================================
# Outer loop: cohorts × outcomes = 6 fits (2 primary + 4 sensitivity cells)
# =============================================================================
outcome_specs <- list(
  SBP_responder = list(col = "SBP_RESPONDER_6MO", baseline = "BASELINE_SBP"),
  DBP_responder = list(col = "DBP_RESPONDER_6MO", baseline = "BASELINE_DBP")
)
cohort_specs <- list(
  Primary    = list(data = primary_cohort,   covariates = covariates_primary),
  FullCohort = list(data = full_cohort,      covariates = covariates_primary),
  HTN_naive  = list(data = htn_naive_cohort, covariates = covariates_htn_naive)
)

fits     <- list()
fit_meta <- tibble(cohort = character(), outcome = character(), label = character())

for (cn in names(cohort_specs)) {
  for (on in names(outcome_specs)) {
    label <- paste0(cn, "__", on)
    fits[[label]] <- fit_logistic(
      outcome_col    = outcome_specs[[on]]$col,
      baseline_covar = outcome_specs[[on]]$baseline,
      data           = cohort_specs[[cn]]$data,
      covariates     = cohort_specs[[cn]]$covariates,
      label          = paste0(cn, " × ", on)
    )
    fit_meta <- bind_rows(
      fit_meta,
      tibble(cohort = cn, outcome = on, label = label)
    )
  }
}

# --- Extract coefficients (OR + 95% CI + p) with cohort + outcome columns ---
extract_or <- function(fit_res, cohort_label, outcome_label) {
  broom::tidy(fit_res$fit, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE) %>%
    mutate(cohort  = cohort_label,
           outcome = outcome_label,
           n       = fit_res$n,
           auc     = fit_res$auc,
           sig_005 = p.value < 0.005) %>%
    rename(OR = estimate, OR_lcl = conf.low, OR_ucl = conf.high) %>%
    select(cohort, outcome, term, OR, std.error, OR_lcl, OR_ucl,
           statistic, p.value, sig_005, n, auc)
}

all_coefs <- pmap_dfr(fit_meta, function(cohort, outcome, label) {
  extract_or(fits[[label]], cohort, outcome)
})

coef_primary   <- all_coefs %>% filter(cohort == "Primary")
coef_fullcoh   <- all_coefs %>% filter(cohort == "FullCohort")
coef_htn_naive <- all_coefs %>% filter(cohort == "HTN_naive")

write_csv(coef_primary,   file.path(results_dir, "6_2_coefficients_primary.csv"))
write_csv(coef_fullcoh,   file.path(results_dir, "6_2_sensitivity_full_cohort.csv"))
write_csv(coef_htn_naive, file.path(results_dir, "6_2_sensitivity_htn_naive.csv"))

# --- Emit standardized sizes for Table A (per SAP review v0.83) ---
.cohort_label_map <- c(Primary = "hypertensive", FullCohort = "full", HTN_naive = "htn_naive")
.n62_ceiling <- c(Primary = nrow(primary_cohort), FullCohort = nrow(full_cohort),
                  HTN_naive = nrow(htn_naive_cohort))   # each cohort's own entry N
sizes_62 <- pmap_dfr(fit_meta, function(cohort, outcome, label) {
  fr <- fits[[label]]
  ceil <- .n62_ceiling[[cohort]]   # before tibble(): its `cohort` column would shadow the argument
  tibble(
    section       = "6.2",
    role          = if (cohort == "Primary") "primary" else "sensitivity",
    analysis      = "WHOOP-signal logistic (dual-responder)",
    cohort        = .cohort_label_map[[cohort]],
    outcome       = gsub("_", " ", outcome),
    predictor_set = "Primary",
    n_users       = as.integer(fr$n),
    n_obs         = as.integer(fr$n),
    n_events      = as.integer(sum(model.response(model.frame(fr$fit)))),
    ceiling_n     = ceil
  )
})
emit_sizes(sizes_62, "6_2", results_dir)

# Retire the old single-outcome DBP-responder sensitivity CSV — per v0.71 DBP-responder
# is co-primary, so its coefficients live in the primary CSV alongside SBP-responder.
old_dbp_sens_csv <- file.path(results_dir, "6_2_sensitivity_dbp_responder.csv")
if (file.exists(old_dbp_sens_csv)) {
  trash_dir   <- path.expand("~/.Trash")
  trash_target <- file.path(trash_dir,
                            paste0("6_2_sensitivity_dbp_responder_retired_v0.71_",
                                   format(Sys.time(), "%Y-%m-%d"), ".csv"))
  rename_ok <- file.rename(old_dbp_sens_csv, trash_target)
  if (isTRUE(rename_ok)) {
    cat("\nRetired legacy ", basename(old_dbp_sens_csv),
        " to ~/.Trash/ (DBP-responder co-primary per v0.71).\n", sep = "")
  } else {
    warning("Could not retire legacy ", basename(old_dbp_sens_csv),
            " to ~/.Trash/ — file.rename returned FALSE. ",
            "Inspect manually: ", old_dbp_sens_csv)
  }
}

cat("\n\n--- Primary fits: OR + 95% CI + p for the 10 WHOOP-signal predictors ---\n")
for (on in names(outcome_specs)) {
  cat("\n### ", on, " (Primary cohort) ###\n", sep = "")
  print(coef_primary %>%
          filter(outcome == on, term %in% all_preds) %>%
          select(term, OR, OR_lcl, OR_ucl, p.value, sig_005),
        n = Inf, width = Inf)
}

# --- Diagnostic: pairwise correlation matrix among the 10 predictors (Primary SBP fit) ---
primary_sbp_fit <- fits[["Primary__SBP_responder"]]
cor_data <- primary_sbp_fit$data %>% select(all_of(all_preds))
cor_mat <- cor(cor_data, use = "pairwise.complete.obs")
cor_long <- as_tibble(cor_mat, rownames = "predictor1") %>%
  pivot_longer(-predictor1, names_to = "predictor2", values_to = "r")
write_csv(cor_long, file.path(results_dir, "6_2_correlation_matrix.csv"))

cat("\n--- Within-metric baseline-vs-study correlations (key flag) ---\n")
print(cor_long %>%
        filter(predictor1 != predictor2) %>%
        mutate(metric1 = sub("_(BASELINE|STUDY)_MEAN$", "", predictor1),
               metric2 = sub("_(BASELINE|STUDY)_MEAN$", "", predictor2)) %>%
        filter(metric1 == metric2) %>%
        select(predictor1, predictor2, r) %>%
        distinct(), n = Inf)

# --- VIF check (primary SBP-responder fit; same predictor structure for DBP) ---
vif_primary <- tryCatch(car::vif(primary_sbp_fit$fit), error = function(e) NULL)
if (!is.null(vif_primary)) {
  if (is.matrix(vif_primary)) vif_primary <- vif_primary[, "GVIF"]
  vif_table_df <- tibble(predictor = names(vif_primary), VIF = as.numeric(vif_primary))
  write_csv(vif_table_df, file.path(results_dir, "6_2_vif_primary.csv"))
  cat("\n--- VIF check (primary SBP-responder fit; structure identical for DBP-responder) ---\n")
  print(vif_table_df)
  vif_flags <- vif_table_df %>% filter(VIF > 5)
  if (nrow(vif_flags) > 0) {
    cat("\nVIF flags (> 5):\n")
    print(vif_flags)
  } else {
    cat("\nNo predictors with VIF > 5.\n")
  }
}

# --- Model summaries across all 6 fits ---
model_summ <- pmap_dfr(fit_meta, function(cohort, outcome, label) {
  fr <- fits[[label]]
  tibble(cohort   = cohort,
         outcome  = outcome,
         n        = fr$n,
         auc      = fr$auc,
         aic      = AIC(fr$fit),
         null_dev = fr$fit$null.deviance,
         res_dev  = fr$fit$deviance,
         pseudo_r2 = 1 - fr$fit$deviance / fr$fit$null.deviance)
})
write_csv(model_summ, file.path(results_dir, "6_2_model_summaries.csv"))
cat("\n--- Model summaries (all 6 fits) ---\n")
print(model_summ)

# --- Calibration data: both primary outcomes (deciles of predicted probability vs observed rate) ---
calib_one <- function(fit_res, outcome_col) {
  fit_res$data %>%
    mutate(pred_decile = ntile(pred_prob, 10)) %>%
    group_by(pred_decile) %>%
    summarise(
      n           = n(),
      mean_pred   = mean(pred_prob),
      mean_obs    = mean(.data[[outcome_col]]),
      .groups = "drop"
    )
}

calib_sbp <- calib_one(fits[["Primary__SBP_responder"]], "SBP_RESPONDER_6MO") %>%
  mutate(outcome = "SBP_responder")
calib_dbp <- calib_one(fits[["Primary__DBP_responder"]], "DBP_RESPONDER_6MO") %>%
  mutate(outcome = "DBP_responder")
calib_data <- bind_rows(calib_sbp, calib_dbp) %>%
  select(outcome, pred_decile, n, mean_pred, mean_obs)
write_csv(calib_data, file.path(results_dir, "6_2_calibration_data.csv"))

# --- Figures ---
# Forest plot: 10 wearable-biometric predictor ORs, both outcomes overlaid.
# Predictors clustered Baseline-window first, then 6-month window, each cluster
# ordered RHR → HRV → Sleep Duration → Sleep Consistency → Physical Activity.
# legendry::guide_axis_nested draws the window grouping as a bracket OUTSIDE the
# axis (left of the metric tick labels). Axis value = "metric&window".
suppressPackageStartupMessages(library(legendry))

metric_labels <- c(RHR = "RHR", HRV = "HRV",
                   SLEEP_DURATION = "Sleep Duration",
                   SLEEP_CONSISTENCY = "Sleep Consistency",
                   STRAIN = "Physical Activity")
metric_order  <- c("RHR", "HRV", "SLEEP_DURATION", "SLEEP_CONSISTENCY", "STRAIN")
# Bottom-to-top factor order is rev(top-to-bottom); the bracket key splits on "&".
order_top_to_bottom <- c(paste(metric_labels[metric_order], "Baseline", sep = "&"),
                         " ",   # non-breaking-space spacer row → blank gap (avoids an "NA" tick)
                         paste(metric_labels[metric_order], "6 Month",  sep = "&"))
# Real (non-spacer) rows → used as the regular axis key so the spacer gets no tick.
real_levels <- c(paste(metric_labels[metric_order], "Baseline", sep = "&"),
                 paste(metric_labels[metric_order], "6 Month",  sep = "&"))

forest_df <- coef_primary %>%
  filter(term %in% all_preds) %>%
  mutate(
    metric     = sub("_(BASELINE|STUDY)_MEAN$", "", term),
    window     = ifelse(grepl("BASELINE", term), "Baseline", "6 Month"),
    term_label = factor(paste(metric_labels[metric], window, sep = "&"),
                        levels = rev(order_top_to_bottom)),
    # SBP listed last → position_dodge draws it at the higher (top) y position;
    # legend order is restored to SBP-first via scale_color_manual(breaks=).
    outcome    = factor(outcome, levels = c("DBP_responder", "SBP_responder")),
    sig_label  = ifelse(sig_005, "p < 0.005", "p ≥ 0.005")
  )

outcome_cols <- c("SBP_responder" = "#00538A", "DBP_responder" = "#5B9BD5")

p_forest <- ggplot(forest_df,
                   aes(x = OR, y = term_label,
                       color = outcome, shape = sig_label)) +
  geom_vline(xintercept = 1, color = "gray60", linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(xmin = OR_lcl, xmax = OR_ucl), width = 0.4,
                linewidth = 0.5, orientation = "y",
                position = position_dodge(width = 0.6)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_x_log10(limits = c(0.7, 1.4), breaks = seq(0.7, 1.4, by = 0.1),
                expand = expansion(0)) +   # plot area capped 0.7–1.4, labels flush to the box
  scale_y_discrete(drop = FALSE,
                   guide = guide_axis_nested(
                     key = key_range_auto(sep = "&"),
                     regular_key = key_manual(aesthetic = real_levels,
                                              label = sub("&.*$", "", real_levels)),
                     bracket = "square",
                     levels_text = list(
                       element_text(face = "plain", angle = 0),                  # inner metric tick labels
                       element_text(face = "bold",  angle = 90, hjust = 0.5)))) + # centred Baseline / 6 Month titles
  scale_color_manual(values = outcome_cols, name = NULL,
                     breaks = c("SBP_responder", "DBP_responder"),
                     labels = c(SBP_responder = "SBP-Responder", DBP_responder = "DBP-Responder")) +
  scale_shape_manual(values = c("p < 0.005" = 16, "p ≥ 0.005" = 1), guide = "none") +
  labs(
    x = "Adjusted OR (per Unit Predictor, Log Scale, 95% CI)",
    y = NULL,
    title    = "Wearable Derived Biometric Predictors of Clinically Meaningful Blood Pressure Reductions",
    subtitle = "Adjusted Multivariable Logistic Regression"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position    = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.4),
    axis.title.x       = element_text(face = "bold", margin = margin(t = 10)),
    axis.text.y        = element_text(face = "plain"),
    plot.title         = element_text(face = "bold", size = 11)
  )

suppressWarnings(  # benign key_range_auto split note from the blank spacer row
  ggsave(file.path(figures_dir, "6_2_forest_primary.png"), p_forest,
         width = 9.5, height = 7, dpi = 300))
cat("\nWrote forest plot:", file.path(figures_dir, "6_2_forest_primary.png"), "\n")

# Calibration plot — facet by outcome
p_calib <- ggplot(calib_data, aes(x = mean_pred, y = mean_obs)) +
  geom_abline(slope = 1, intercept = 0, color = "gray60", linetype = "dashed") +
  geom_point(aes(color = outcome), size = 3) +
  geom_line(aes(color = outcome), linewidth = 0.6) +
  facet_wrap(~ outcome, scales = "free") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = outcome_cols, name = NULL,
                     labels = c(SBP_responder = "SBP-Responder", DBP_responder = "DBP-Responder")) +
  labs(
    x = "Mean predicted probability (per decile)",
    y = "Observed proportion of responders",
    title    = "§6.2 calibration — predicted vs observed responder probability (dual-primary)",
    subtitle = "10 deciles of predicted probability; dashed = ideal calibration"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(figures_dir, "6_2_calibration_primary.png"), p_calib,
       width = 11, height = 5.5, dpi = 300)
cat("Wrote calibration plot:", file.path(figures_dir, "6_2_calibration_primary.png"), "\n")

# --- Paper-ready markdown (Table 5) ---
fmt_or_md <- function(or, lcl, ucl, p) {
  sig <- ifelse(p < 0.005, "**", "")
  sprintf("%s%.3f (%.3f, %.3f)%s [p = %s]",
          sig, or, lcl, ucl, sig,
          ifelse(p < 1e-10, formatC(p, format = "e", digits = 1),
                 ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                        sprintf("%.3f", p))))
}

md_lines <- c(
  "# §6.2 Table 5 — WHOOP-signal predictors of BP-responder (dual-primary per v0.71)",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  paste0("Primary cohort: §2.2 + sex restriction + ≥8 baseline-window WHOOP days. ",
         "SBP-responder n = ", format(fits[["Primary__SBP_responder"]]$n, big.mark = ","),
         ", AUC = ", round(fits[["Primary__SBP_responder"]]$auc, 3), "; ",
         "DBP-responder n = ", format(fits[["Primary__DBP_responder"]]$n, big.mark = ","),
         ", AUC = ", round(fits[["Primary__DBP_responder"]]$auc, 3), "."),
  "",
  "Bold = significant at α = 0.005 per SAP §9.",
  ""
)

for (on in names(outcome_specs)) {
  md_lines <- c(md_lines, paste0("## Primary fit — ", on, " outcome"), "")
  for (window in c("Baseline", "Study")) {
    md_lines <- c(md_lines, paste0("### ", window, "-window means"))
    win_preds <- if (window == "Baseline") baseline_preds else study_preds
    for (pred in win_preds) {
      r <- coef_primary %>% filter(outcome == on, term == pred)
      if (nrow(r) == 0) next
      pred_clean <- sub("_(BASELINE|STUDY)_MEAN$", "", pred)
      md_lines <- c(md_lines, sprintf(
        "- **%s**: %s",
        pred_clean, fmt_or_md(r$OR, r$OR_lcl, r$OR_ucl, r$p.value)
      ))
    }
    md_lines <- c(md_lines, "")
  }
}

# Sensitivities comparison — both outcomes × both sensitivities for predictors sig in any primary fit
sig_preds <- coef_primary %>%
  filter(term %in% all_preds, sig_005) %>%
  distinct(term) %>% pull(term)

md_lines <- c(md_lines,
  "## Sensitivity analyses (predictors significant in either primary outcome)",
  ""
)
if (length(sig_preds) > 0) {
  for (pred in sig_preds) {
    pred_clean <- gsub("_", " ", sub("_(BASELINE|STUDY)_MEAN$", "", pred)) %>%
      paste0(" (", ifelse(grepl("BASELINE", pred), "Baseline", "Study"), ")")
    md_lines <- c(md_lines, paste0("**", pred_clean, "**"))
    for (cohort_label in c("Primary", "FullCohort", "HTN_naive")) {
      for (on in c("SBP_responder", "DBP_responder")) {
        r <- all_coefs %>% filter(cohort == cohort_label, outcome == on, term == pred)
        if (nrow(r) == 0) next
        md_lines <- c(md_lines, sprintf(
          "- %s × %s (n = %s): %s",
          cohort_label, on, format(r$n, big.mark = ","),
          fmt_or_md(r$OR, r$OR_lcl, r$OR_ucl, r$p.value)
        ))
      }
    }
    md_lines <- c(md_lines, "")
  }
} else {
  md_lines <- c(md_lines, "_(No predictors reached α = 0.005 in either primary fit; sensitivities reported in main coefficient CSVs.)_", "")
}

md_lines <- c(md_lines,
  "## Notes",
  "",
  "- **Model class**: logistic regression `BP_RESPONDER ~ baseline_BP + 10 WHOOP signals + age + sex + baseline_BMI + baseline_HTN_status`. Per SAP v0.83, `baseline_HTN_status` (on-at-baseline indicator: `stable_on_meds` ∪ `recent_initiator` = on; `never_logged` + `mid_study_initiator` = not-on) added as a main-effect covariate. Per SAP v0.71, both SBP-responder (ΔSBP ≤ −5 mmHg) and DBP-responder (ΔDBP ≤ −5 mmHg) are primary outcomes, fit in parallel on the same predictor set with outcome-matched baseline_BP covariate (baseline_SBP for SBP-responder fits, baseline_DBP for DBP-responder fits). `baseline_HTN_status` is **omitted** in the HTN-naive sensitivity (constant — all `never_logged`).",
  "- **Cohort**: §2.2 hypertensive (Elevated + Stage 1 + Stage 2) + sex restriction + ≥8 days of WHOOP coverage in the 30-day pre-baseline window per SAP §6.2 (and §11). Per-fit complete-case n is reported per cell.",
  "- **Predictors**: 5 daily WHOOP metrics (sleep duration, sleep consistency, strain, HRV, RHR) computed as both 30-day pre-baseline mean and 6-month study-window mean. Per SAP §6.2, both forms retained by design to capture starting-state behavior and study-period behavior independently.",
  "- **Significance**: α = 0.005 per SAP §9.",
  "- **Pairwise correlations**: see `6_2_correlation_matrix.csv`. Within-metric baseline-vs-study correlations are reported. VIF check at `6_2_vif_primary.csv` (computed on the primary SBP-responder fit; structure identical for DBP-responder fits, so VIF values transfer).",
  "- **Sensitivities** (each applied to both primary outcomes per v0.71): full-cohort (4 strata) tests robustness to the hypertensive-only inclusion criterion; HTN-medication-naive subset (`never_logged`) tests whether predictor associations are driven by medication-treated users.",
  "- Pairs with: §6.1 LE8 multivariable (continuous-outcome + dual-responder logistic companion); §5.2 clinical-threshold proportions (definition of responder); §7 daily within/between LMM (within-person responsiveness)."
)

writeLines(md_lines, file.path(results_dir, "6_2_paper_table5.md"))
cat("\nWrote markdown:", file.path(results_dir, "6_2_paper_table5.md"), "\n")

cat("\n=== §6.2 logistic done (dual-primary per v0.71) ===\n")
