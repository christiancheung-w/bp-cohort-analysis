# 8_engagement_logistic.R
# SAP §8: Engagement predictors of clinically meaningful BP response.
#
# Per SAP v0.56:
#   PRIMARY FRAME — 7 separate logistic regressions per outcome.
#     Each model: 1 engagement predictor + age + sex + BMI + outcome-specific
#     baseline_BP. Each fit on its own complete-case cohort (§2.2 restricted
#     to non-missing values on the predictor + covariates).
#     Total primary fits: 7 predictors × 2 outcomes = 14.
#
#   SECONDARY FRAME — 1 joint logistic per outcome.
#     All 7 engagement predictors entered simultaneously on the complete-case
#     subset (where all 7 are non-missing). Reported as a mutually-adjusted
#     sensitivity with VIFs. (Joint cohort is now much larger than the prior
#     ~7.5k: AVG_SLEEP_NEED_PCT replaced SLEEP_PCT_OF_TARGET, lifting the
#     binding sleep-coverage constraint from 38% to ~100%.)
#     Total secondary fits: 2.
#
# Both frames refit on the HTN-medication-naive subset (never_logged) as the
# §8 sensitivity (per SAP §8 + §12.8).
#
# Per SAP v0.54: SBP-responder AND DBP-responder are parallel primary outcomes
# (dual-outcome logistic frame, symmetric with §6.1 v0.53).
#
# Engagement predictors (7):
#   - PCT_DAYS_JOURNALED
#   - LOG_AVG_DAILY_TIME_IN_APP_MIN
#   - PCT_DAYS_ENGAGED_NOTIFICATIONS
#   - PCT_DAYS_USED_COACH
#   - AVG_SLEEP_NEED_PCT            (mean % of nightly sleep-need fulfilled; ~100% §2.2 coverage)
#   - PCT_DAYS_MET_OPTIMAL_STRAIN
#   - LOG_N_BP_READINGS_STUDY
#
# Per SAP §9: α = 0.005 formal-test threshold; 95% CIs retained as effect-
# size convention.
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/8_engagement/derived/):
#   results/8_coefficients_separate_primary.csv   — per-predictor OR + 95% CI + p (primary, 14 rows)
#   results/8_coefficients_separate_htn_naive.csv — per-predictor OR + 95% CI + p (HTN-naive, 14 rows)
#   results/8_coefficients_joint_primary.csv      — joint-model OR + 95% CI + p (secondary, 14 rows)
#   results/8_coefficients_joint_htn_naive.csv    — joint-model OR + 95% CI + p (HTN-naive, 14 rows)
#   results/8_cohort_sizes.csv                    — per-predictor and joint cohort sizes
#   results/8_correlation_matrix.csv              — 7x7 pairwise correlations among predictors
#   results/8_vif_joint.csv                       — VIFs from joint secondary fits
#   results/8_model_summaries.csv                 — AUC, AIC, pseudo-R², n per fit
#   results/8_calibration_data.csv                — calibration deciles for all 32 fits (14 separate primary + 2 joint primary + 14 separate HTN-naive + 2 joint HTN-naive) × 10 deciles = 320 rows. Per v0.65 closure (was previously joint-SBP-only).
#   results/8_paper.md                            — paper-ready text
#   results/8.log                                 — sink log (Caveat: source(echo=TRUE) + on.exit + sink interaction captures only initial echo line; analytic outputs unaffected. Same caveat as §5.1 / §6.2 / §7.1 / §7.2.)
#   (violin PNGs 8_figure_5_violin_{sbp,dbp}.png and forest PNGs
#    8_forest_separate_primary.png / 8_forest_joint_secondary.png are rendered by
#    the standalone 8_figure_5.R and 8_figure_forest.R, not this script)
#   figures/8_calibration_separate_sbp.png        — calibration plot, 7 separate-model SBP-responder fits
#   figures/8_calibration_separate_dbp.png        — calibration plot, 7 separate-model DBP-responder fits
#   figures/8_calibration_joint.png               — calibration plot, joint models (2 panels: SBP-responder + DBP-responder)

library(tidyverse)
library(broom)
library(car)        # for vif()
library(pROC)       # for AUC
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "8_engagement", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(out_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "8.log")
sink(log_path, split = TRUE)
cat("§8 Engagement logistic — separate primary + joint secondary —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read frame, build cohort, construct DBP_RESPONDER_6MO ---
user_frame <- read_rds(input_path)
cohort <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female")) %>%
  mutate(DBP_RESPONDER_6MO = as.integer((MO6_DBP - BASELINE_DBP) <= -5),
         SEX = factor(SEX),
         # Per-10-pp rescaling for percentage predictors (SAP_CHANGELOG v0.61).
         # ORs reported per +10 pp rather than per +1 pp for interpretability —
         # the per-1-pp ORs (1.001–1.010) were uninformatively close to unity.
         # Raw 0–100 PCT values retained for the Figure 5 violin plot.
         PCT_DAYS_JOURNALED_per10             = PCT_DAYS_JOURNALED / 10,
         PCT_DAYS_ENGAGED_NOTIFICATIONS_per10 = PCT_DAYS_ENGAGED_NOTIFICATIONS / 10,
         PCT_DAYS_USED_COACH_per10            = PCT_DAYS_USED_COACH / 10,
         AVG_SLEEP_NEED_PCT_per10             = AVG_SLEEP_NEED_PCT / 10,
         PCT_DAYS_MET_OPTIMAL_STRAIN_per10    = PCT_DAYS_MET_OPTIMAL_STRAIN / 10) %>%
  mutate(HTN_BASELINE_ON = htn_baseline_on(HTN_MED_BASELINE))
cat("§2.2 cohort + sex restriction:", nrow(cohort), "users\n")

# --- Predictor + covariate definitions ---
engagement_preds <- c(
  "PCT_DAYS_JOURNALED_per10",
  "LOG_AVG_DAILY_TIME_IN_APP_MIN",
  "PCT_DAYS_ENGAGED_NOTIFICATIONS_per10",
  "PCT_DAYS_USED_COACH_per10",
  "AVG_SLEEP_NEED_PCT_per10",
  "PCT_DAYS_MET_OPTIMAL_STRAIN_per10",
  "LOG_N_BP_READINGS_STUDY"
)
covariates_primary  <- c("AGE_AT_START", "SEX", "BASELINE_BMI", "HTN_BASELINE_ON")
covariates_htn_naive <- c("AGE_AT_START", "SEX", "BASELINE_BMI")

# Pretty labels for predictors
pred_label_map <- c(
  "PCT_DAYS_JOURNALED_per10"             = "% days journaled (per +10 pp)",
  "LOG_AVG_DAILY_TIME_IN_APP_MIN"        = "log(avg daily app min + 1)",
  "PCT_DAYS_ENGAGED_NOTIFICATIONS_per10" = "% days engaged w/ notifications (per +10 pp)",
  "PCT_DAYS_USED_COACH_per10"            = "% days used AI Coach (per +10 pp)",
  "AVG_SLEEP_NEED_PCT_per10"             = "Average sleep need achieved (per +10 pp)",
  "PCT_DAYS_MET_OPTIMAL_STRAIN_per10"    = "% days met optimal strain (per +10 pp)",
  "LOG_N_BP_READINGS_STUDY"              = "log(N BP readings + 1)"
)

# --- Cohort-builder helper (complete-case on given predictor set) ---
build_cohort <- function(d, preds, outcome_col, baseline_covar) {
  required <- c(preds, covariates_primary, baseline_covar, outcome_col)
  d %>% filter(if_all(all_of(required), ~ !is.na(.)))
}

# --- Cohort sizes: per-predictor (separate) and joint ---
sep_sizes <- crossing(
  predictor = engagement_preds,
  outcome   = c("SBP", "DBP")
) %>%
  rowwise() %>%
  mutate(
    baseline_covar = paste0("BASELINE_", outcome),
    outcome_col    = paste0(outcome, "_RESPONDER_6MO"),
    n_complete     = nrow(build_cohort(cohort, predictor, outcome_col, baseline_covar))
  ) %>%
  ungroup() %>%
  mutate(frame = "separate") %>%
  select(frame, predictor, outcome, n_complete)

joint_sizes <- crossing(
  outcome = c("SBP", "DBP")
) %>%
  rowwise() %>%
  mutate(
    baseline_covar = paste0("BASELINE_", outcome),
    outcome_col    = paste0(outcome, "_RESPONDER_6MO"),
    n_complete     = nrow(build_cohort(cohort, engagement_preds, outcome_col, baseline_covar))
  ) %>%
  ungroup() %>%
  mutate(frame = "joint", predictor = "(all 7)") %>%
  select(frame, predictor, outcome, n_complete)

cohort_sizes <- bind_rows(sep_sizes, joint_sizes)
write_csv(cohort_sizes, file.path(results_dir, "8_cohort_sizes.csv"))

cat("\n--- Cohort sizes by frame × predictor × outcome ---\n")
print(cohort_sizes, n = Inf)

# --- Pairwise correlation matrix (on joint complete-case cohort, SBP outcome) ---
cohort_for_cor <- build_cohort(cohort, engagement_preds, "SBP_RESPONDER_6MO", "BASELINE_SBP")
cor_data <- cohort_for_cor %>% select(all_of(engagement_preds))
cor_mat <- cor(cor_data, use = "pairwise.complete.obs")
cor_long <- as_tibble(cor_mat, rownames = "predictor1") %>%
  pivot_longer(-predictor1, names_to = "predictor2", values_to = "r")
write_csv(cor_long, file.path(results_dir, "8_correlation_matrix.csv"))

cat("\n--- High pairwise correlations among engagement predictors (|r| > 0.3) ---\n")
print(cor_long %>%
        filter(predictor1 != predictor2, abs(r) > 0.3) %>%
        arrange(desc(abs(r))) %>%
        head(20))

# --- Helper: fit logistic + extract OR + AUC + pseudo-R² ---
fit_logistic <- function(outcome_col, baseline_covar, predictors, data, label,
                         covariates = covariates_primary) {
  cat(sprintf("\n=== %s (n = %d) ===\n", label, nrow(data)))
  rhs <- paste(c(baseline_covar, predictors, covariates), collapse = " + ")
  f <- as.formula(paste0(outcome_col, " ~ ", rhs))
  fit <- glm(f, family = binomial(), data = data)
  cat("  n_obs:", nobs(fit), "\n")
  d <- data
  d$pred_prob <- predict(fit, type = "response")
  auc_obj <- pROC::roc(d[[outcome_col]], d$pred_prob, quiet = TRUE)
  auc_val <- as.numeric(auc_obj$auc)
  cat("  AUC:", round(auc_val, 4), "\n")
  pseudo_r2 <- 1 - fit$deviance / fit$null.deviance
  cat("  Pseudo-R²:", round(pseudo_r2, 4), "\n")
  list(fit = fit, n = nobs(fit), auc = auc_val, pseudo_r2 = pseudo_r2,
       label = label, data_with_preds = d)
}

# --- Fit loop ---
# Key scheme: "<flavor>__<frame>__<predictor-or-ALL>__<outcome>"
#   flavor    ∈ {Primary, HTNNaive}
#   frame     ∈ {sep, joint}
#   predictor = engagement predictor name (for sep), or "ALL" (for joint)
#   outcome   ∈ {SBP, DBP}
fits <- list()

for (otc in c("SBP", "DBP")) {
  outcome_col   <- paste0(otc, "_RESPONDER_6MO")
  baseline_covar <- paste0("BASELINE_", otc)

  # ----- SEPARATE frame: 7 single-predictor fits -----
  for (p in engagement_preds) {
    d_pri <- build_cohort(cohort, p, outcome_col, baseline_covar)
    cat("\n----------------------------------------\n")
    cat(sprintf("PRIMARY sep: %s × %s (n = %d)\n", p, otc, nrow(d_pri)))
    cat("----------------------------------------\n")
    fits[[paste0("Primary__sep__", p, "__", otc)]] <-
      fit_logistic(outcome_col, baseline_covar, p, d_pri,
                   sprintf("sep %s × %s primary", p, otc))

    d_htn <- d_pri %>% filter(HTN_MED_BASELINE == "never_logged")
    cat(sprintf("\nHTN-naive sep: %s × %s (n = %d)\n", p, otc, nrow(d_htn)))
    fits[[paste0("HTNNaive__sep__", p, "__", otc)]] <-
      fit_logistic(outcome_col, baseline_covar, p, d_htn,
                   sprintf("sep %s × %s HTN-naive", p, otc),
                   covariates = covariates_htn_naive)
  }

  # ----- JOINT frame: 1 multivariable fit -----
  d_joint <- build_cohort(cohort, engagement_preds, outcome_col, baseline_covar)
  cat("\n----------------------------------------\n")
  cat(sprintf("SECONDARY joint × %s (n = %d)\n", otc, nrow(d_joint)))
  cat("----------------------------------------\n")
  fits[[paste0("Primary__joint__ALL__", otc)]] <-
    fit_logistic(outcome_col, baseline_covar, engagement_preds, d_joint,
                 sprintf("joint × %s primary", otc))

  d_joint_htn <- d_joint %>% filter(HTN_MED_BASELINE == "never_logged")
  cat(sprintf("\nHTN-naive joint × %s (n = %d)\n", otc, nrow(d_joint_htn)))
  fits[[paste0("HTNNaive__joint__ALL__", otc)]] <-
    fit_logistic(outcome_col, baseline_covar, engagement_preds, d_joint_htn,
                 sprintf("joint × %s HTN-naive", otc),
                 covariates = covariates_htn_naive)
}

# --- VIF check on joint SBP primary fit ---
vif_joint_sbp <- tryCatch(car::vif(fits[["Primary__joint__ALL__SBP"]]$fit),
                          error = function(e) NULL)
vif_joint_dbp <- tryCatch(car::vif(fits[["Primary__joint__ALL__DBP"]]$fit),
                          error = function(e) NULL)
vif_rows <- bind_rows(
  if (!is.null(vif_joint_sbp)) tibble(outcome = "SBP",
                                      predictor = names(vif_joint_sbp),
                                      VIF = as.numeric(vif_joint_sbp)),
  if (!is.null(vif_joint_dbp)) tibble(outcome = "DBP",
                                      predictor = names(vif_joint_dbp),
                                      VIF = as.numeric(vif_joint_dbp))
)
if (nrow(vif_rows) > 0) {
  cat("\n--- VIF check (joint secondary fits) ---\n")
  print(vif_rows, n = Inf)
  vif_flags <- vif_rows %>% filter(VIF > 5)
  if (nrow(vif_flags) > 0) {
    cat("\nVIF flags (> 5):\n"); print(vif_flags)
  } else {
    cat("\nNo predictors with VIF > 5 in joint fits.\n")
  }
  write_csv(vif_rows, file.path(results_dir, "8_vif_joint.csv"))
}

# --- Extract coefficients ---
extract_or <- function(fr, key) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  flavor    <- parts[1]
  frame     <- parts[2]
  predictor <- parts[3]
  outcome   <- parts[4]
  broom::tidy(fr$fit, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE) %>%
    mutate(flavor    = flavor,
           frame     = frame,
           predictor = predictor,
           outcome   = outcome,
           model     = fr$label,
           n         = fr$n,
           auc       = fr$auc,
           sig_005   = p.value < 0.005) %>%
    rename(OR = estimate, OR_lcl = conf.low, OR_ucl = conf.high) %>%
    select(flavor, frame, predictor, outcome, term, OR, std.error, OR_lcl, OR_ucl,
           statistic, p.value, sig_005, n, auc, model)
}

coefs_all <- imap_dfr(fits, extract_or)

# Build separate-model engagement-OR table (one row per predictor: the
# predictor's own coefficient row from its single-predictor model).
coefs_sep <- coefs_all %>%
  filter(frame == "sep", term == predictor)
coefs_sep_primary <- coefs_sep %>% filter(flavor == "Primary")
coefs_sep_htn     <- coefs_sep %>% filter(flavor == "HTNNaive")
write_csv(coefs_sep_primary, file.path(results_dir, "8_coefficients_separate_primary.csv"))
write_csv(coefs_sep_htn,     file.path(results_dir, "8_coefficients_separate_htn_naive.csv"))

# Build joint-model engagement-OR table (7 rows per outcome × flavor: all
# engagement-predictor coefficient rows from the joint fit).
coefs_joint <- coefs_all %>%
  filter(frame == "joint", term %in% engagement_preds)
coefs_joint_primary <- coefs_joint %>% filter(flavor == "Primary")
coefs_joint_htn     <- coefs_joint %>% filter(flavor == "HTNNaive")
write_csv(coefs_joint_primary, file.path(results_dir, "8_coefficients_joint_primary.csv"))
write_csv(coefs_joint_htn,     file.path(results_dir, "8_coefficients_joint_htn_naive.csv"))

cat("\n\n--- Separate-model primary ORs (one per predictor × outcome) ---\n")
print(coefs_sep_primary %>%
        select(predictor, outcome, OR, OR_lcl, OR_ucl, p.value, sig_005, n),
      n = Inf, width = Inf)

cat("\n--- Separate-model HTN-naive ORs ---\n")
print(coefs_sep_htn %>%
        select(predictor, outcome, OR, OR_lcl, OR_ucl, p.value, sig_005, n),
      n = Inf, width = Inf)

cat("\n--- Joint-model primary ORs (7 predictors × 2 outcomes) ---\n")
print(coefs_joint_primary %>%
        select(term, outcome, OR, OR_lcl, OR_ucl, p.value, sig_005, n),
      n = Inf, width = Inf)

cat("\n--- Joint-model HTN-naive ORs ---\n")
print(coefs_joint_htn %>%
        select(term, outcome, OR, OR_lcl, OR_ucl, p.value, sig_005, n),
      n = Inf, width = Inf)

# --- Model summaries ---
model_summ <- imap_dfr(fits, function(fr, key) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  tibble(flavor    = parts[1],
         frame     = parts[2],
         predictor = parts[3],
         outcome   = parts[4],
         n         = fr$n, auc = fr$auc, aic = AIC(fr$fit),
         pseudo_r2 = fr$pseudo_r2)
})
write_csv(model_summ, file.path(results_dir, "8_model_summaries.csv"))

cat("\n--- Model summaries (all fits) ---\n")
print(model_summ, n = Inf)

# --- emit_sizes: §8 ---
ceiling_8 <- as.integer(nrow(cohort))
ceiling_8_htn <- as.integer(sum(cohort$HTN_MED_BASELINE == "never_logged"))   # HTN-naive entry N
get_n_events <- function(key) {
  fr <- fits[[key]]
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  outcome_col <- paste0(parts[4], "_RESPONDER_6MO")
  as.integer(sum(fr$data_with_preds[[outcome_col]], na.rm = TRUE))
}
sizes_8_sep_pri <- crossing(predictor = engagement_preds, outcome = c("SBP", "DBP")) %>%
  rowwise() %>%
  mutate(
    key    = paste0("Primary__sep__", predictor, "__", outcome),
    n_comp = fits[[key]]$n,
    n_ev   = get_n_events(key)
  ) %>%
  ungroup() %>%
  transmute(
    section = "8", role = "primary", analysis = "Engagement separate logistic",
    cohort = "hypertensive", outcome = paste0(outcome, " responder"),
    predictor_set = predictor, n_users = as.integer(n_comp), n_obs = as.integer(n_comp),
    n_events = as.integer(n_ev), ceiling_n = ceiling_8
  )
sizes_8_joint_pri <- tibble(outcome = c("SBP", "DBP")) %>%
  rowwise() %>%
  mutate(
    key    = paste0("Primary__joint__ALL__", outcome),
    n_comp = fits[[key]]$n,
    n_ev   = get_n_events(key)
  ) %>%
  ungroup() %>%
  transmute(
    section = "8", role = "primary", analysis = "Engagement joint logistic",
    cohort = "hypertensive", outcome = paste0(outcome, " responder"),
    predictor_set = "joint", n_users = as.integer(n_comp), n_obs = as.integer(n_comp),
    n_events = as.integer(n_ev), ceiling_n = ceiling_8
  )
sizes_8_sep_htn <- crossing(predictor = engagement_preds, outcome = c("SBP", "DBP")) %>%
  rowwise() %>%
  mutate(
    key    = paste0("HTNNaive__sep__", predictor, "__", outcome),
    n_comp = fits[[key]]$n,
    n_ev   = get_n_events(key)
  ) %>%
  ungroup() %>%
  transmute(
    section = "8", role = "sensitivity", analysis = "Engagement separate logistic",
    cohort = "htn_naive", outcome = paste0(outcome, " responder"),
    predictor_set = predictor, n_users = as.integer(n_comp), n_obs = as.integer(n_comp),
    n_events = as.integer(n_ev), ceiling_n = ceiling_8_htn
  )
sizes_8_joint_htn <- tibble(outcome = c("SBP", "DBP")) %>%
  rowwise() %>%
  mutate(
    key    = paste0("HTNNaive__joint__ALL__", outcome),
    n_comp = fits[[key]]$n,
    n_ev   = get_n_events(key)
  ) %>%
  ungroup() %>%
  transmute(
    section = "8", role = "sensitivity", analysis = "Engagement joint logistic",
    cohort = "htn_naive", outcome = paste0(outcome, " responder"),
    predictor_set = "joint", n_users = as.integer(n_comp), n_obs = as.integer(n_comp),
    n_events = as.integer(n_ev), ceiling_n = ceiling_8_htn
  )
sizes_8 <- bind_rows(sizes_8_sep_pri, sizes_8_joint_pri, sizes_8_sep_htn, sizes_8_joint_htn)
emit_sizes(sizes_8, "8", results_dir)

# --- Calibration data (all 32 fits: Primary + HTNNaive × separate + joint × SBP + DBP) ---
# Per SAP §8 Diagnostic spec — calibration plot (predicted vs observed responder
# probability by decile of predicted probability). Computed for every fit; long-
# format CSV covers all 32 fits, PNGs rendered only for the 16 Primary fits below.
build_calib <- function(fr, key) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  flavor    <- parts[1]
  frame     <- parts[2]
  predictor <- parts[3]
  outcome   <- parts[4]
  outcome_col <- paste0(outcome, "_RESPONDER_6MO")
  fr$data_with_preds %>%
    mutate(pred_decile = ntile(pred_prob, 10)) %>%
    group_by(pred_decile) %>%
    summarise(n         = n(),
              mean_pred = mean(pred_prob),
              mean_obs  = mean(.data[[outcome_col]]),
              .groups   = "drop") %>%
    mutate(flavor    = flavor,
           frame     = frame,
           predictor = predictor,
           outcome   = outcome) %>%
    select(flavor, frame, predictor, outcome, pred_decile, n, mean_pred, mean_obs)
}

calib_data <- imap_dfr(fits, build_calib)
write_csv(calib_data, file.path(results_dir, "8_calibration_data.csv"))

cat("\n--- Calibration computed for all fits (long-format CSV) ---\n")
cat("Rows:", nrow(calib_data), "(expected 32 fits × 10 deciles = 320)\n")

# --- Figure 5 (violins) ---
# The 7-panel violin PNGs (8_figure_5_violin_{sbp,dbp}.png) are rendered by the
# standalone, house-styled 8_figure_5.R — NOT here. Keeping the violin code out of
# this analysis script prevents it from overwriting the polished, tail-clipped
# figures on every re-run. Run 8_figure_5.R after this script.

# --- Forest plots ---
# Forest PNGs (8_forest_separate_primary.png, 8_forest_joint_secondary.png) are
# rendered by the standalone, house-styled 8_figure_forest.R — NOT here. Keeping
# the forest code out of this analysis script prevents it from overwriting the
# polished figures on every re-run. Run 8_figure_forest.R after this script.

# --- Calibration plots: 3 PNGs (separates × SBP, separates × DBP, joint × {SBP, DBP}) ---
# PNGs rendered for the 16 Primary fits only. HTN-naive calibration data is
# available in the CSV but no PNG rendered (out of scope per Tier 4 / Item 2).

# Per-outcome 7-panel grid for primary separates
for (otc in c("SBP", "DBP")) {
  d_sep <- calib_data %>%
    filter(flavor == "Primary", frame == "sep", outcome == otc) %>%
    mutate(predictor_label = factor(pred_label_map[predictor],
                                    levels = pred_label_map))
  p_sep <- ggplot(d_sep, aes(x = mean_pred, y = mean_obs)) +
    geom_abline(slope = 1, intercept = 0, color = "gray60", linetype = "dashed") +
    geom_point(size = 2.2, color = "#1f78b4") +
    geom_line(color = "#1f78b4", linewidth = 0.5) +
    facet_wrap(~ predictor_label, ncol = 4, scales = "free") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Mean predicted probability (per decile)",
      y = "Observed proportion of responders",
      title    = paste0("§8 calibration — ", otc, "-responder primary separates (7 models)"),
      subtitle = "10 deciles per model; dashed = ideal calibration",
      caption  = "Each panel = one of 7 separate per-predictor logistic regressions (primary frame, SAP v0.56)."
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title   = element_text(face = "bold"),
          strip.text   = element_text(face = "bold", size = 9),
          plot.caption = element_text(size = 8, color = "gray40", hjust = 0))
  fname <- paste0("8_calibration_separate_", tolower(otc), ".png")
  out_png <- file.path(figures_dir, fname)
  ggsave(out_png, p_sep, width = 12, height = 7, dpi = 300)
  cat("Wrote calibration (separate ", otc, "):\n  ", out_png, "\n", sep = "")
}

# 2-panel joint (secondary frame) — SBP + DBP side-by-side
d_joint <- calib_data %>%
  filter(flavor == "Primary", frame == "joint") %>%
  mutate(outcome_label = factor(outcome,
                                levels = c("SBP", "DBP"),
                                labels = c("SBP-responder (joint)", "DBP-responder (joint)")))
p_joint <- ggplot(d_joint, aes(x = mean_pred, y = mean_obs)) +
  geom_abline(slope = 1, intercept = 0, color = "gray60", linetype = "dashed") +
  geom_point(size = 3, color = "#1f78b4") +
  geom_line(color = "#1f78b4", linewidth = 0.6) +
  facet_wrap(~ outcome_label, ncol = 2, scales = "free") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Mean predicted probability (per decile)",
    y = "Observed proportion of responders",
    title    = "§8 calibration — joint-model secondary frame",
    subtitle = "All 7 engagement predictors entered jointly; 10 deciles of predicted probability",
    caption  = "Secondary frame per SAP v0.56 — mutually-adjusted multivariable logistic."
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title   = element_text(face = "bold"),
        strip.text   = element_text(face = "bold"),
        plot.caption = element_text(size = 9, color = "gray40", hjust = 0))
out_png <- file.path(figures_dir, "8_calibration_joint.png")
ggsave(out_png, p_joint, width = 12, height = 6, dpi = 300)
cat("Wrote calibration (joint, 2-panel):\n  ", out_png, "\n")

# --- Paper-ready markdown ---
fmt_or <- function(or, lcl, ucl, p) {
  sig <- ifelse(p < 0.005, "**", "")
  sprintf("%s%.3f (%.3f, %.3f)%s [p = %s]",
          sig, or, lcl, ucl, sig,
          ifelse(p < 1e-10, formatC(p, format = "e", digits = 1),
                 ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                        sprintf("%.3f", p))))
}

md_lines <- c(
  "# §8 Engagement predictors of BP-responder status",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "Per SAP v0.56: 7 separate per-predictor logistic models as the primary frame; a joint multivariable model as the secondary mutually-adjusted sensitivity. Per SAP v0.54, SBP-responder and DBP-responder are parallel primary outcomes.",
  "",
  "Bold = significant at α = 0.005 per SAP §9.",
  "",
  "## Cohort sizes",
  ""
)
for (i in seq_len(nrow(cohort_sizes))) {
  r <- cohort_sizes[i, ]
  md_lines <- c(md_lines, sprintf("- **%s** | %s × %s: n = %s",
                                   r$frame, r$predictor, r$outcome,
                                   format(r$n_complete, big.mark = ",")))
}

md_lines <- c(md_lines, "", "## Primary frame — separate logistic models (marginal OR per predictor)", "")
for (otc in c("SBP", "DBP")) {
  md_lines <- c(md_lines, paste0("### ", otc, "-responder"))
  rows <- coefs_sep_primary %>% filter(outcome == otc)
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    md_lines <- c(md_lines, sprintf("- **%s** (n = %s): %s",
                                     pred_label_map[r$predictor],
                                     format(r$n, big.mark = ","),
                                     fmt_or(r$OR, r$OR_lcl, r$OR_ucl, r$p.value)))
  }
  md_lines <- c(md_lines, "")
}

md_lines <- c(md_lines, "## Secondary frame — joint logistic model (mutually-adjusted OR)", "")
for (otc in c("SBP", "DBP")) {
  md_lines <- c(md_lines, paste0("### ", otc, "-responder (joint, n = ",
                                  format(joint_sizes$n_complete[joint_sizes$outcome == otc],
                                         big.mark = ","), ")"))
  rows <- coefs_joint_primary %>% filter(outcome == otc)
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    md_lines <- c(md_lines, sprintf("- **%s**: %s",
                                     pred_label_map[r$term],
                                     fmt_or(r$OR, r$OR_lcl, r$OR_ucl, r$p.value)))
  }
  md_lines <- c(md_lines, "")
}

md_lines <- c(md_lines, "## HTN-medication-naive sensitivity — separate models", "")
for (otc in c("SBP", "DBP")) {
  md_lines <- c(md_lines, paste0("### ", otc, "-responder (HTN-naive)"))
  rows <- coefs_sep_htn %>% filter(outcome == otc)
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    md_lines <- c(md_lines, sprintf("- **%s** (n = %s): %s",
                                     pred_label_map[r$predictor],
                                     format(r$n, big.mark = ","),
                                     fmt_or(r$OR, r$OR_lcl, r$OR_ucl, r$p.value)))
  }
  md_lines <- c(md_lines, "")
}

md_lines <- c(md_lines, "## HTN-medication-naive sensitivity — joint model", "")
for (otc in c("SBP", "DBP")) {
  htn_n <- model_summ %>%
    filter(flavor == "HTNNaive", frame == "joint", outcome == otc) %>%
    pull(n)
  md_lines <- c(md_lines, paste0("### ", otc, "-responder (joint, HTN-naive, n = ",
                                  format(htn_n, big.mark = ","), ")"))
  rows <- coefs_joint_htn %>% filter(outcome == otc)
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    md_lines <- c(md_lines, sprintf("- **%s**: %s",
                                     pred_label_map[r$term],
                                     fmt_or(r$OR, r$OR_lcl, r$OR_ucl, r$p.value)))
  }
  md_lines <- c(md_lines, "")
}

md_lines <- c(md_lines,
  "## Notes",
  "",
  "- **Model class**: logistic regression, fit by maximum likelihood. Per SAP v0.54, both SBP-responder and DBP-responder are primary frames.",
  "- **Primary frame (SAP v0.56)**: 7 separate single-predictor models per outcome — each fit on its own complete-case cohort (engagement predictor + age/sex/BMI/baseline_BP). Marginal interpretation: each OR reflects the association between that engagement metric and response, adjusted only for demographics + baseline BP.",
  "- **Secondary frame**: single joint logistic per outcome with all 7 predictors entered simultaneously on the ~7.5k complete-case subset. Mutually-adjusted interpretation: each OR is conditional on the other 6 engagement predictors. VIFs reported in `8_vif_joint.csv`.",
  "- **Cohort**: §2.2 hypertensive + sex restriction + complete-case on the relevant predictor set.",
  "- **Significance**: α = 0.005 per SAP §9.",
  "- **Sensitivity**: HTN-medication-naive subset (`never_logged` HTN-meds) refit on both frames.",
  "- **Figure 5**: 7-panel violin plots by responder group — parallel SBP-responder (`8_figure_5_violin_sbp.png`) and DBP-responder (`8_figure_5_violin_dbp.png`) versions matching the dual-outcome primary structure.",
  "- **Calibration diagnostics (per SAP §8 spec)**: predicted-vs-observed responder probability by deciles. Long-format CSV at `derived/results/8_calibration_data.csv` covers all 32 fits (Primary + HTN-naive × separate + joint × 2 outcomes). PNGs (Primary frame only): per-outcome 7-panel grids at `derived/figures/8_calibration_separate_{sbp,dbp}.png` for the primary separate frame; 2-panel joint at `8_calibration_joint.png` for the secondary frame.",
  "- Pairs with: §6.1 LE8 predictors; §6.2 WHOOP-signal predictors; §7.1 trajectories; §7.2 W/B daily LMM."
)

writeLines(md_lines, file.path(results_dir, "8_paper.md"))
cat("\nWrote markdown:\n  ", file.path(results_dir, "8_paper.md"), "\n")

cat("\n=== §8 engagement logistic done ===\n")
sink()
