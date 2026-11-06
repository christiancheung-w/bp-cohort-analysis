# mid_excluded_sensitivity.R
# -----------------------------------------------------------------------------
# Sensitivity analysis excluding participants who first logged antihypertensive
# medication after baseline (HTN_MED_BASELINE == "mid_study_initiator").
#
# Built from mid_excluded_sensitivity.R with the cohort filter swapped; model specs
# are the primary specs unchanged (HTN_BASELINE_ON still varies and is kept).
#   - Scope: §5.1, §6.2, §7.2, §8 (the analyses retained for npj).
#   - Discordance rule: significance flip at alpha = 0.005 OR sign reversal vs
#     the primary cohort.
#   - Fidelity gate: each primary re-run must reproduce the frozen primary CSVs.
#
# Run from repo root:  Rscript PRIMARY/Analyses/mid_excluded_sensitivity/mid_excluded_sensitivity.R
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(performance)
  library(broom)
  library(car)
  library(pROC)
  library(clubSandwich)
})
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")

ALPHA <- 0.005

project_dir  <- "./PRIMARY"
frame_path   <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
reading_path <- file.path(project_dir, "Build", "derived", "qc", "bp_reading_analytic_frame_1X.rds")

out_dir     <- file.path(project_dir, "Analyses", "mid_excluded_sensitivity", "derived")
results_dir <- file.path(out_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Frozen primary result roots (for fidelity gate)
p51 <- file.path(project_dir, "Analyses", "5_1", "derived", "results")

# --- Load frame ---------------------------------------------------------------
user_frame <- read_rds(frame_path)
stopifnot(!any(is.na(user_frame$HTN_MED_BASELINE)))
cat("User frame:", nrow(user_frame), "rows;  mid_study_initiator:",
    sum(user_frame$HTN_MED_BASELINE == "mid_study_initiator"), "\n")

# Cohort builders. `base` applies the §2.2 + sex filter; the two cohorts differ
# ONLY by the mid-study-initiator exclusion.
cohort_variants <- function(base_df) {
  list(primary    = base_df,
       mid_excluded = base_df %>% filter(!(HTN_MED_BASELINE %in% "mid_study_initiator")))
}

# --- Discordance machinery ---------------------------------------------------
# Standardize any result table to (section, family, key, stat, p). `stat` is a
# signed statistic (t/z/estimate) used ONLY for sign; `p` drives significance.
standardize <- function(df, key_cols, stat_col, p_col, section, family) {
  df %>%
    transmute(
      section = section,
      family  = family,
      key     = apply(across(all_of(key_cols)), 1, paste, collapse = " | "),
      stat    = .data[[stat_col]],
      p       = .data[[p_col]]
    )
}

# Compare naive vs primary standardized tables -> discordance rows.
discordance <- function(std_prim, std_naive) {
  full_join(
    std_prim  %>% rename(prim_stat = stat,  prim_p = p),
    std_naive %>% rename(naive_stat = stat, naive_p = p),
    by = c("section", "family", "key")
  ) %>%
    mutate(
      prim_sig  = prim_p  < ALPHA,
      naive_sig = naive_p < ALPHA,
      sig_flip  = (prim_sig != naive_sig) %in% TRUE,
      sign_flip = (sign(prim_stat) != sign(naive_stat)) %in% TRUE,
      discordant = sig_flip | sign_flip
    )
}

disc_all <- list()   # accumulate discordance rows across chapters

# =============================================================================
# §5.1 — BP-trajectory LMM
# =============================================================================
cat("\n===== §5.1 BP-trajectory LMM =====\n")

s51_base <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female"))

s51_long <- function(df) {
  df %>%
    select(STUDY_ID, BASELINE_AHA, AGE_AT_START, SEX, BASELINE_BMI, HTN_MED_BASELINE,
           BASELINE_SBP, MO3_SBP, MO6_SBP,
           BASELINE_DBP, MO3_DBP, MO6_DBP,
           BASELINE_TIME_OF_DAY, MO3_TIME_OF_DAY, MO6_TIME_OF_DAY,
           BASELINE_SEASON, MO3_SEASON, MO6_SEASON) %>%
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
      SEASON       = factor(SEASON, levels = c("winter", "spring", "summer", "fall")),
      HTN_BASELINE_ON = factor(
        if_else(HTN_MED_BASELINE %in% c("stable_on_meds", "recent_initiator"),
                "on_at_baseline", "not_on_at_baseline"),
        levels = c("not_on_at_baseline", "on_at_baseline"))
    )
}

s51_fit_lmm <- function(outcome_col, adjusted, data) {
  if (adjusted) {
    f <- as.formula(paste0(
      outcome_col, " ~ Time * BASELINE_AHA",
      " + AGE_AT_START + SEX + BASELINE_BMI + HTN_BASELINE_ON",
      " + sin(2 * pi * TIME_OF_DAY / 24) + cos(2 * pi * TIME_OF_DAY / 24)",
      " + SEASON",
      " + Time:AGE_AT_START + Time:SEX + Time:BASELINE_BMI + Time:HTN_BASELINE_ON",
      " + (1 | STUDY_ID)"))
    d <- data %>% filter(!is.na(.data[[outcome_col]]), !is.na(BASELINE_AHA),
                         !is.na(AGE_AT_START), !is.na(SEX), !is.na(BASELINE_BMI),
                         !is.na(HTN_BASELINE_ON), !is.na(TIME_OF_DAY), !is.na(SEASON))
  } else {
    f <- as.formula(paste0(outcome_col, " ~ Time * BASELINE_AHA + (1 | STUDY_ID)"))
    d <- data %>% filter(!is.na(.data[[outcome_col]]), !is.na(BASELINE_AHA))
  }
  lmer(f, data = d, REML = TRUE)
}

s51_extract <- function(fit) {
  emm <- suppressMessages(emmeans(fit, ~ Time | BASELINE_AHA))
  delta <- suppressMessages(contrast(emm, method = "trt.vs.ctrl", adjust = "none"))
  cross <- suppressMessages(contrast(delta, "pairwise", by = "contrast", adjust = "none"))
  list(within = tidy_contrast(delta), cross = tidy_contrast(cross))
}

s51_run <- function(df) {
  long <- s51_long(df)
  models <- list(sbp_unadj = s51_fit_lmm("SBP", FALSE, long),
                 sbp_adj   = s51_fit_lmm("SBP", TRUE,  long),
                 dbp_unadj = s51_fit_lmm("DBP", FALSE, long),
                 dbp_adj   = s51_fit_lmm("DBP", TRUE,  long))
  within <- imap_dfr(models, ~ s51_extract(.x)$within %>% mutate(model = .y))
  cross  <- imap_dfr(models, ~ s51_extract(.x)$cross  %>% mutate(model = .y))
  n_by_stratum <- imap_dfr(models, ~ model.frame(.x) %>% distinct(STUDY_ID, BASELINE_AHA) %>%
                             count(BASELINE_AHA, name = "n_users") %>% mutate(model = .y))
  list(within = within, cross = cross, n_by_stratum = n_by_stratum,
       n_users = map_int(models, ~ as.integer(lme4::ngrps(.x))))
}

s51 <- map(cohort_variants(s51_base), s51_run)
cat("§5.1 n_users (sbp_adj) — primary:", s51$primary$n_users[["sbp_adj"]],
    " mid_excluded:", s51$mid_excluded$n_users[["sbp_adj"]], "\n")

# Fidelity gate: our primary within-stratum Δ vs frozen 5_1_delta_contrasts.csv
frozen_51 <- read_csv(file.path(p51, "5_1_delta_contrasts.csv"), show_col_types = FALSE)
fid_51 <- s51$primary$within %>%
  select(model, contrast, BASELINE_AHA, est_ours = estimate) %>%
  inner_join(frozen_51 %>% select(model, contrast, BASELINE_AHA, est_frozen = estimate),
             by = c("model", "contrast", "BASELINE_AHA")) %>%
  mutate(absdiff = abs(est_ours - est_frozen))
cat(sprintf("§5.1 FIDELITY: matched %d/%d within-stratum contrasts; max |Δestimate| = %.2e\n",
            nrow(fid_51), nrow(frozen_51), max(fid_51$absdiff)))
stopifnot(nrow(fid_51) == nrow(frozen_51), max(fid_51$absdiff) < 1e-6)

# Write mid-excluded contrast tables
write_csv(s51$mid_excluded$within, file.path(results_dir, "5_1_mid_excluded_delta_contrasts.csv"))
write_csv(s51$mid_excluded$cross,  file.path(results_dir, "5_1_mid_excluded_cross_stratum_diffs.csv"))
write_csv(bind_rows(s51$primary$n_by_stratum %>% mutate(cohort = "primary"),
                    s51$mid_excluded$n_by_stratum %>% mutate(cohort = "mid_excluded")),
          file.path(results_dir, "5_1_mid_excluded_n_by_stratum.csv"))

# Discordance (within-stratum Δ family; keyed on model|contrast|stratum; Holm p)
disc_all[["5_1_within"]] <- discordance(
  standardize(s51$primary$within,    c("model","contrast","BASELINE_AHA"), "estimate", "p.holm", "5.1", "within_stratum_delta"),
  standardize(s51$mid_excluded$within, c("model","contrast","BASELINE_AHA"), "estimate", "p.holm", "5.1", "within_stratum_delta"))
disc_all[["5_1_cross"]] <- discordance(
  standardize(s51$primary$cross,    c("model","contrast1","contrast"), "estimate", "p.holm", "5.1", "cross_stratum_diff"),
  standardize(s51$mid_excluded$cross, c("model","contrast1","contrast"), "estimate", "p.holm", "5.1", "cross_stratum_diff"))

cat("§5.1 discordant (within):", sum(disc_all[["5_1_within"]]$discordant, na.rm = TRUE),
    " (cross):", sum(disc_all[["5_1_cross"]]$discordant, na.rm = TRUE), "\n")

AHA3 <- c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")

# Standardize a coefficient table (focal terms only) for discordance.
std_coef <- function(df, key_cols, section, family, focal, pcol = "p.value") {
  df %>% filter(term %in% focal) %>%
    standardize(key_cols, "estimate", pcol, section, family)
}
# Fidelity helper: compare our-primary vs frozen on a numeric column.
fidelity <- function(ours, frozen, keys, col_ours, col_frozen, label, tol = 1e-5) {
  m <- inner_join(ours   %>% select(all_of(keys), .o = all_of(col_ours)),
                  frozen %>% select(all_of(keys), .f = all_of(col_frozen)),
                  by = keys)
  md <- max(abs(m$.o - m$.f), na.rm = TRUE)
  cat(sprintf("%s FIDELITY: matched %d rows; max |Δ| = %.2e\n", label, nrow(m), md))
  stopifnot(nrow(m) > 0, md < tol)
}

# =============================================================================
# §6.2 — WHOOP-signal responder logistic (10 predictors, dual outcome)
# =============================================================================
cat("\n===== §6.2 WHOOP responder logistic =====\n")
p62 <- file.path(project_dir, "Analyses", "6_2", "derived", "results")
preds62 <- c("SLEEP_DURATION_BASELINE_MEAN", "SLEEP_CONSISTENCY_BASELINE_MEAN",
             "STRAIN_BASELINE_MEAN", "HRV_BASELINE_MEAN", "RHR_BASELINE_MEAN",
             "SLEEP_DURATION_STUDY_MEAN", "SLEEP_CONSISTENCY_STUDY_MEAN",
             "STRAIN_STUDY_MEAN", "HRV_STUDY_MEAN", "RHR_STUDY_MEAN")
cov62 <- c("AGE_AT_START", "SEX", "BASELINE_BMI", "HTN_BASELINE_ON")

s62_base <- user_frame %>%
  mutate(DBP_RESPONDER_6MO = as.integer((MO6_DBP - BASELINE_DBP) <= -5),
         HTN_BASELINE_ON   = htn_baseline_on(HTN_MED_BASELINE)) %>%
  filter(BASELINE_AHA %in% AHA3, SEX %in% c("male", "female"), BASELINE_N_DAYS >= 8)

s62_run <- function(df) {
  out <- list()
  for (otc in c("SBP", "DBP")) {
    base_col <- paste0("BASELINE_", otc); resp <- paste0(otc, "_RESPONDER_6MO")
    rhs <- paste(c(base_col, preds62, cov62), collapse = " + ")
    d <- df %>% filter(if_all(all_of(c(preds62, cov62, base_col, resp)), ~ !is.na(.)))
    fit <- glm(as.formula(paste0(resp, " ~ ", rhs)), family = binomial(), data = d)
    out[[otc]] <- broom::tidy(fit) %>%
      mutate(OR = exp(estimate), outcome = paste0(otc, "_responder"), n = nrow(d))
  }
  bind_rows(out)
}
s62 <- map(cohort_variants(s62_base), s62_run)
fidelity(s62$primary,
         read_csv(file.path(p62, "6_2_coefficients_primary.csv"), show_col_types = FALSE) %>% filter(cohort == "Primary"),
         c("outcome", "term"), "OR", "OR", "§6.2 logistic")
write_csv(s62$mid_excluded, file.path(results_dir, "6_2_mid_excluded_coefficients.csv"))
disc_all[["6_2"]] <- discordance(
  std_coef(s62$primary,    c("outcome", "term"), "6.2", "whoop_responder", preds62),
  std_coef(s62$mid_excluded, c("outcome", "term"), "6.2", "whoop_responder", preds62))

# =============================================================================
# §7.2 — Within/between daily LMM (cluster-robust CR1)
# =============================================================================
cat("\n===== §7.2 within/between daily LMM (this is the slow one) =====\n")
p72 <- file.path(project_dir, "Analyses", "7_2", "derived", "results")
daily72 <- c(RHR = "RHR_LOOKBACK_MEAN", HRV = "HRV_LOOKBACK_MEAN",
             Strain = "STRAIN_LOOKBACK_MEAN", SleepDuration = "SLEEP_DURATION_LOOKBACK_MEAN",
             SleepConsistency = "SLEEP_CONSISTENCY_LOOKBACK_MEAN")
wb_terms72 <- as.vector(rbind(paste0(names(daily72), "_b"), paste0(names(daily72), "_w")))
bp_frame_j <- read_rds(reading_path)
stopifnot(!any(is.na(bp_frame_j$HTN_MED_BASELINE)))

s72_run <- function(df) {
  cohort <- df %>%
    filter(BASELINE_AHA %in% AHA3, SEX %in% c("male", "female"),
           !IS_BASELINE_READING, PERIOD == "DURING") %>%
    group_by(STUDY_ID) %>% filter(n() >= 2) %>% ungroup()
  for (m in names(daily72)) {
    cohort <- cohort %>%
      group_by(STUDY_ID) %>%
      mutate(!!paste0(m, "_b") := mean(.data[[daily72[[m]]]], na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(!!paste0(m, "_w") := .data[[daily72[[m]]]] - .data[[paste0(m, "_b")]])
  }
  cohort <- cohort %>%
    mutate(HTN_BASELINE_ON = htn_baseline_on(HTN_MED_BASELINE)) %>%
    filter(if_all(all_of(c(wb_terms72, "AGE_AT_START", "SEX", "BASELINE_BMI",
                           "HTN_BASELINE_ON", "SBP_MINUS_BASELINE", "DBP_MINUS_BASELINE")), ~ !is.na(.))) %>%
    mutate(SEX = factor(SEX))
  for (m in names(daily72)) {
    b_col <- paste0(m, "_b"); w_col <- paste0(m, "_w")
    ub <- distinct(cohort[, c("STUDY_ID", b_col)])
    cohort[[b_col]] <- (cohort[[b_col]] - mean(ub[[b_col]], na.rm = TRUE)) / sd(ub[[b_col]], na.rm = TRUE)
    cohort[[w_col]] <- as.numeric(scale(cohort[[w_col]]))
  }
  rhs <- paste(c("AGE_AT_START", "SEX", "BASELINE_BMI", "HTN_BASELINE_ON", wb_terms72), collapse = " + ")
  out <- list()
  for (otc in c("SBP", "DBP")) {
    oc <- paste0(otc, "_MINUS_BASELINE")
    fit <- lmer(as.formula(paste0(oc, " ~ ", rhs, " + (1 | STUDY_ID)")), data = cohort, REML = TRUE)
    vcr <- clubSandwich::vcovCR(fit, cluster = cohort$STUDY_ID, type = "CR1")
    ct <- clubSandwich::coef_test(fit, vcov = vcr, test = "z", coefs = "All")
    out[[otc]] <- tibble(outcome = paste0(otc, "-baseline"), term = rownames(ct),
                         estimate = ct$beta, z_stat = ct$tstat, p_value = ct$p_z)
  }
  bind_rows(out)
}
s72 <- map(cohort_variants(bp_frame_j), s72_run)
fidelity(s72$primary, read_csv(file.path(p72, "7_2_coefficients.csv"), show_col_types = FALSE),
         c("outcome", "term"), "estimate", "estimate", "§7.2 coefficients")
write_csv(s72$mid_excluded, file.path(results_dir, "7_2_mid_excluded_coefficients.csv"))
disc_all[["7_2"]] <- discordance(
  std_coef(s72$primary,    c("outcome", "term"), "7.2", "within_between", wb_terms72, pcol = "p_value"),
  std_coef(s72$mid_excluded, c("outcome", "term"), "7.2", "within_between", wb_terms72, pcol = "p_value"))

# =============================================================================
# §8 — Digital engagement responder logistic (separate + joint)
# =============================================================================
cat("\n===== §8 engagement logistic =====\n")
p8 <- file.path(project_dir, "Analyses", "8_engagement", "derived", "results")
eng_preds <- c("PCT_DAYS_JOURNALED_per10", "LOG_AVG_DAILY_TIME_IN_APP_MIN",
               "PCT_DAYS_ENGAGED_NOTIFICATIONS_per10", "PCT_DAYS_USED_COACH_per10",
               "AVG_SLEEP_NEED_PCT_per10", "PCT_DAYS_MET_OPTIMAL_STRAIN_per10",
               "LOG_N_BP_READINGS_STUDY")
cov8 <- c("AGE_AT_START", "SEX", "BASELINE_BMI", "HTN_BASELINE_ON")
s8_base <- user_frame %>%
  filter(BASELINE_AHA %in% AHA3, SEX %in% c("male", "female")) %>%
  mutate(DBP_RESPONDER_6MO = as.integer((MO6_DBP - BASELINE_DBP) <= -5),
         SEX = factor(SEX),
         PCT_DAYS_JOURNALED_per10             = PCT_DAYS_JOURNALED / 10,
         PCT_DAYS_ENGAGED_NOTIFICATIONS_per10 = PCT_DAYS_ENGAGED_NOTIFICATIONS / 10,
         PCT_DAYS_USED_COACH_per10            = PCT_DAYS_USED_COACH / 10,
         AVG_SLEEP_NEED_PCT_per10             = AVG_SLEEP_NEED_PCT / 10,
         PCT_DAYS_MET_OPTIMAL_STRAIN_per10    = PCT_DAYS_MET_OPTIMAL_STRAIN / 10,
         HTN_BASELINE_ON = htn_baseline_on(HTN_MED_BASELINE))

s8_run <- function(df) {
  sep <- list(); joint <- list()
  for (otc in c("SBP", "DBP")) {
    base_col <- paste0("BASELINE_", otc); resp <- paste0(otc, "_RESPONDER_6MO")
    for (pred in eng_preds) {
      d <- df %>% filter(if_all(all_of(c(pred, cov8, base_col, resp)), ~ !is.na(.)))
      fit <- glm(as.formula(paste0(resp, " ~ ", paste(c(base_col, pred, cov8), collapse = " + "))),
                 family = binomial(), data = d)
      sep[[paste(otc, pred)]] <- broom::tidy(fit) %>% filter(term == pred) %>%
        mutate(OR = exp(estimate), frame = "sep", predictor = pred, outcome = otc, n = nrow(d))
    }
    dj <- df %>% filter(if_all(all_of(c(eng_preds, cov8, base_col, resp)), ~ !is.na(.)))
    fitj <- glm(as.formula(paste0(resp, " ~ ", paste(c(base_col, eng_preds, cov8), collapse = " + "))),
                family = binomial(), data = dj)
    joint[[otc]] <- broom::tidy(fitj) %>% filter(term %in% eng_preds) %>%
      mutate(OR = exp(estimate), frame = "joint", predictor = "ALL", outcome = otc, n = nrow(dj))
  }
  list(sep = bind_rows(sep), joint = bind_rows(joint))
}
s8 <- map(cohort_variants(s8_base), s8_run)
fidelity(s8$primary$sep, read_csv(file.path(p8, "8_coefficients_separate_primary.csv"), show_col_types = FALSE),
         c("predictor", "outcome", "term"), "OR", "OR", "§8 separate")
fidelity(s8$primary$joint, read_csv(file.path(p8, "8_coefficients_joint_primary.csv"), show_col_types = FALSE),
         c("outcome", "term"), "OR", "OR", "§8 joint")
write_csv(s8$mid_excluded$sep,   file.path(results_dir, "8_mid_excluded_coefficients_separate.csv"))
write_csv(s8$mid_excluded$joint, file.path(results_dir, "8_mid_excluded_coefficients_joint.csv"))
disc_all[["8_sep"]] <- discordance(
  std_coef(s8$primary$sep,    c("outcome", "term"), "8", "engagement_separate", eng_preds),
  std_coef(s8$mid_excluded$sep, c("outcome", "term"), "8", "engagement_separate", eng_preds))
disc_all[["8_joint"]] <- discordance(
  std_coef(s8$primary$joint,    c("outcome", "term"), "8", "engagement_joint", eng_preds),
  std_coef(s8$mid_excluded$joint, c("outcome", "term"), "8", "engagement_joint", eng_preds))

# =============================================================================
# Consolidated discordance report
# =============================================================================
cat("\n===== Building consolidated discordance report =====\n")
disc_report <- bind_rows(disc_all) %>%
  mutate(
    both_nonsig = !(prim_sig %in% TRUE) & !(naive_sig %in% TRUE),
    # substantive = a real significance change; sign flips on terms that are
    # non-significant in BOTH cohorts are near-zero estimates wobbling sign (noise).
    flag_type = case_when(
      !discordant                              ~ "concordant",
      sig_flip &  (naive_sig %in% TRUE)        ~ "gained significance",
      sig_flip & !(naive_sig %in% TRUE)        ~ "lost significance",
      sign_flip & both_nonsig                  ~ "sign flip (both non-sig; near-zero)",
      sign_flip                                ~ "sign flip",
      TRUE                                     ~ "other"),
    substantive = discordant & flag_type %in% c("gained significance", "lost significance", "sign flip")
  ) %>%
  select(section, family, key, prim_stat, prim_p, prim_sig,
         naive_stat, naive_p, naive_sig, sig_flip, sign_flip, discordant, flag_type, substantive)
write_csv(disc_report, file.path(results_dir, "mid_excluded_discordance_report.csv"))

n_terms <- nrow(disc_report)
n_disc  <- sum(disc_report$discordant,  na.rm = TRUE)
n_subst <- sum(disc_report$substantive, na.rm = TRUE)
cat(sprintf("\n=== Mid-study-initiator-excluded sensitivity complete ===\n"))
cat(sprintf("Compared %d focal terms/contrasts across §5.1, §6.2, §7.2, §8.\n", n_terms))
cat(sprintf("DISCORDANT by rule (sig-flip at a=%.3f OR sign reversal): %d\n", ALPHA, n_disc))
cat(sprintf("  of which SUBSTANTIVE (>=1 cohort significant): %d\n", n_subst))
cat(sprintf("  sign-flips on null terms (both non-sig; noise):  %d\n", n_disc - n_subst))
cat("\nDiscordance count by section:\n")
print(disc_report %>% group_by(section) %>%
        summarise(n = n(), discordant = sum(discordant, na.rm = TRUE),
                  substantive = sum(substantive, na.rm = TRUE), .groups = "drop"))
if (n_disc > 0) {
  cat("\n--- All discordant findings (most meaningful first) ---\n")
  print(disc_report %>% filter(discordant) %>%
          arrange(desc(substantive), section) %>%
          mutate(prim_p = signif(prim_p, 3), naive_p = signif(naive_p, 3)) %>%
          select(section, family, key, prim_p, prim_sig, naive_p, naive_sig, flag_type),
        n = Inf, width = Inf)
}
cat("\nWrote:", file.path(results_dir, "mid_excluded_discordance_report.csv"), "\n")
