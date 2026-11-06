# glp1_naive_sensitivity.R
# -----------------------------------------------------------------------------
# GLP-1-naive sensitivity analysis (single-script mirror of the antihypertensive
# HTN-naive sensitivity, but EXCLUDING GLP-1 users).
#
# Decision record (2026-07-20):
#   - Cohort DIRECTION: EXCLUDE GLP-1 users. The GLP-1-naive cohort = the §2.2
#     primary cohort with GLP1_EVER == TRUE users removed. This mirrors the
#     htn_naive design (which removes med users), so the model specs are the
#     PRIMARY specs unchanged (HTN_BASELINE_ON still varies and is KEPT; GLP1 is
#     constant/FALSE and is not a covariate). The only change vs primary is the
#     added filter(!GLP1_EVER).
#   - SCOPE (per written Methods umbrella + Results): §5.1, §6.1, §6.2, §7.1,
#     §7.2, §8, plus a descriptive Table 1.
#   - DISCORDANCE rule: a finding is discordant if a term/contrast flips
#     significance at alpha = 0.005 OR reverses sign vs the primary cohort.
#
# Design: each chapter defines ONE fit routine, run on (a) the primary §2.2
# cohort and (b) the GLP-1-naive cohort. A fidelity gate checks that our primary
# re-run reproduces the frozen primary CSVs (max abs diff on key estimates); if
# it matches, the GLP-1-naive run using the identical code path is trustworthy.
# Discordance is then our-naive vs our-primary (identical to naive-vs-frozen).
#
# Run from repo root:  Rscript PRIMARY/Analyses/glp1_naive_sensitivity/glp1_naive_sensitivity.R
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
glp1_path    <- file.path(project_dir, "Build", "derived", "09_glp1_flag.rds")

out_dir     <- file.path(project_dir, "Analyses", "glp1_naive_sensitivity", "derived")
results_dir <- file.path(out_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Frozen primary result roots (for fidelity gate)
p51 <- file.path(project_dir, "Analyses", "5_1", "derived", "results")

# --- Load frames + attach GLP1_EVER -----------------------------------------
glp1 <- read_rds(glp1_path) %>% select(STUDY_ID, GLP1_EVER)

user_frame <- read_rds(frame_path) %>%
  left_join(glp1, by = "STUDY_ID")
stopifnot(!any(is.na(user_frame$GLP1_EVER)))   # roster universe fully covered

cat("User frame:", nrow(user_frame), "rows;  GLP1_EVER TRUE:",
    sum(user_frame$GLP1_EVER), " FALSE:", sum(!user_frame$GLP1_EVER), "\n")

# Cohort builders. `base` applies the §2.2 + sex filter; the two cohorts differ
# ONLY by the GLP-1 exclusion.
cohort_variants <- function(base_df) {
  list(primary    = base_df,
       glp1_naive = base_df %>% filter(!GLP1_EVER))
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
  list(within = within, cross = cross,
       n_users = map_int(models, ~ as.integer(lme4::ngrps(.x))))
}

s51 <- map(cohort_variants(s51_base), s51_run)
cat("§5.1 n_users (sbp_adj) — primary:", s51$primary$n_users[["sbp_adj"]],
    " glp1_naive:", s51$glp1_naive$n_users[["sbp_adj"]], "\n")

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

# Write GLP-1-naive contrast tables
write_csv(s51$glp1_naive$within, file.path(results_dir, "5_1_glp1_naive_delta_contrasts.csv"))
write_csv(s51$glp1_naive$cross,  file.path(results_dir, "5_1_glp1_naive_cross_stratum_diffs.csv"))

# Discordance (within-stratum Δ family; keyed on model|contrast|stratum; Holm p)
disc_all[["5_1_within"]] <- discordance(
  standardize(s51$primary$within,    c("model","contrast","BASELINE_AHA"), "estimate", "p.holm", "5.1", "within_stratum_delta"),
  standardize(s51$glp1_naive$within, c("model","contrast","BASELINE_AHA"), "estimate", "p.holm", "5.1", "within_stratum_delta"))
disc_all[["5_1_cross"]] <- discordance(
  standardize(s51$primary$cross,    c("model","contrast1","contrast"), "estimate", "p.holm", "5.1", "cross_stratum_diff"),
  standardize(s51$glp1_naive$cross, c("model","contrast1","contrast"), "estimate", "p.holm", "5.1", "cross_stratum_diff"))

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
# §6.1 — LE8 lifestyle regressions (continuous ANCOVA + responder logistic)
# =============================================================================
cat("\n===== §6.1 LE8 regressions =====\n")
p61 <- file.path(project_dir, "Analyses", "6_1", "derived", "results")
le8_primary <- c("SLEEP_DURATION_STUDY_MEAN", "STRAIN_STUDY_MEAN", "LOG_TOBACCO_NICOTINE_USE_DAYS")
cov61 <- c("AGE_AT_START", "SEX", "BASELINE_BMI", "HTN_BASELINE_ON",
           "SLEEP_DURATION_BASELINE_MEAN", "STRAIN_BASELINE_MEAN")
variants61 <- list(Primary     = le8_primary,
                   PlusDweight = c(le8_primary, "DELTA_WEIGHT_KG"),
                   PlusDHDL    = c(le8_primary, "DELTA_HDL_CHOLESTEROL"),
                   PlusDNonHDL = c(le8_primary, "DELTA_NON_HDL_CHOLESTEROL"),
                   PlusDTG     = c(le8_primary, "DELTA_TRIGLYCERIDES"))
focal61 <- c(le8_primary, "DELTA_WEIGHT_KG", "DELTA_HDL_CHOLESTEROL",
             "DELTA_NON_HDL_CHOLESTEROL", "DELTA_TRIGLYCERIDES")

s61_base <- user_frame %>%
  filter(BASELINE_AHA %in% AHA3, SEX %in% c("male", "female")) %>%
  mutate(DBP_RESPONDER_6MO = as.integer((MO6_DBP - BASELINE_DBP) <= -5),
         HTN_BASELINE_ON   = htn_baseline_on(HTN_MED_BASELINE))

s61_run <- function(df) {
  cont <- list(); logi <- list()
  for (vn in names(variants61)) {
    preds <- variants61[[vn]]
    for (otc in c("SBP", "DBP")) {
      base_col <- paste0("BASELINE_", otc); out_col <- paste0("MO6_", otc)
      rhs <- paste(c(base_col, preds, cov61), collapse = " + ")
      dC <- df %>% filter(if_all(all_of(c(preds, cov61, base_col, out_col)), ~ !is.na(.)))
      fitC <- lm(as.formula(paste0(out_col, " ~ ", rhs)), data = dC)
      cont[[paste(otc, vn, sep = "__")]] <- broom::tidy(fitC) %>%
        mutate(outcome = otc, variant = vn, n = nrow(dC))
      resp <- paste0(otc, "_RESPONDER_6MO")
      dL <- df %>% filter(if_all(all_of(c(preds, cov61, base_col, resp)), ~ !is.na(.)))
      fitL <- glm(as.formula(paste0(resp, " ~ ", rhs)), family = binomial(), data = dL)
      logi[[paste(otc, vn, sep = "__")]] <- broom::tidy(fitL) %>%
        mutate(OR = exp(estimate), responder_outcome = paste0(otc, "_responder"), variant = vn, n = nrow(dL))
    }
  }
  list(continuous = bind_rows(cont), logistic = bind_rows(logi))
}
s61 <- map(cohort_variants(s61_base), s61_run)
fidelity(s61$primary$continuous, read_csv(file.path(p61, "6_1_coefficients_continuous.csv"), show_col_types = FALSE),
         c("outcome", "variant", "term"), "estimate", "estimate", "§6.1 continuous")
fidelity(s61$primary$logistic, read_csv(file.path(p61, "6_1_coefficients_logistic.csv"), show_col_types = FALSE),
         c("responder_outcome", "variant", "term"), "OR", "OR", "§6.1 logistic")
write_csv(s61$glp1_naive$continuous, file.path(results_dir, "6_1_glp1_naive_coefficients_continuous.csv"))
write_csv(s61$glp1_naive$logistic,   file.path(results_dir, "6_1_glp1_naive_coefficients_logistic.csv"))
disc_all[["6_1_cont"]] <- discordance(
  std_coef(s61$primary$continuous,    c("outcome", "variant", "term"), "6.1", "le8_continuous", focal61),
  std_coef(s61$glp1_naive$continuous, c("outcome", "variant", "term"), "6.1", "le8_continuous", focal61))
disc_all[["6_1_logi"]] <- discordance(
  std_coef(s61$primary$logistic,    c("responder_outcome", "variant", "term"), "6.1", "le8_logistic", focal61),
  std_coef(s61$glp1_naive$logistic, c("responder_outcome", "variant", "term"), "6.1", "le8_logistic", focal61))

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
write_csv(s62$glp1_naive, file.path(results_dir, "6_2_glp1_naive_coefficients.csv"))
disc_all[["6_2"]] <- discordance(
  std_coef(s62$primary,    c("outcome", "term"), "6.2", "whoop_responder", preds62),
  std_coef(s62$glp1_naive, c("outcome", "term"), "6.2", "whoop_responder", preds62))

# =============================================================================
# §7.1 — WHOOP-metric trajectory LMM by responder status (dual-primary)
# =============================================================================
cat("\n===== §7.1 WHOOP trajectory LMM =====\n")
p71 <- file.path(project_dir, "Analyses", "7_1", "derived", "results")
metrics71 <- c(RHR = "RHR", SleepDuration = "SLEEP_DURATION",
               SleepConsistency = "SLEEP_CONSISTENCY", Strain = "STRAIN", HRV = "HRV")
uf71 <- user_frame %>%
  mutate(DBP_RESPONDER_6MO = as.integer((MO6_DBP - BASELINE_DBP) <= -5),
         HTN_BASELINE_ON   = htn_baseline_on(HTN_MED_BASELINE))

s71_long <- function(prefix, d) {
  cols <- paste0(prefix, c("_BASELINE_MEAN", "_MO3_MEAN", "_MO6_MEAN"))
  d %>%
    select(STUDY_ID, responder_group, AGE_AT_START, SEX, BASELINE_BMI, HTN_BASELINE_ON, all_of(cols)) %>%
    rename(baseline = !!sym(cols[1]), `3mo` = !!sym(cols[2]), `6mo` = !!sym(cols[3])) %>%
    pivot_longer(c(baseline, `3mo`, `6mo`), names_to = "Time", values_to = "value") %>%
    mutate(Time = factor(Time, levels = c("baseline", "3mo", "6mo")), SEX = factor(SEX)) %>%
    filter(!is.na(value), !is.na(AGE_AT_START), !is.na(BASELINE_BMI))
}
s71_run <- function(df) {
  delta_rows <- list(); anova_rows <- list()
  for (rd in c("SBP", "DBP")) {
    resp_col <- paste0(rd, "_RESPONDER_6MO")
    cohort <- df %>%
      filter(BASELINE_AHA %in% AHA3, SEX %in% c("male", "female"), !is.na(.data[[resp_col]])) %>%
      mutate(responder_group = factor(ifelse(.data[[resp_col]] == 1L, "responder", "non_responder"),
                                      levels = c("non_responder", "responder")))
    for (m in names(metrics71)) {
      long_d <- s71_long(metrics71[[m]], cohort)
      fit <- lmer(value ~ Time * responder_group + AGE_AT_START + SEX + BASELINE_BMI +
                    HTN_BASELINE_ON + (1 | STUDY_ID), data = long_d, REML = TRUE)
      emm <- suppressMessages(emmeans(fit, ~ Time | responder_group))
      delta_rows[[paste(rd, m)]] <- tidy_contrast(suppressMessages(
        contrast(emm, method = "trt.vs.ctrl", adjust = "none"))) %>%
        mutate(responder_def = rd, metric = m)
      a <- suppressMessages(anova(fit, type = 3))
      anova_rows[[paste(rd, m)]] <- as_tibble(a, rownames = "term") %>%
        filter(term == "Time:responder_group") %>%
        mutate(responder_def = rd, metric = m)
    }
  }
  list(delta = bind_rows(delta_rows), anova = bind_rows(anova_rows))
}
s71 <- map(cohort_variants(uf71 %>% mutate(.dummy = TRUE) %>% select(-.dummy)), s71_run)
fidelity(s71$primary$delta, read_csv(file.path(p71, "7_1_delta_contrasts.csv"), show_col_types = FALSE),
         c("responder_def", "metric", "contrast", "responder_group"), "estimate", "estimate", "§7.1 delta")
write_csv(s71$glp1_naive$delta, file.path(results_dir, "7_1_glp1_naive_delta_contrasts.csv"))
write_csv(s71$glp1_naive$anova, file.path(results_dir, "7_1_glp1_naive_anova.csv"))
disc_all[["7_1_delta"]] <- discordance(
  standardize(s71$primary$delta,    c("responder_def", "metric", "contrast", "responder_group"), "estimate", "p.holm", "7.1", "trajectory_delta"),
  standardize(s71$glp1_naive$delta, c("responder_def", "metric", "contrast", "responder_group"), "estimate", "p.holm", "7.1", "trajectory_delta"))
disc_all[["7_1_interaction"]] <- discordance(
  standardize(s71$primary$anova,    c("responder_def", "metric"), "F value", "Pr(>F)", "7.1", "time_x_responder"),
  standardize(s71$glp1_naive$anova, c("responder_def", "metric"), "F value", "Pr(>F)", "7.1", "time_x_responder"))

# =============================================================================
# §7.2 — Within/between daily LMM (cluster-robust CR1)
# =============================================================================
cat("\n===== §7.2 within/between daily LMM (this is the slow one) =====\n")
p72 <- file.path(project_dir, "Analyses", "7_2", "derived", "results")
daily72 <- c(RHR = "RHR_LOOKBACK_MEAN", HRV = "HRV_LOOKBACK_MEAN",
             Strain = "STRAIN_LOOKBACK_MEAN", SleepDuration = "SLEEP_DURATION_LOOKBACK_MEAN",
             SleepConsistency = "SLEEP_CONSISTENCY_LOOKBACK_MEAN")
wb_terms72 <- as.vector(rbind(paste0(names(daily72), "_b"), paste0(names(daily72), "_w")))
bp_frame_j <- read_rds(reading_path) %>% left_join(glp1, by = "STUDY_ID")
stopifnot(!any(is.na(bp_frame_j$GLP1_EVER)))

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
write_csv(s72$glp1_naive, file.path(results_dir, "7_2_glp1_naive_coefficients.csv"))
disc_all[["7_2"]] <- discordance(
  std_coef(s72$primary,    c("outcome", "term"), "7.2", "within_between", wb_terms72, pcol = "p_value"),
  std_coef(s72$glp1_naive, c("outcome", "term"), "7.2", "within_between", wb_terms72, pcol = "p_value"))

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
write_csv(s8$glp1_naive$sep,   file.path(results_dir, "8_glp1_naive_coefficients_separate.csv"))
write_csv(s8$glp1_naive$joint, file.path(results_dir, "8_glp1_naive_coefficients_joint.csv"))
disc_all[["8_sep"]] <- discordance(
  std_coef(s8$primary$sep,    c("outcome", "term"), "8", "engagement_separate", eng_preds),
  std_coef(s8$glp1_naive$sep, c("outcome", "term"), "8", "engagement_separate", eng_preds))
disc_all[["8_joint"]] <- discordance(
  std_coef(s8$primary$joint,    c("outcome", "term"), "8", "engagement_joint", eng_preds),
  std_coef(s8$glp1_naive$joint, c("outcome", "term"), "8", "engagement_joint", eng_preds))

# =============================================================================
# Table 1 — descriptive between-stratum tests (no sex restriction)
# =============================================================================
cat("\n===== Table 1 descriptives =====\n")
t1_prim  <- user_frame %>% filter(BASELINE_AHA %in% AHA3)
t1_naive <- t1_prim %>% filter(!GLP1_EVER)
write_csv(make_char_table(t1_naive, AHA3), file.path(results_dir, "table_1_glp1_naive.csv"))
pv_prim  <- stratum_pvals(t1_prim)
pv_naive <- stratum_pvals(t1_naive)
disc_all[["table1"]] <- tibble(
  section = "Table1", family = "between_stratum_p", key = names(pv_prim),
  prim_stat = 1, prim_p = as.numeric(pv_prim), naive_stat = 1, naive_p = as.numeric(pv_naive)) %>%
  mutate(prim_sig = prim_p < ALPHA, naive_sig = naive_p < ALPHA,
         sig_flip = (prim_sig != naive_sig) %in% TRUE, sign_flip = FALSE,
         discordant = sig_flip)

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
write_csv(disc_report, file.path(results_dir, "glp1_naive_discordance_report.csv"))

n_terms <- nrow(disc_report)
n_disc  <- sum(disc_report$discordant,  na.rm = TRUE)
n_subst <- sum(disc_report$substantive, na.rm = TRUE)
cat(sprintf("\n=== GLP-1-naive sensitivity complete ===\n"))
cat(sprintf("Compared %d focal terms/contrasts across §5.1, §6.1, §6.2, §7.1, §7.2, §8, Table 1.\n", n_terms))
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
cat("\nWrote:", file.path(results_dir, "glp1_naive_discordance_report.csv"), "\n")
