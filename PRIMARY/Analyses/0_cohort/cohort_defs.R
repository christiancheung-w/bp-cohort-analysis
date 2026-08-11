# cohort_defs.R
# -----------------------------------------------------------------------------
# Central cohort definitions + Table-1 characteristic builder + between-stratum
# tests for the BP dID cohort chapter.
#
# This is a SOURCED LIBRARY — it defines objects/functions and writes nothing.
# Run the 0_cohort table scripts (which source this) to produce outputs.
#
# Rationale: the AHA-stratum cohort predicate was previously copy-pasted across
# ~15 analysis scripts. This module is the single canonical definition. NOTE
# (per review decision 2026-06): we are NOT retrofitting the existing frozen
# chapters to consume this — it backs the new 0_cohort reporting only.
#
# Exposes:
#   AHA_HYPERTENSIVE / AHA_FULL  ordered stratum vocab
#   COHORTS                      named registry: label, filter, strata
#   apply_cohort(d, name)        subset the user frame by a named cohort
#   msd(), n_pct()               cell formatters
#   p_continuous(), p_categorical(), fmt_p()   between-stratum tests + formatter
#   compute_one_stratum(d)       per-group characteristic row
#   make_char_table(d, strata)   Characteristic x {Overall, strata, P-value}
# -----------------------------------------------------------------------------

library(tidyverse)

# --- AHA stratum vocabulary (display order) ----------------------------------
AHA_HYPERTENSIVE <- c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")
AHA_FULL         <- c("Normal", AHA_HYPERTENSIVE)

# --- Cohort registry ---------------------------------------------------------
# Each entry: human label, a filter predicate over the user frame, and the
# ordered stratum levels shown as columns. Table 1 (paper) = `hypertensive`;
# `full` (incl. Normal) and `htn_naive` are the sensitivity cohorts.
COHORTS <- list(
  hypertensive = list(
    label  = "Hypertensive (§2.2 primary)",
    filter = function(d) filter(d, BASELINE_AHA %in% AHA_HYPERTENSIVE),
    strata = AHA_HYPERTENSIVE
  ),
  full = list(
    label  = "Full anchored (§2.1, incl. Normal)",
    filter = function(d) filter(d, BASELINE_AHA %in% AHA_FULL),
    strata = AHA_FULL
  ),
  htn_naive = list(
    label  = "HTN-naive (hypertensive + never_logged)",
    filter = function(d) filter(d, BASELINE_AHA %in% AHA_HYPERTENSIVE,
                                   HTN_MED_BASELINE == "never_logged"),
    strata = AHA_HYPERTENSIVE
  )
)

apply_cohort <- function(d, name) {
  if (!name %in% names(COHORTS)) stop("Unknown cohort: ", name)
  COHORTS[[name]]$filter(d)
}

# --- Cell formatters (carried over from the now-retired cohort_characteristics.R) ----
msd <- function(x, d = 1) {
  sprintf("%.*f ± %.*f", d, mean(x, na.rm = TRUE), d, sd(x, na.rm = TRUE))
}
n_pct <- function(n, total, d = 1) {
  sprintf("%s (%.*f)", format(n, big.mark = ","), d, 100 * n / total)
}

# --- HTN-med 3-level collapse (matches the reported rows; tests must match) ---
htn3 <- function(x) {
  case_when(
    x %in% c("stable_on_meds", "recent_initiator") ~ "On at baseline",
    x == "mid_study_initiator"                     ~ "Initiated during study",
    x == "never_logged"                            ~ "Never logged",
    TRUE                                           ~ NA_character_
  )
}

# --- Between-stratum tests ----------------------------------------------------
# Continuous: one-way ANOVA F-test (classic; pairs with mean +/- SD reporting).
p_continuous <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  if (length(unique(g[ok])) < 2) return(NA_real_)
  summary(aov(x[ok] ~ factor(g[ok])))[[1L]][["Pr(>F)"]][1L]
}

# Categorical: chi-square, auto-falling back to Fisher's exact when any expected
# cell count < 5. For the current cohorts every expected count is in the
# hundreds, so chi-square is always used; the fallback is defensive for any
# future sparse sub-cohort. set.seed() makes the simulated-p path deterministic
# if it is ever exercised (keeps the pipeline reproducible).
p_categorical <- function(f, g) {
  ok  <- !is.na(f) & !is.na(g)
  tab <- table(f[ok], g[ok])
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (any(dim(tab) < 2)) return(NA_real_)
  expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  if (any(expected < 5)) {
    set.seed(1)
    fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value
  } else {
    suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
  }
}

fmt_p <- function(p) {
  ifelse(is.na(p), "—",
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

# --- Per-stratum characteristic row ------------------------------------------
compute_one_stratum <- function(d) {
  n_total       <- nrow(d)
  male_n        <- sum(d$SEX == "male",   na.rm = TRUE)
  female_n      <- sum(d$SEX == "female", na.rm = TRUE)
  on_baseline_n <- sum(d$HTN_MED_BASELINE %in% c("stable_on_meds", "recent_initiator"), na.rm = TRUE)
  mid_study_n   <- sum(d$HTN_MED_BASELINE == "mid_study_initiator", na.rm = TRUE)
  never_log_n   <- sum(d$HTN_MED_BASELINE == "never_logged",        na.rm = TRUE)

  tibble(
    `N`                                       = format(n_total, big.mark = ","),
    `Age, mean ± SD, y`                  = msd(d$AGE_AT_START),
    `Sex: Male, n (%)`                        = n_pct(male_n,   n_total),
    `Sex: Female, n (%)`                      = n_pct(female_n, n_total),
    `BMI, mean ± SD, kg/m²`         = msd(d$BASELINE_BMI),
    `Baseline SBP, mean ± SD, mmHg`      = msd(d$BASELINE_SBP),
    `Baseline DBP, mean ± SD, mmHg`      = msd(d$BASELINE_DBP),
    `HTN meds: On at baseline, n (%)`         = n_pct(on_baseline_n, n_total),
    `HTN meds: Initiated during study, n (%)` = n_pct(mid_study_n,   n_total),
    `HTN meds: Never logged, n (%)`           = n_pct(never_log_n,   n_total)
  )
}

# --- Between-stratum p-value per characteristic group ------------------------
# Sex test = Male vs Female only (non-binary / prefer-not-to-say excluded from
# the sex comparison, consistent with the Table 1 sex-row footnote). HTN test
# uses the 3-level collapse so it matches the reported categories.
stratum_pvals <- function(d) {
  sex_ok <- d %>% filter(SEX %in% c("male", "female"))
  c(
    age = p_continuous(d$AGE_AT_START, d$BASELINE_AHA),
    sex = p_categorical(sex_ok$SEX, sex_ok$BASELINE_AHA),
    bmi = p_continuous(d$BASELINE_BMI, d$BASELINE_AHA),
    sbp = p_continuous(d$BASELINE_SBP, d$BASELINE_AHA),
    dbp = p_continuous(d$BASELINE_DBP, d$BASELINE_AHA),
    htn = p_categorical(htn3(d$HTN_MED_BASELINE), d$BASELINE_AHA)
  )
}

# --- Display label for stratum columns (drop "Hypertension" suffix) ----------
stratum_label <- function(s) {
  case_when(
    s == "Stage 1 Hypertension" ~ "Stage 1",
    s == "Stage 2 Hypertension" ~ "Stage 2",
    TRUE                        ~ s
  )
}

# --- Characteristic table builder: Characteristic x {Overall, strata, P} -----
make_char_table <- function(d, strata_levels) {
  overall <- compute_one_stratum(d) %>% mutate(stratum = "Overall", .before = 1)

  per_stratum <- d %>%
    mutate(BASELINE_AHA = factor(BASELINE_AHA, levels = strata_levels)) %>%
    group_by(BASELINE_AHA) %>%
    group_modify(~ compute_one_stratum(.x)) %>%
    ungroup() %>%
    rename(stratum = BASELINE_AHA) %>%
    mutate(stratum = as.character(stratum))

  wide <- bind_rows(overall, per_stratum) %>%
    mutate(stratum = stratum_label(stratum)) %>%
    pivot_longer(-stratum, names_to = "Characteristic", values_to = "value") %>%
    pivot_wider(names_from = stratum, values_from = value)

  # Attach the between-stratum p-value to the primary row of each variable
  # (multi-row variables Sex / HTN meds carry their p on the first sub-row).
  pv <- stratum_pvals(d)
  wide$`P-value` <- case_when(
    wide$Characteristic == "Age, mean ± SD, y"             ~ fmt_p(pv[["age"]]),
    wide$Characteristic == "Sex: Male, n (%)"                  ~ fmt_p(pv[["sex"]]),
    wide$Characteristic == "BMI, mean ± SD, kg/m²"   ~ fmt_p(pv[["bmi"]]),
    wide$Characteristic == "Baseline SBP, mean ± SD, mmHg" ~ fmt_p(pv[["sbp"]]),
    wide$Characteristic == "Baseline DBP, mean ± SD, mmHg" ~ fmt_p(pv[["dbp"]]),
    wide$Characteristic == "HTN meds: On at baseline, n (%)"   ~ fmt_p(pv[["htn"]]),
    TRUE                                                       ~ ""
  )

  wide
}

# --- Table A size-emitter contract -------------------------------------------
# Each modeled chapter computes its own size table FROM its fitted objects
# (approach A2: numbers read off the actual models, never hand-typed) and calls
# emit_sizes() to write a uniform <chapter>_sizes_std.csv. The aggregator
# (0_cohort/03_modeled_n_supp.R) globs these into Table A. Standard schema:
#   section       SAP/chapter section, e.g. "5.1"
#   role          "primary" or "sensitivity" — routes the row to Table A vs A.2
#   analysis      short label, e.g. "Primary BP-trajectory LMM"
#   cohort        hypertensive / full / htn_naive / ...
#   outcome       SBP / DBP / SBP responder / Δ HbA1c / ...
#   predictor_set model/variant label, e.g. "unadjusted" / "adjusted"
#   n_users       distinct users in the fitted model
#   n_obs         model rows (= n_users cross-sectionally; user×timepoint for LMMs)
#   n_events      events for logistic outcomes (NA for continuous)
#   ceiling_n     cohort entry N (max possible) — the Table B denominator;
#                 n_users / ceiling_n = % of cohort retained
SIZES_SCHEMA <- c("section", "role", "analysis", "cohort", "outcome", "predictor_set",
                  "n_users", "n_obs", "n_events", "ceiling_n")

emit_sizes <- function(df, chapter, results_dir) {
  missing_cols <- setdiff(SIZES_SCHEMA, names(df))
  if (length(missing_cols) > 0)
    stop("emit_sizes: missing required column(s): ", paste(missing_cols, collapse = ", "))
  df <- df[, SIZES_SCHEMA]   # enforce column order; drop any extras
  path <- file.path(results_dir, paste0(chapter, "_sizes_std.csv"))
  readr::write_csv(df, path)
  message("emit_sizes: wrote ", path, " (", nrow(df), " rows)")
  invisible(path)
}

# --- Shared §5.1-family covariate + contrast helpers -------------------------
# Baseline-HTN-status covariate (baseline-only trait; mid-study initiation
# collapsed to "not on" to avoid mediator-adjustment). Restores the v0.5-locked
# covariate per CHANGELOG v0.79. (The §5.1 primary script carries an identical
# inline definition; these helpers standardize it for the sensitivity scripts.)
htn_baseline_on <- function(x) {
  factor(if_else(x %in% c("stable_on_meds", "recent_initiator"),
                 "on_at_baseline", "not_on_at_baseline"),
         levels = c("not_on_at_baseline", "on_at_baseline"))
}

# Tidy an emmeans contrast set: attach 95% CIs (infer) + within-family
# Holm-adjusted p. Shared so the §5.1 universal-Holm convention is identical
# across the primary and every sensitivity.
tidy_contrast <- function(cobj) {
  df <- as_tibble(summary(cobj, infer = c(TRUE, TRUE), adjust = "none"))
  df %>% mutate(p.holm = p.adjust(p.value, method = "holm"))
}
