# 98_qc_user_frame.R
# Phase 1 / Slice 98: centralized QC pass over the 8 substrate files.
# Applies per-variable plausibility bounds (NA-outs implausible values) and
# re-derives any downstream values that depend on bounded raw inputs.
#
# Inputs:  Build/derived/*.rds (slices 01–08, raw substrate outputs)
# Outputs: Build/derived/qc/*.rds + .csv (QC'd substrate outputs, mirroring filenames)
#
# Architecture: raw substrates in derived/ remain unchanged (auditable history);
# QC'd substrates in derived/qc/ are the inputs to slice 99 (final assembly).

library(tidyverse)

project_dir <- "./PRIMARY"
derived_dir <- file.path(project_dir, "Build", "derived")
qc_dir      <- file.path(derived_dir, "qc")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

# --- QC helper: NA-out values outside [lo, hi]; report counts. ---
# Used both directly (qc_bound(x, lo, hi, "NAME")) and inside across(... ~ qc_bound(.x, lo, hi, cur_column())).
qc_bound <- function(x, lo, hi, name = "") {
  bad <- !is.na(x) & (x < lo | x > hi)
  n_bad <- sum(bad)
  if (n_bad > 0) {
    cat(sprintf("  %-40s  NA'd %5d values outside [%g, %g]\n", name, n_bad, lo, hi))
  }
  if_else(!is.na(x) & x >= lo & x <= hi, x, NA_real_)
}

# --- Unit-conversion helper for lab analytes with documented mg/dL vs alternate-unit ambiguity. ---
# Applied BEFORE qc_bound; converts plausibly-alternate-unit values to the target unit (mg/dL or
# NGSP %) and leaves the rest unchanged. The subsequent qc_bound then NA's anything still outside
# the unit-correct plausible range (catches real outliers + any unit-mismatched values our heuristic
# didn't catch).
#
# direction:
#   "divide_above":    if x > threshold, divide by factor  (e.g., creatinine μmol/L → mg/dL: x/88.4)
#   "multiply_below":  if x < threshold, multiply by factor (e.g., HDL mmol/L → mg/dL: x*38.67)
#   "hba1c_above":     if x > threshold, apply IFCC mmol/mol → NGSP % conversion (x/10.929 + 2.15)
#
# Documented surfacing 2026-05-22: ~78% of creatinine NA'd values were in plausible μmol/L range,
# despite database labeling as mg/dL. Similar mixed-unit risk applies to all 6 SAP labs +
# Vitamin D + ApoB (each has a well-documented alternate unit). eGFR has no unit ambiguity.
convert_units <- function(x, threshold, factor, direction, name = "") {
  if (direction == "divide_above") {
    mask <- !is.na(x) & x > threshold
    x_new <- if_else(mask, x / factor, x)
  } else if (direction == "multiply_below") {
    mask <- !is.na(x) & x < threshold
    x_new <- if_else(mask, x * factor, x)
  } else if (direction == "hba1c_above") {
    mask <- !is.na(x) & x > threshold
    x_new <- if_else(mask, (x / 10.929) + 2.15, x)
  } else {
    stop("Unknown direction: ", direction)
  }
  n_converted <- sum(mask)
  if (n_converted > 0) {
    cat(sprintf("  %-40s  Converted %5d values (%s, threshold=%g)\n",
                name, n_converted, direction, threshold))
  }
  x_new
}

# --- Read raw substrates ---
bp_timepoints <- read_rds(file.path(derived_dir, "01_bp_timepoints.rds"))
bp_outcomes   <- read_rds(file.path(derived_dir, "02_bp_outcomes.rds"))
demo_clinical <- read_rds(file.path(derived_dir, "03_demographics_clinical.rds"))
lab_weight    <- read_rds(file.path(derived_dir, "04_lab_weight_outcomes.rds"))
whoop_means   <- read_rds(file.path(derived_dir, "05_whoop_window_means.rds"))
le8_toba_nic  <- read_rds(file.path(derived_dir, "06_le8_tobacco_nicotine.rds"))
engagement    <- read_rds(file.path(derived_dir, "07_engagement.rds"))
lmm_cov       <- read_rds(file.path(derived_dir, "08_lmm_covariates.rds"))

cat("Loaded substrates.\n")

# =====================================================================
# Slice 01 — BP value bounds
# =====================================================================
cat("\n=== Slice 01 (bp_timepoints) — BP value bounds ===\n")
bp_timepoints_qc <- bp_timepoints %>%
  mutate(
    across(c(BASELINE_SBP, MO3_SBP, MO6_SBP, MO12_SBP),
           ~ qc_bound(.x, 70, 250, cur_column())),
    across(c(BASELINE_DBP, MO3_DBP, MO6_DBP, MO12_DBP),
           ~ qc_bound(.x, 40, 150, cur_column()))
  )

# =====================================================================
# Slice 02 — re-derive from QC'd slice 01 (AHA, Δ, responder, BPCAT)
# =====================================================================
cat("\n=== Slice 02 (bp_outcomes) — re-derive from QC'd slice 01 ===\n")

classify_aha <- function(sbp, dbp) {
  case_when(
    is.na(sbp) | is.na(dbp) ~ NA_character_,
    sbp >= 140 | dbp >= 90  ~ "Stage 2 Hypertension",
    sbp >= 130 | dbp >= 80  ~ "Stage 1 Hypertension",
    sbp >= 120 & dbp <  80  ~ "Elevated",
    TRUE                    ~ "Normal"
  )
}
bpcat_levels <- c("improved", "stable", "worsened")
bin_bp_change <- function(delta) {
  factor(
    case_when(
      is.na(delta) ~ NA_character_,
      delta <= -5  ~ "improved",
      delta >=  5  ~ "worsened",
      TRUE         ~ "stable"
    ),
    levels = bpcat_levels
  )
}

bp_outcomes_qc <- bp_timepoints_qc %>%
  mutate(
    BASELINE_AHA   = classify_aha(BASELINE_SBP, BASELINE_DBP),
    MO3_AHA        = classify_aha(MO3_SBP,      MO3_DBP),
    MO6_AHA        = classify_aha(MO6_SBP,      MO6_DBP),
    MO12_AHA       = classify_aha(MO12_SBP,     MO12_DBP),
    DELTA_3MO_SBP  = MO3_SBP  - BASELINE_SBP,
    DELTA_3MO_DBP  = MO3_DBP  - BASELINE_DBP,
    DELTA_6MO_SBP  = MO6_SBP  - BASELINE_SBP,
    DELTA_6MO_DBP  = MO6_DBP  - BASELINE_DBP,
    DELTA_12MO_SBP = MO12_SBP - BASELINE_SBP,
    DELTA_12MO_DBP = MO12_DBP - BASELINE_DBP,
    SBP_RESPONDER_6MO = DELTA_6MO_SBP <= -5,
    BPCAT_6MO_SBP  = bin_bp_change(DELTA_6MO_SBP),
    BPCAT_6MO_DBP  = bin_bp_change(DELTA_6MO_DBP),
    BPCAT_12MO_SBP = bin_bp_change(DELTA_12MO_SBP),
    BPCAT_12MO_DBP = bin_bp_change(DELTA_12MO_DBP)
  ) %>%
  select(STUDY_ID,
         BASELINE_AHA, MO3_AHA, MO6_AHA, MO12_AHA,
         DELTA_3MO_SBP,  DELTA_3MO_DBP,
         DELTA_6MO_SBP,  DELTA_6MO_DBP,
         DELTA_12MO_SBP, DELTA_12MO_DBP,
         SBP_RESPONDER_6MO,
         BPCAT_6MO_SBP,  BPCAT_6MO_DBP,
         BPCAT_12MO_SBP, BPCAT_12MO_DBP)

# Cross-check: per-stratum n's after QC vs original
cat("BASELINE_AHA counts — original vs QC'd:\n")
print(full_join(
  bp_outcomes   %>% count(BASELINE_AHA) %>% rename(n_raw = n),
  bp_outcomes_qc %>% count(BASELINE_AHA) %>% rename(n_qc = n),
  by = "BASELINE_AHA"
))

# =====================================================================
# Slice 03 — passthrough (slice 3 anthropometrics already bounded in source)
# =====================================================================
cat("\n=== Slice 03 (demo_clinical) — passthrough ===\n")
demo_clinical_qc <- demo_clinical

# =====================================================================
# Slice 04 — unit conversion (heuristic) + bounds + Δ re-derivation
# =====================================================================
cat("\n=== Slice 04 (lab_weight_outcomes) — unit conversion + bounds + Δ re-derivation ===\n")

cat("\n  Step 1: heuristic unit conversion to mg/dL (NGSP % for HbA1c)\n")
cat("  Values clearly outside plausible mg/dL range are assumed to be in the alternate documented\n")
cat("  unit and converted; subsequent qc_bound NA's any residual outliers.\n")

lab_weight_qc <- lab_weight %>%
  mutate(
    # Creatinine: μmol/L → mg/dL (÷ 88.4) if value > 20 (mg/dL ≤ 15 plausible; >20 = μmol/L)
    PRE_CREATININE         = convert_units(PRE_CREATININE,          20,  88.40, "divide_above",   "PRE_CREATININE"),
    POST_CREATININE        = convert_units(POST_CREATININE,         20,  88.40, "divide_above",   "POST_CREATININE"),
    # Total cholesterol: mmol/L → mg/dL (× 38.67) if value < 20 (mmol/L max ~13; mg/dL min ~50)
    PRE_TOTAL_CHOLESTEROL  = convert_units(PRE_TOTAL_CHOLESTEROL,   20,  38.67, "multiply_below", "PRE_TOTAL_CHOLESTEROL"),
    POST_TOTAL_CHOLESTEROL = convert_units(POST_TOTAL_CHOLESTEROL,  20,  38.67, "multiply_below", "POST_TOTAL_CHOLESTEROL"),
    # HDL: mmol/L → mg/dL (× 38.67) if value < 10
    PRE_HDL_CHOLESTEROL    = convert_units(PRE_HDL_CHOLESTEROL,     10,  38.67, "multiply_below", "PRE_HDL_CHOLESTEROL"),
    POST_HDL_CHOLESTEROL   = convert_units(POST_HDL_CHOLESTEROL,    10,  38.67, "multiply_below", "POST_HDL_CHOLESTEROL"),
    # LDL: mmol/L → mg/dL (× 38.67) if value < 10
    PRE_LDL_CHOLESTEROL    = convert_units(PRE_LDL_CHOLESTEROL,     10,  38.67, "multiply_below", "PRE_LDL_CHOLESTEROL"),
    POST_LDL_CHOLESTEROL   = convert_units(POST_LDL_CHOLESTEROL,    10,  38.67, "multiply_below", "POST_LDL_CHOLESTEROL"),
    # Triglycerides: mmol/L → mg/dL (× 88.6) if value < 20
    PRE_TRIGLYCERIDES      = convert_units(PRE_TRIGLYCERIDES,       20,  88.60, "multiply_below", "PRE_TRIGLYCERIDES"),
    POST_TRIGLYCERIDES     = convert_units(POST_TRIGLYCERIDES,      20,  88.60, "multiply_below", "POST_TRIGLYCERIDES"),
    # Glucose: mmol/L → mg/dL (× 18.02) if value < 30
    PRE_GLUCOSE            = convert_units(PRE_GLUCOSE,             30,  18.02, "multiply_below", "PRE_GLUCOSE"),
    POST_GLUCOSE           = convert_units(POST_GLUCOSE,            30,  18.02, "multiply_below", "POST_GLUCOSE"),
    # HbA1c: IFCC mmol/mol → NGSP % via formula if value > 20
    PRE_HBA1C              = convert_units(PRE_HBA1C,               20,  NA,    "hba1c_above",    "PRE_HBA1C"),
    POST_HBA1C             = convert_units(POST_HBA1C,              20,  NA,    "hba1c_above",    "POST_HBA1C"),
    # Vitamin D: nmol/L → ng/mL (÷ 2.5) if value > 200
    PRE_VITAMIN_D          = convert_units(PRE_VITAMIN_D,          200,   2.5,  "divide_above",   "PRE_VITAMIN_D"),
    POST_VITAMIN_D         = convert_units(POST_VITAMIN_D,         200,   2.5,  "divide_above",   "POST_VITAMIN_D"),
    # ApoB: g/L → mg/dL (× 100) if value < 5
    PRE_APOLIPOPROTEIN_B   = convert_units(PRE_APOLIPOPROTEIN_B,    5,  100,    "multiply_below", "PRE_APOLIPOPROTEIN_B"),
    POST_APOLIPOPROTEIN_B  = convert_units(POST_APOLIPOPROTEIN_B,   5,  100,    "multiply_below", "POST_APOLIPOPROTEIN_B")
    # eGFR: no unit ambiguity (always mL/min/1.73m²) — no conversion step.
  ) %>%
  # Step 2: bounds (after conversion)
  mutate(
    WEIGHT_BASELINE_MEASURED_KG = qc_bound(WEIGHT_BASELINE_MEASURED_KG, 30, 250, "WEIGHT_BASELINE_MEASURED_KG"),
    WEIGHT_FOLLOWUP_MEASURED_KG = qc_bound(WEIGHT_FOLLOWUP_MEASURED_KG, 30, 250, "WEIGHT_FOLLOWUP_MEASURED_KG"),
    PRE_TOTAL_CHOLESTEROL  = qc_bound(PRE_TOTAL_CHOLESTEROL,   50,  500, "PRE_TOTAL_CHOLESTEROL"),
    POST_TOTAL_CHOLESTEROL = qc_bound(POST_TOTAL_CHOLESTEROL,  50,  500, "POST_TOTAL_CHOLESTEROL"),
    PRE_HDL_CHOLESTEROL    = qc_bound(PRE_HDL_CHOLESTEROL,     10,  150, "PRE_HDL_CHOLESTEROL"),
    POST_HDL_CHOLESTEROL   = qc_bound(POST_HDL_CHOLESTEROL,    10,  150, "POST_HDL_CHOLESTEROL"),
    PRE_LDL_CHOLESTEROL    = qc_bound(PRE_LDL_CHOLESTEROL,     10,  400, "PRE_LDL_CHOLESTEROL"),
    POST_LDL_CHOLESTEROL   = qc_bound(POST_LDL_CHOLESTEROL,    10,  400, "POST_LDL_CHOLESTEROL"),
    PRE_TRIGLYCERIDES      = qc_bound(PRE_TRIGLYCERIDES,       20, 2000, "PRE_TRIGLYCERIDES"),
    POST_TRIGLYCERIDES     = qc_bound(POST_TRIGLYCERIDES,      20, 2000, "POST_TRIGLYCERIDES"),
    PRE_APOLIPOPROTEIN_B   = qc_bound(PRE_APOLIPOPROTEIN_B,    20,  250, "PRE_APOLIPOPROTEIN_B"),
    POST_APOLIPOPROTEIN_B  = qc_bound(POST_APOLIPOPROTEIN_B,   20,  250, "POST_APOLIPOPROTEIN_B"),
    PRE_GLUCOSE            = qc_bound(PRE_GLUCOSE,             30,  600, "PRE_GLUCOSE"),
    POST_GLUCOSE           = qc_bound(POST_GLUCOSE,            30,  600, "POST_GLUCOSE"),
    PRE_HBA1C              = qc_bound(PRE_HBA1C,                3,   20, "PRE_HBA1C"),
    POST_HBA1C             = qc_bound(POST_HBA1C,               3,   20, "POST_HBA1C"),
    PRE_VITAMIN_D          = qc_bound(PRE_VITAMIN_D,            3,  200, "PRE_VITAMIN_D"),
    POST_VITAMIN_D         = qc_bound(POST_VITAMIN_D,           3,  200, "POST_VITAMIN_D"),
    PRE_CREATININE         = qc_bound(PRE_CREATININE,         0.2,   15, "PRE_CREATININE"),
    POST_CREATININE        = qc_bound(POST_CREATININE,        0.2,   15, "POST_CREATININE"),
    PRE_EGFR               = qc_bound(PRE_EGFR,                 5,  200, "PRE_EGFR"),
    POST_EGFR              = qc_bound(POST_EGFR,                5,  200, "POST_EGFR")
  ) %>%
  # Re-derive non-HDL from QC'd PRE/POST TC and HDL (overwrites slice 04's raw
  # non-HDL, which was computed before unit conversion)
  mutate(
    PRE_NON_HDL_CHOLESTEROL  = PRE_TOTAL_CHOLESTEROL  - PRE_HDL_CHOLESTEROL,
    POST_NON_HDL_CHOLESTEROL = POST_TOTAL_CHOLESTEROL - POST_HDL_CHOLESTEROL
  ) %>%
  # Bound non-HDL post-derivation — catches negative values (HDL > TC, impossible
  # by definition) and implausibly high values
  mutate(
    PRE_NON_HDL_CHOLESTEROL  = qc_bound(PRE_NON_HDL_CHOLESTEROL,  10, 450, "PRE_NON_HDL_CHOLESTEROL (derived)"),
    POST_NON_HDL_CHOLESTEROL = qc_bound(POST_NON_HDL_CHOLESTEROL, 10, 450, "POST_NON_HDL_CHOLESTEROL (derived)")
  ) %>%
  # Re-derive Δs from QC'd PRE/POST (incl. derived non-HDL)
  mutate(
    DELTA_WEIGHT_KG           = WEIGHT_FOLLOWUP_MEASURED_KG - WEIGHT_BASELINE_MEASURED_KG,
    DELTA_TOTAL_CHOLESTEROL   = POST_TOTAL_CHOLESTEROL   - PRE_TOTAL_CHOLESTEROL,
    DELTA_HDL_CHOLESTEROL     = POST_HDL_CHOLESTEROL     - PRE_HDL_CHOLESTEROL,
    DELTA_LDL_CHOLESTEROL     = POST_LDL_CHOLESTEROL     - PRE_LDL_CHOLESTEROL,
    DELTA_NON_HDL_CHOLESTEROL = POST_NON_HDL_CHOLESTEROL - PRE_NON_HDL_CHOLESTEROL,
    DELTA_TRIGLYCERIDES       = POST_TRIGLYCERIDES       - PRE_TRIGLYCERIDES,
    DELTA_GLUCOSE             = POST_GLUCOSE             - PRE_GLUCOSE,
    DELTA_HBA1C               = POST_HBA1C               - PRE_HBA1C
  ) %>%
  # Final bound on Δweight
  mutate(
    DELTA_WEIGHT_KG = qc_bound(DELTA_WEIGHT_KG, -50, 50, "DELTA_WEIGHT_KG (post-arithmetic)")
  )

# =====================================================================
# Slice 05 — WHOOP window-mean bounds
# =====================================================================
cat("\n=== Slice 05 (whoop_means) ===\n")
whoop_means_qc <- whoop_means %>%
  mutate(
    across(starts_with("RHR_")               & ends_with("_MEAN"), ~ qc_bound(.x,  30, 130, cur_column())),
    across(starts_with("HRV_")               & ends_with("_MEAN"), ~ qc_bound(.x,   5, 300, cur_column())),
    across(starts_with("STRAIN_")            & ends_with("_MEAN"), ~ qc_bound(.x,   0,  21, cur_column())),
    across(starts_with("SLEEP_DURATION_")    & ends_with("_MEAN"), ~ qc_bound(.x,   0,  15, cur_column())),
    across(starts_with("SLEEP_CONSISTENCY_") & ends_with("_MEAN"), ~ qc_bound(.x,   0, 100, cur_column()))
  )

# =====================================================================
# Slice 06 — passthrough (counts bounded by window construction)
# =====================================================================
cat("\n=== Slice 06 (le8_tobacco_nicotine) — passthrough ===\n")
le8_toba_nic_qc <- le8_toba_nic

# =====================================================================
# Slice 07 — engagement bounds; re-derive logs
# =====================================================================
cat("\n=== Slice 07 (engagement) — bounds + log re-derivation ===\n")
engagement_qc <- engagement %>%
  mutate(
    PCT_DAYS_JOURNALED             = qc_bound(PCT_DAYS_JOURNALED,             0, 102, "PCT_DAYS_JOURNALED"),
    PCT_DAYS_ENGAGED_NOTIFICATIONS = qc_bound(PCT_DAYS_ENGAGED_NOTIFICATIONS, 0, 102, "PCT_DAYS_ENGAGED_NOTIFICATIONS"),
    PCT_DAYS_USED_COACH            = qc_bound(PCT_DAYS_USED_COACH,            0, 102, "PCT_DAYS_USED_COACH"),
    PCT_DAYS_MET_OPTIMAL_STRAIN    = qc_bound(PCT_DAYS_MET_OPTIMAL_STRAIN,    0, 102, "PCT_DAYS_MET_OPTIMAL_STRAIN"),
    SLEEP_PCT_OF_TARGET            = qc_bound(SLEEP_PCT_OF_TARGET,            0, 200, "SLEEP_PCT_OF_TARGET"),
    AVG_SLEEP_NEED_PCT             = qc_bound(AVG_SLEEP_NEED_PCT,             0, 100, "AVG_SLEEP_NEED_PCT"),
    AVG_DAILY_TIME_IN_APP_MIN      = qc_bound(AVG_DAILY_TIME_IN_APP_MIN,      0, 240, "AVG_DAILY_TIME_IN_APP_MIN"),
    N_BP_READINGS_STUDY            = qc_bound(N_BP_READINGS_STUDY,            0, 543, "N_BP_READINGS_STUDY")
  ) %>%
  # Re-derive logs from QC'd raw values (NA propagates)
  mutate(
    LOG_AVG_DAILY_TIME_IN_APP_MIN = log(AVG_DAILY_TIME_IN_APP_MIN + 1),
    LOG_N_BP_READINGS_STUDY       = log(N_BP_READINGS_STUDY + 1)
  )

# =====================================================================
# Slice 08 — passthrough
# =====================================================================
cat("\n=== Slice 08 (lmm_covariates) — passthrough ===\n")
lmm_cov_qc <- lmm_cov

# =====================================================================
# Write QC'd substrates (rds + csv)
# =====================================================================
out <- list(
  "01_bp_timepoints"          = bp_timepoints_qc,
  "02_bp_outcomes"            = bp_outcomes_qc,
  "03_demographics_clinical"  = demo_clinical_qc,
  "04_lab_weight_outcomes"    = lab_weight_qc,
  "05_whoop_window_means"     = whoop_means_qc,
  "06_le8_tobacco_nicotine"   = le8_toba_nic_qc,
  "07_engagement"             = engagement_qc,
  "08_lmm_covariates"         = lmm_cov_qc
)

cat("\n=== Writing QC'd substrate files to derived/qc/ ===\n")
for (name in names(out)) {
  rds_path <- file.path(qc_dir, paste0(name, ".rds"))
  csv_path <- file.path(qc_dir, paste0(name, ".csv"))
  write_rds(out[[name]], rds_path)
  write_csv(out[[name]], csv_path)
  cat(sprintf("  wrote %s (%d rows × %d cols)\n", name, nrow(out[[name]]), ncol(out[[name]])))
}

cat("\nDONE. QC'd substrates ready as inputs for slice 99 (final assembly).\n")
