# 01_bp_timepoints.R
# Derive per-user BP timepoint readings (baseline / 3mo / 6mo / 12mo) from
# raw BP_MEASURES per SAP §3, and validate against the roster's pre-computed
# BASELINE_*, FOLLOWUP_*, FOLLOWUP_MIDPOINT_DAY, DELTA_* columns.
#
# Phase 1 / Slice 1 of the user-level frame build (see Build/BUILD_TODO.md).
#
# SAP §3 windowing (uniform tiebreak: closest to target day; earlier reading_ts):
#   Baseline: target 0,   single first manual-cuff reading on STUDY_START_DATE
#   3 months: target 90,  window 60–120 d   (±30 d)
#   6 months: target 180, window 150–240 d  (−30 / +60 d)   PRIMARY follow-up
#   12 months: target 365, window 335–390 d (−30 / +25 d)   SAP §5.5 secondary
#
# BP_MEASURES is already manual-cuff only (confirmed 2026-05-21); no source-type
# filter needed. PERIOD column is decoration — filter on DAY_FROM_START only.

library(tidyverse)

project_dir <- "./PRIMARY"
data_dir    <- file.path(project_dir, "Query Results")
roster_dir  <- file.path(project_dir, "Roster")
build_dir   <- file.path(project_dir, "Build")
derived_dir <- file.path(build_dir, "derived")
pull_dir    <- file.path(data_dir, "2026-05-12 Pulls")

dir.create(derived_dir, showWarnings = FALSE, recursive = TRUE)

# Tolerance for BP-value equality checks. Roster columns and BP_MEASURES columns
# both round-trip through the warehouse → CSV → R at finite precision; arithmetic on
# the R side (FOLLOWUP - BASELINE) can land on a double a few bits away from the
# directly-read DELTA. 1e-4 mmHg is well below clinical precision (1 mmHg) and
# orders of magnitude above any plausible CSV roundtrip noise.
mmhg_tol <- 1e-4

# --- Read ---
roster <- read_csv(file.path(roster_dir, "ROSTER_1X_2026-05-12.csv"))
bp     <- read_csv(file.path(pull_dir, "BP_MEASURES_1X_2026-05-12.csv"))

# --- Pre-flight ---
cat("Roster rows:", nrow(roster), " users:", n_distinct(roster$STUDY_ID), "\n")
cat("BP_MEASURES rows:", nrow(bp), " users:", n_distinct(bp$STUDY_ID), "\n")
cat("STUDY_ID type — roster:", class(roster$STUDY_ID),
    " bp:", class(bp$STUDY_ID), "\n")
cat("STUDY_START_DATE class:", class(roster$STUDY_START_DATE), "\n")
cat("DATE_TIME class:", class(bp$DATE_TIME), "\n")

stopifnot(!any(duplicated(roster$STUDY_ID)))

bp_not_in_roster <- setdiff(unique(bp$STUDY_ID), roster$STUDY_ID)
roster_not_in_bp <- setdiff(roster$STUDY_ID, unique(bp$STUDY_ID))
cat("Users in BP_MEASURES but not roster:", length(bp_not_in_roster), "\n")
cat("Users in roster but no BP_MEASURES rows:", length(roster_not_in_bp), "\n")

# --- Join STUDY_START_DATE and compute day-offset ---
bp <- bp %>%
  inner_join(roster %>% select(STUDY_ID, STUDY_START_DATE), by = "STUDY_ID") %>%
  mutate(
    READING_DATE   = as.Date(DATE_TIME),
    DAY_FROM_START = as.integer(READING_DATE - STUDY_START_DATE)
  )

cat("BP rows after roster-scoped join:", nrow(bp), "\n")
cat("DAY_FROM_START range: [", min(bp$DAY_FROM_START), ",",
    max(bp$DAY_FROM_START), "]\n")

# --- Reusable timepoint derivation (SAP §3) ---
# Window filter → arrange by |day − target| asc, DATE_TIME asc → first row per user.
# Returns one row per user with the chosen reading's SBP, DBP, DATE_TIME, day-from-start.
derive_timepoint <- function(bp, target_day, window_lo, window_hi) {
  bp %>%
    filter(DAY_FROM_START >= window_lo, DAY_FROM_START <= window_hi) %>%
    group_by(STUDY_ID) %>%
    arrange(abs(DAY_FROM_START - target_day), DATE_TIME, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(STUDY_ID, SBP, DBP, DATE_TIME, DAY_FROM_START)
}

baseline_r <- derive_timepoint(bp, target_day = 0,   window_lo = 0,   window_hi = 0)
mo3_r      <- derive_timepoint(bp, target_day = 90,  window_lo = 60,  window_hi = 120)
mo6_r      <- derive_timepoint(bp, target_day = 180, window_lo = 150, window_hi = 240)
mo12_r     <- derive_timepoint(bp, target_day = 365, window_lo = 335, window_hi = 390)

cat("\nR-derived timepoint coverage (users with a qualifying reading):\n")
cat("  Baseline (day 0):     ", nrow(baseline_r), "\n")
cat("  3mo (60–120):         ", nrow(mo3_r),      "\n")
cat("  6mo (150–240):        ", nrow(mo6_r),      "\n")
cat("  12mo (335–390):       ", nrow(mo12_r),     "\n")

# --- Validation: baseline vs. roster BASELINE_SBP / BASELINE_DBP ---
check_base <- roster %>%
  select(STUDY_ID, BASELINE_SBP, BASELINE_DBP) %>%
  left_join(
    baseline_r %>% rename(BASELINE_SBP_R = SBP, BASELINE_DBP_R = DBP,
                          BASELINE_READING_TS = DATE_TIME,
                          BASELINE_DAY = DAY_FROM_START),
    by = "STUDY_ID"
  )

n_base_sbp_disagree <- sum(!near(check_base$BASELINE_SBP, check_base$BASELINE_SBP_R, tol = mmhg_tol), na.rm = TRUE)
n_base_dbp_disagree <- sum(!near(check_base$BASELINE_DBP, check_base$BASELINE_DBP_R, tol = mmhg_tol), na.rm = TRUE)
n_base_missing      <- sum(is.na(check_base$BASELINE_SBP_R))
cat("\nBaseline validation vs. roster:\n")
cat("  SBP disagreements:", n_base_sbp_disagree, "\n")
cat("  DBP disagreements:", n_base_dbp_disagree, "\n")
cat("  Users with no R-derived baseline:", n_base_missing, "\n")
if (n_base_sbp_disagree > 0 || n_base_dbp_disagree > 0) {
  cat("  First few disagreements:\n")
  print(head(check_base %>%
               filter(!near(BASELINE_SBP, BASELINE_SBP_R, tol = mmhg_tol) |
                      !near(BASELINE_DBP, BASELINE_DBP_R, tol = mmhg_tol)), 10))
}
stopifnot(n_base_sbp_disagree == 0)
stopifnot(n_base_dbp_disagree == 0)
stopifnot(n_base_missing == 0)

# --- Validation: 6mo vs. roster FOLLOWUP_SBP / FOLLOWUP_DBP / FOLLOWUP_MIDPOINT_DAY ---
check_fu <- roster %>%
  select(STUDY_ID, FOLLOWUP_SBP, FOLLOWUP_DBP, FOLLOWUP_MIDPOINT_DAY) %>%
  left_join(
    mo6_r %>% rename(FOLLOWUP_SBP_R = SBP, FOLLOWUP_DBP_R = DBP,
                     FOLLOWUP_READING_TS = DATE_TIME,
                     FOLLOWUP_MIDPOINT_DAY_R = DAY_FROM_START),
    by = "STUDY_ID"
  )

n_fu_sbp_disagree <- sum(!near(check_fu$FOLLOWUP_SBP, check_fu$FOLLOWUP_SBP_R, tol = mmhg_tol), na.rm = TRUE)
n_fu_dbp_disagree <- sum(!near(check_fu$FOLLOWUP_DBP, check_fu$FOLLOWUP_DBP_R, tol = mmhg_tol), na.rm = TRUE)
n_fu_day_disagree <- sum(check_fu$FOLLOWUP_MIDPOINT_DAY != check_fu$FOLLOWUP_MIDPOINT_DAY_R, na.rm = TRUE)
n_fu_missing      <- sum(is.na(check_fu$FOLLOWUP_SBP_R))
cat("\n6-month (followup) validation vs. roster:\n")
cat("  SBP disagreements:         ", n_fu_sbp_disagree, "\n")
cat("  DBP disagreements:         ", n_fu_dbp_disagree, "\n")
cat("  MIDPOINT_DAY disagreements:", n_fu_day_disagree, "\n")
cat("  Users with no R-derived 6mo reading:", n_fu_missing, "\n")
if (n_fu_sbp_disagree > 0 || n_fu_dbp_disagree > 0 || n_fu_day_disagree > 0) {
  cat("  First few disagreements:\n")
  print(head(check_fu %>%
               filter(!near(FOLLOWUP_SBP, FOLLOWUP_SBP_R, tol = mmhg_tol) |
                      !near(FOLLOWUP_DBP, FOLLOWUP_DBP_R, tol = mmhg_tol) |
                      FOLLOWUP_MIDPOINT_DAY != FOLLOWUP_MIDPOINT_DAY_R), 10))
}
stopifnot(n_fu_sbp_disagree == 0)
stopifnot(n_fu_dbp_disagree == 0)
stopifnot(n_fu_day_disagree == 0)
stopifnot(n_fu_missing == 0)

# --- Validation: roster DELTA_SBP / DELTA_DBP from R-derived (6mo − baseline) ---
delta_check <- roster %>%
  select(STUDY_ID, DELTA_SBP, DELTA_DBP) %>%
  left_join(
    check_base %>% select(STUDY_ID, BASELINE_SBP_R, BASELINE_DBP_R),
    by = "STUDY_ID"
  ) %>%
  left_join(
    check_fu %>% select(STUDY_ID, FOLLOWUP_SBP_R, FOLLOWUP_DBP_R),
    by = "STUDY_ID"
  ) %>%
  mutate(
    DELTA_SBP_R = FOLLOWUP_SBP_R - BASELINE_SBP_R,
    DELTA_DBP_R = FOLLOWUP_DBP_R - BASELINE_DBP_R
  )

n_delta_sbp_disagree <- sum(!near(delta_check$DELTA_SBP, delta_check$DELTA_SBP_R, tol = mmhg_tol), na.rm = TRUE)
n_delta_dbp_disagree <- sum(!near(delta_check$DELTA_DBP, delta_check$DELTA_DBP_R, tol = mmhg_tol), na.rm = TRUE)
cat("\nDelta validation vs. roster (DELTA = followup − baseline):\n")
cat("  SBP delta disagreements:", n_delta_sbp_disagree, "\n")
cat("  DBP delta disagreements:", n_delta_dbp_disagree, "\n")
stopifnot(n_delta_sbp_disagree == 0)
stopifnot(n_delta_dbp_disagree == 0)

# --- Diagnostic (informational, non-fatal): roster STUDY_START_DATE vs. ---
# min(reading date) per user in BP_MEASURES. Per Christian: should be equivalent
# in principle but minor disagreements possible due to source-data structure
# (e.g., readings logged before formal study enrollment). Not expected to be
# substantial; surfaced here so we can verify.
derived_start <- bp %>%
  group_by(STUDY_ID) %>%
  summarise(MIN_READING_DATE = min(READING_DATE), .groups = "drop")

start_diag <- roster %>%
  select(STUDY_ID, STUDY_START_DATE) %>%
  left_join(derived_start, by = "STUDY_ID") %>%
  mutate(differ = STUDY_START_DATE != MIN_READING_DATE)
n_start_differ <- sum(start_diag$differ, na.rm = TRUE)
n_start_missing <- sum(is.na(start_diag$MIN_READING_DATE))
cat("\nSTUDY_START_DATE diagnostic (info only — does not gate):\n")
cat("  Users where roster STUDY_START_DATE != min(BP reading date):", n_start_differ, "\n")
cat("  Users with no BP_MEASURES rows at all:", n_start_missing, "\n")
if (n_start_differ > 0) {
  cat("  Distribution of (min_reading_date − STUDY_START_DATE) in days:\n")
  print(start_diag %>%
          filter(differ) %>%
          mutate(diff_days = as.integer(MIN_READING_DATE - STUDY_START_DATE)) %>%
          count(diff_days) %>%
          arrange(diff_days))
}

# --- Assemble per-user 4-timepoint frame (one row per roster user) ---
bp_timepoints <- roster %>%
  select(STUDY_ID, STUDY_START_DATE) %>%
  left_join(
    baseline_r %>% rename(BASELINE_SBP = SBP, BASELINE_DBP = DBP,
                          BASELINE_READING_TS = DATE_TIME,
                          BASELINE_DAY = DAY_FROM_START),
    by = "STUDY_ID"
  ) %>%
  left_join(
    mo3_r %>% rename(MO3_SBP = SBP, MO3_DBP = DBP,
                     MO3_READING_TS = DATE_TIME, MO3_DAY = DAY_FROM_START),
    by = "STUDY_ID"
  ) %>%
  left_join(
    mo6_r %>% rename(MO6_SBP = SBP, MO6_DBP = DBP,
                     MO6_READING_TS = DATE_TIME, MO6_DAY = DAY_FROM_START),
    by = "STUDY_ID"
  ) %>%
  left_join(
    mo12_r %>% rename(MO12_SBP = SBP, MO12_DBP = DBP,
                      MO12_READING_TS = DATE_TIME, MO12_DAY = DAY_FROM_START),
    by = "STUDY_ID"
  )

stopifnot(nrow(bp_timepoints) == nrow(roster))
stopifnot(!any(duplicated(bp_timepoints$STUDY_ID)))

cat("\nbp_timepoints frame:", nrow(bp_timepoints), "rows,",
    ncol(bp_timepoints), "cols\n")
cat("  Have baseline reading:", sum(!is.na(bp_timepoints$BASELINE_SBP)), "\n")
cat("  Have 3mo reading:     ", sum(!is.na(bp_timepoints$MO3_SBP)),      "\n")
cat("  Have 6mo reading:     ", sum(!is.na(bp_timepoints$MO6_SBP)),      "\n")
cat("  Have 12mo reading:    ", sum(!is.na(bp_timepoints$MO12_SBP)),     "\n")

# --- Write ---
write_rds(bp_timepoints, file.path(derived_dir, "01_bp_timepoints.rds"))
write_csv(bp_timepoints, file.path(derived_dir, "01_bp_timepoints.csv"))
cat("\nWrote:\n")
cat("  ", file.path(derived_dir, "01_bp_timepoints.rds"), "\n")
cat("  ", file.path(derived_dir, "01_bp_timepoints.csv"), "\n")

