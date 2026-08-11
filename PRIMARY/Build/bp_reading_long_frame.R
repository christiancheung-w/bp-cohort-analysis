# bp_reading_long_frame.R
# Phase 2 build: BP-reading-level long frame for SAP §7.2 (within/between-decomposed
# daily WHOOP signals predicting same-day BP). Sibling to the Phase 1 user-level
# analytic frame; lives in the same QC'd-output directory.
#
# Grain: one row per BP reading per user (INCLUDES baseline readings; analyses
#        filter `IS_BASELINE_READING == FALSE` for the §7.2 LMM. The baseline-
#        reading row's `*_MINUS_BASELINE` is structurally 0; its 7-day lookback
#        uses pre-study daily signals — useful for descriptive / sensitivity work,
#        not for the §7.2 within-user slope estimate.).
#
# Inputs:
#   Roster CSV
#   BP_MEASURES (all readings, 2026-05-12 pull)
#   DAILY_SLEEP_ACTIVITIES_RECOVERY_1X (2026-05-13 pull, col_select 7 cols)
#   user_analytic_frame_1X.rds (Phase 1 output — source of user-level covariates
#                                AND of BASELINE_SBP/DBP, ensuring consistency with
#                                Phase 1's already-QC'd values)
#
# Output: Build/derived/qc/bp_reading_analytic_frame_1X.{rds, csv}
#
# Key columns (per row, ~28 total):
#   Identifiers       : STUDY_ID, READING_TS (POSIXct UTC), READING_DATE,
#                       STUDY_END_DATE (user constant),
#                       DAY_FROM_BASELINE (integer), IS_BASELINE_READING (lgl),
#                       PERIOD ("DURING" or "AFTER" — from source BP_MEASURES.PERIOD)
#   BP outcomes       : SBP, DBP (raw, QC-bounded)
#                       BASELINE_SBP, BASELINE_DBP (user constants from Phase 1)
#                       SBP_MINUS_BASELINE, DBP_MINUS_BASELINE
#   §7.2 predictors   : 7-day lookback [READING_DATE − 7, READING_DATE − 1]
#                       {RHR, HRV, STRAIN, SLEEP_DURATION, SLEEP_CONSISTENCY}_LOOKBACK_MEAN
#                       {...}_LOOKBACK_N_DAYS  (informational; the mean is NA'd if N < 5)
#   User covariates   : AGE_AT_START, SEX, BASELINE_BMI, BASELINE_AHA, HTN_MED_BASELINE
#
# >>> PERIOD-filter warning <<<
# The §7.2 LMM is the ONLY SAP section that needs `PERIOD == "DURING"` filtering.
# Other SAP sections (§5.1, §5.5, §6.1, §6.2, §7.1, §8) use day-offset windowing
# (per SAP §3) that intentionally includes some PERIOD=AFTER readings (e.g., the
# §5.5 MO12 window [335, 390] is ENTIRELY post-study). Filtering by PERIOD for
# those analyses would silently drop SAP-prescribed data. The user-level frame
# (Phase 1) deliberately omits PERIOD for this reason.
#
# QC applied inline:
#   - BP bounds: SBP [70, 250]; DBP [40, 150]   (same as Phase 1 slice 01 + 98)
#   - Daily WHOOP bounds before lookback aggregation: RHR [30, 130]; HRV [5, 300];
#     STRAIN [0, 21]; SLEEP_DURATION (raw minutes) [0, 900] (≤ 15 hr); SLEEP_CONSISTENCY [0, 100]
#   - §11 ≥5-of-7-day per-predictor coverage threshold: mean NA'd if N_DAYS < 5
#
# Implementation: data.table non-equi join for the 7-day lookback over 8.8M daily
# rows × ~289K BP readings. Runs in ~30s on a modern Mac.

library(tidyverse)
library(data.table)

project_dir   <- "./PRIMARY"
data_dir      <- file.path(project_dir, "Query Results")
roster_dir    <- file.path(project_dir, "Roster")
derived_dir   <- file.path(project_dir, "Build", "derived")
qc_dir        <- file.path(derived_dir, "qc")
pull_dir_0512 <- file.path(data_dir, "2026-05-12 Pulls")
pull_dir_0513 <- file.path(data_dir, "2026-05-13 Pulls")

dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read inputs ---
roster     <- read_csv(file.path(roster_dir, "ROSTER_1X_2026-05-12.csv"))
user_frame <- read_rds(file.path(qc_dir,    "user_analytic_frame_1X.rds"))

bp <- read_csv(file.path(pull_dir_0512,
                         "BP_MEASURES_1X_2026-05-12.csv"))

cat("Reading study daily pull (8.8M rows — col_select to 7 cols)...\n")
daily <- read_csv(
  file.path(pull_dir_0513, "DAILY_SLEEP_ACTIVITIES_RECOVERY_1X_2026-05-13.csv"),
  col_select = c("STUDY_ID", "LOWER_DAYS", "STUDY_START_DATE",
                 "RESTING_HEART_RATE", "HRV", "SCALED_DAY_STRAIN",
                 "MINUTES_OF_SLEEP", "SLEEP_CONSISTENCY")
)

# --- Pre-flight ---
cat("\nRoster users:    ", n_distinct(roster$STUDY_ID),    "\n")
cat("BP_MEASURES rows:", nrow(bp),    " users:", n_distinct(bp$STUDY_ID),    "\n")
cat("Daily rows:      ", nrow(daily), " users:", n_distinct(daily$STUDY_ID), "\n")
stopifnot(all(bp$STUDY_ID    %in% roster$STUDY_ID))
stopifnot(all(daily$STUDY_ID %in% roster$STUDY_ID))

# --- Defensive dedupe of daily on (STUDY_ID, LOWER_DAYS) ---
metric_cols <- c("RESTING_HEART_RATE", "HRV", "SCALED_DAY_STRAIN",
                 "MINUTES_OF_SLEEP", "SLEEP_CONSISTENCY")
n_before <- nrow(daily)
daily <- daily %>%
  group_by(STUDY_ID, LOWER_DAYS, STUDY_START_DATE) %>%
  summarise(across(all_of(metric_cols), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(across(all_of(metric_cols), ~ if_else(is.nan(.x), NA_real_, .x)))
cat("Dedupe daily — collapsed", n_before - nrow(daily), "duplicate (user, day) pairs\n")

# --- Apply QC bounds to daily WHOOP signals (BEFORE lookback aggregation) ---
qc_inrange <- function(x, lo, hi) if_else(!is.na(x) & x >= lo & x <= hi, x, NA_real_)
daily <- daily %>%
  mutate(
    RESTING_HEART_RATE = qc_inrange(RESTING_HEART_RATE, 30, 130),
    HRV                = qc_inrange(HRV,                 5, 300),
    SCALED_DAY_STRAIN  = qc_inrange(SCALED_DAY_STRAIN,   0,  21),
    MINUTES_OF_SLEEP   = qc_inrange(MINUTES_OF_SLEEP,    0, 900),   # ≤ 15 hr
    SLEEP_CONSISTENCY  = qc_inrange(SLEEP_CONSISTENCY,   0, 100)
  )

# --- Build BP-reading long base ---
bp <- bp %>%
  inner_join(roster %>% select(STUDY_ID, STUDY_START_DATE, STUDY_END_DATE), by = "STUDY_ID") %>%
  mutate(
    READING_TS         = DATE_TIME,
    READING_DATE       = as.Date(DATE_TIME),
    DAY_FROM_BASELINE  = as.integer(READING_DATE - STUDY_START_DATE),
    IS_BASELINE_READING = DAY_FROM_BASELINE == 0,
    SBP                = qc_inrange(SBP, 70, 250),
    DBP                = qc_inrange(DBP, 40, 150)
  )

# Attach baseline_SBP/DBP from QC'd Phase 1 user frame
baseline_per_user <- user_frame %>%
  select(STUDY_ID,
         BASELINE_SBP = BASELINE_SBP,
         BASELINE_DBP = BASELINE_DBP)

bp <- bp %>%
  left_join(baseline_per_user, by = "STUDY_ID") %>%
  mutate(
    SBP_MINUS_BASELINE = SBP - BASELINE_SBP,
    DBP_MINUS_BASELINE = DBP - BASELINE_DBP
  )

cat("\nBP reading frame:", nrow(bp), "rows\n")
cat("  Baseline rows:    ", sum(bp$IS_BASELINE_READING), "\n")
cat("  Non-baseline rows:", sum(!bp$IS_BASELINE_READING), "\n")
cat("  DAY_FROM_BASELINE range: [", min(bp$DAY_FROM_BASELINE), ",",
    max(bp$DAY_FROM_BASELINE), "]\n")

# --- 7-day lookback aggregation via data.table non-equi join ---
# For each BP reading at READING_DATE, find daily rows with
# LOWER_DAYS in [READING_DATE − 7, READING_DATE − 1] for the same user.
# Aggregate per BP reading: per-predictor mean + n_days (non-NA count).

bp_dt    <- as.data.table(bp)
daily_dt <- as.data.table(daily)

bp_dt[, reading_id := .I]
bp_dt[, `:=`(win_start = READING_DATE - 7,
             win_end   = READING_DATE - 1)]
daily_dt[, LOWER_DAYS := as.Date(LOWER_DAYS)]
setkey(daily_dt, STUDY_ID, LOWER_DAYS)

cat("\nComputing 7-day lookback per reading (data.table non-equi join)...\n")
lookback_agg <- daily_dt[
  bp_dt,
  on = .(STUDY_ID, LOWER_DAYS >= win_start, LOWER_DAYS <= win_end),
  .(reading_id            = i.reading_id,
    RHR_LOOKBACK_N_DAYS               = sum(!is.na(RESTING_HEART_RATE)),
    HRV_LOOKBACK_N_DAYS               = sum(!is.na(HRV)),
    STRAIN_LOOKBACK_N_DAYS            = sum(!is.na(SCALED_DAY_STRAIN)),
    SLEEP_DURATION_LOOKBACK_N_DAYS    = sum(!is.na(MINUTES_OF_SLEEP)),
    SLEEP_CONSISTENCY_LOOKBACK_N_DAYS = sum(!is.na(SLEEP_CONSISTENCY)),
    RHR_LOOKBACK_MEAN                 = mean(RESTING_HEART_RATE, na.rm = TRUE),
    HRV_LOOKBACK_MEAN                 = mean(HRV,                 na.rm = TRUE),
    STRAIN_LOOKBACK_MEAN              = mean(SCALED_DAY_STRAIN,   na.rm = TRUE),
    SLEEP_DURATION_LOOKBACK_MEAN      = mean(MINUTES_OF_SLEEP / 60, na.rm = TRUE),
    SLEEP_CONSISTENCY_LOOKBACK_MEAN   = mean(SLEEP_CONSISTENCY,   na.rm = TRUE)
  ),
  by = .EACHI
]

# Result has one row per BP reading (by .EACHI). The on= clause produces a
# data.table with STUDY_ID + TWO LOWER_DAYS columns (one for each bound) — that
# duplicate name fails as_tibble(). Subset to wanted columns in data.table land
# (which tolerates dup names) before converting.
keep_cols <- c("reading_id",
               "RHR_LOOKBACK_N_DAYS",                "HRV_LOOKBACK_N_DAYS",
               "STRAIN_LOOKBACK_N_DAYS",             "SLEEP_DURATION_LOOKBACK_N_DAYS",
               "SLEEP_CONSISTENCY_LOOKBACK_N_DAYS",
               "RHR_LOOKBACK_MEAN",                  "HRV_LOOKBACK_MEAN",
               "STRAIN_LOOKBACK_MEAN",               "SLEEP_DURATION_LOOKBACK_MEAN",
               "SLEEP_CONSISTENCY_LOOKBACK_MEAN")
lookback_agg <- as_tibble(lookback_agg[, ..keep_cols])

cat("Lookback aggregation complete:", nrow(lookback_agg), "rows.\n")

# Convert NaN means (from readings with zero matching daily rows) → NA
# AND apply ≥5-of-7-day per-predictor coverage threshold:
# if N_DAYS < 5, set MEAN to NA.
nan_to_na <- function(x) if_else(is.nan(x), NA_real_, x)
lookback_agg <- lookback_agg %>%
  mutate(across(ends_with("_MEAN"), nan_to_na)) %>%
  mutate(
    RHR_LOOKBACK_MEAN               = if_else(RHR_LOOKBACK_N_DAYS               < 5, NA_real_, RHR_LOOKBACK_MEAN),
    HRV_LOOKBACK_MEAN               = if_else(HRV_LOOKBACK_N_DAYS               < 5, NA_real_, HRV_LOOKBACK_MEAN),
    STRAIN_LOOKBACK_MEAN            = if_else(STRAIN_LOOKBACK_N_DAYS            < 5, NA_real_, STRAIN_LOOKBACK_MEAN),
    SLEEP_DURATION_LOOKBACK_MEAN    = if_else(SLEEP_DURATION_LOOKBACK_N_DAYS    < 5, NA_real_, SLEEP_DURATION_LOOKBACK_MEAN),
    SLEEP_CONSISTENCY_LOOKBACK_MEAN = if_else(SLEEP_CONSISTENCY_LOOKBACK_N_DAYS < 5, NA_real_, SLEEP_CONSISTENCY_LOOKBACK_MEAN)
  )

# --- Join lookback aggregates onto BP base ---
bp_with_lookback <- as_tibble(bp_dt) %>%
  select(-win_start, -win_end) %>%
  left_join(lookback_agg, by = "reading_id")

stopifnot(nrow(bp_with_lookback) == nrow(bp))

# --- Join user-level covariates from Phase 1 user frame ---
user_covars <- user_frame %>%
  select(STUDY_ID, AGE_AT_START, SEX, BASELINE_BMI, BASELINE_AHA, HTN_MED_BASELINE)

bp_reading_frame <- bp_with_lookback %>%
  left_join(user_covars, by = "STUDY_ID") %>%
  select(STUDY_ID, READING_TS, READING_DATE, STUDY_END_DATE,
         DAY_FROM_BASELINE, IS_BASELINE_READING, PERIOD,
         SBP, DBP, BASELINE_SBP, BASELINE_DBP,
         SBP_MINUS_BASELINE, DBP_MINUS_BASELINE,
         RHR_LOOKBACK_MEAN,               RHR_LOOKBACK_N_DAYS,
         HRV_LOOKBACK_MEAN,               HRV_LOOKBACK_N_DAYS,
         STRAIN_LOOKBACK_MEAN,            STRAIN_LOOKBACK_N_DAYS,
         SLEEP_DURATION_LOOKBACK_MEAN,    SLEEP_DURATION_LOOKBACK_N_DAYS,
         SLEEP_CONSISTENCY_LOOKBACK_MEAN, SLEEP_CONSISTENCY_LOOKBACK_N_DAYS,
         AGE_AT_START, SEX, BASELINE_BMI, BASELINE_AHA, HTN_MED_BASELINE)

stopifnot(nrow(bp_reading_frame) == nrow(bp))
stopifnot(setequal(bp_reading_frame$STUDY_ID, roster$STUDY_ID))

# --- Summary prints ---
cat("\n=== bp_reading_analytic_frame_1X ===\n")
cat("Rows:", nrow(bp_reading_frame), "  Cols:", ncol(bp_reading_frame), "\n")
cat("Users represented:", n_distinct(bp_reading_frame$STUDY_ID), "\n")
cat("Baseline rows:    ", sum(bp_reading_frame$IS_BASELINE_READING), "\n")
cat("Non-baseline rows:", sum(!bp_reading_frame$IS_BASELINE_READING), "\n")

cat("\nPer-predictor: n readings with lookback mean computable (N_DAYS ≥ 5):\n")
for (p in c("RHR", "HRV", "STRAIN", "SLEEP_DURATION", "SLEEP_CONSISTENCY")) {
  mean_col <- paste0(p, "_LOOKBACK_MEAN")
  ndays_col <- paste0(p, "_LOOKBACK_N_DAYS")
  n_ok    <- sum(!is.na(bp_reading_frame[[mean_col]]))
  n_below <- sum(bp_reading_frame[[ndays_col]] > 0 & bp_reading_frame[[ndays_col]] < 5)
  n_zero  <- sum(bp_reading_frame[[ndays_col]] == 0, na.rm = TRUE)
  n_na    <- sum(is.na(bp_reading_frame[[ndays_col]]))
  cat(sprintf("  %-22s   n_pass = %6d   below_thresh = %5d   zero_days = %5d   na = %5d\n",
              p, n_ok, n_below, n_zero, n_na))
}

cat("\nLookback-mean distributions (passing readings only):\n")
for (p in c("RHR", "HRV", "STRAIN", "SLEEP_DURATION", "SLEEP_CONSISTENCY")) {
  col <- paste0(p, "_LOOKBACK_MEAN")
  cat("  ", col, ":\n", sep = "")
  print(summary(bp_reading_frame[[col]]))
}

cat("\nBP outcome distributions:\n")
cat("  SBP:\n");  print(summary(bp_reading_frame$SBP))
cat("  DBP:\n");  print(summary(bp_reading_frame$DBP))
cat("  SBP_MINUS_BASELINE:\n");  print(summary(bp_reading_frame$SBP_MINUS_BASELINE))
cat("  DBP_MINUS_BASELINE:\n");  print(summary(bp_reading_frame$DBP_MINUS_BASELINE))

# --- Write ---
write_rds(bp_reading_frame, file.path(qc_dir, "bp_reading_analytic_frame_1X.rds"))
write_csv(bp_reading_frame, file.path(qc_dir, "bp_reading_analytic_frame_1X.csv"))
cat("\nWrote:\n")
cat("  ", file.path(qc_dir, "bp_reading_analytic_frame_1X.rds"), "\n")
cat("  ", file.path(qc_dir, "bp_reading_analytic_frame_1X.csv"), "\n")

cat("\n=== Phase 2 build complete ===\n")
cat("BP-reading-level analytic frame ready for §7.2.\n")
cat("Analysts: filter `IS_BASELINE_READING == FALSE` + §2.2 cohort + ≥2 non-baseline\n")
cat("readings + apply within/between decomposition at fit time.\n")
