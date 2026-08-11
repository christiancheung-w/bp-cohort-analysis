# 05_whoop_window_means.R
# Phase 1 / Slice 5 of the user-level frame build (see Build/BUILD_TODO.md).
#
# Inputs:
#   Roster CSV
#   Query Results/2026-05-21 Pulls/DAILY_SLEEP_ACTIVITIES_RECOVERY_BASELINE_1X — days [-30, -1]
#   Query Results/2026-05-13 Pulls/DAILY_SLEEP_ACTIVITIES_RECOVERY_1X          — study window (day 0+)
#
# Derives (one row per roster user, n = 28,290):
#   5 metrics × 4 windows = 20 mean columns + 4 per-window day-counts:
#     Metrics: RHR, HRV, STRAIN, SLEEP_DURATION (hours), SLEEP_CONSISTENCY
#     Windows (per SAP §7.1 spec table + §6.2 / §6.1):
#       BASELINE: days [-30, -1] relative to STUDY_START_DATE       (source: BASELINE pull)
#       MO3:      days [75, 105]  (30-day window centered on day 90)  (source: study pull)
#       MO6:      days [165, 195] (30-day window centered on day 180) (source: study pull)
#       STUDY:    days [0, 180]   (full 6-month follow-up window)     (source: study pull)
#
# Coverage thresholds applied at compute time (per Christian's decision —
# daily granularity is lost after aggregation, so the threshold must be baked
# in here rather than deferred to the QC slice):
#   30-day windows (BASELINE, MO3, MO6): ≥8 of 30 days (§11 25% rule)
#   Study window: ≥45 of ~180 days (§11 25% rule)
# Below threshold → that window's 5 means are NA; N_DAYS_<window> retains the
# actual day count (preserves "below threshold" vs. "no data" distinction).
#
# Defensive: dedupe (STUDY_ID, LOWER_DAYS) duplicates in both pulls
# (BUILD_TODO 2026-05-21 notes 2 dup pairs in baseline pull).
#
# Validation: no roster gold-standard. Pre-flight (no infiltration), baseline
# pull day-range sanity check (all rows ∈ [-30, -1]), coverage prints.

library(tidyverse)

project_dir   <- "./PRIMARY"
data_dir      <- file.path(project_dir, "Query Results")
roster_dir    <- file.path(project_dir, "Roster")
build_dir     <- file.path(project_dir, "Build")
derived_dir   <- file.path(build_dir, "derived")
pull_dir_0513 <- file.path(data_dir, "2026-05-13 Pulls")
pull_dir_0521 <- file.path(data_dir, "2026-05-21 Pulls")

dir.create(derived_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read inputs ---
roster <- read_csv(file.path(roster_dir, "ROSTER_1X_2026-05-12.csv"))

cat("Reading baseline-window daily pull...\n")
baseline_daily <- read_csv(file.path(pull_dir_0521,
                                     "DAILY_SLEEP_ACTIVITIES_RECOVERY_BASELINE_1X_2026-05-21.csv"))

cat("Reading study-window daily pull (8.8M rows — may take ~20s)...\n")
study_daily <- read_csv(file.path(pull_dir_0513,
                                  "DAILY_SLEEP_ACTIVITIES_RECOVERY_1X_2026-05-13.csv"))

# --- Pre-flight ---
cat("\nBaseline daily rows:", nrow(baseline_daily),
    "  users:", n_distinct(baseline_daily$STUDY_ID), "\n")
cat("Study daily rows:   ", nrow(study_daily),
    "  users:", n_distinct(study_daily$STUDY_ID), "\n")
cat("Roster users:", n_distinct(roster$STUDY_ID), "\n")

stopifnot(all(baseline_daily$STUDY_ID %in% roster$STUDY_ID))
stopifnot(all(study_daily$STUDY_ID    %in% roster$STUDY_ID))

# --- Dedupe (STUDY_ID, LOWER_DAYS) — average duplicates ---
metric_cols <- c("RESTING_HEART_RATE", "HRV", "SCALED_DAY_STRAIN",
                 "MINUTES_OF_SLEEP", "SLEEP_CONSISTENCY")

dedupe_daily <- function(daily) {
  daily %>%
    select(STUDY_ID, LOWER_DAYS, STUDY_START_DATE, all_of(metric_cols)) %>%
    group_by(STUDY_ID, LOWER_DAYS, STUDY_START_DATE) %>%
    summarise(across(all_of(metric_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

n_baseline_before <- nrow(baseline_daily)
n_study_before    <- nrow(study_daily)
baseline_daily <- dedupe_daily(baseline_daily)
study_daily    <- dedupe_daily(study_daily)
cat("\nDedupe — baseline:", n_baseline_before, "→", nrow(baseline_daily),
    " (collapsed", n_baseline_before - nrow(baseline_daily), "duplicate (user, day) pairs)\n")
cat("Dedupe — study:   ", n_study_before, "→", nrow(study_daily),
    " (collapsed", n_study_before - nrow(study_daily), "duplicate (user, day) pairs)\n")

# --- Day-from-start ---
baseline_daily <- baseline_daily %>%
  mutate(DAY_FROM_START = as.integer(LOWER_DAYS - STUDY_START_DATE))
study_daily <- study_daily %>%
  mutate(DAY_FROM_START = as.integer(LOWER_DAYS - STUDY_START_DATE))

# Sanity: baseline pull should be entirely in [-30, -1]
n_baseline_out <- sum(baseline_daily$DAY_FROM_START < -30 |
                      baseline_daily$DAY_FROM_START >  -1)
cat("\nBaseline pull day-range sanity (expect [-30, -1]):\n")
cat("  Range observed:", min(baseline_daily$DAY_FROM_START), "to",
    max(baseline_daily$DAY_FROM_START), "\n")
cat("  Rows outside [-30, -1]:", n_baseline_out, "\n")
stopifnot(n_baseline_out == 0)

cat("Study pull day-range observed:",
    min(study_daily$DAY_FROM_START), "to",
    max(study_daily$DAY_FROM_START), "\n")

# --- Window-mean helper ---
# Filters daily frame to [window_lo, window_hi], groups by STUDY_ID, computes
# n_days + 5 metric means (with na.rm). If n_days < min_days, all means set to
# NA; n_days retains its actual value. Returns per-user frame with N_DAYS_<label>
# and 5 *_<label>_MEAN columns. SLEEP_DURATION is in HOURS (MINUTES_OF_SLEEP / 60).
compute_window_means <- function(daily, window_lo, window_hi, min_days, label) {
  out <- daily %>%
    filter(DAY_FROM_START >= window_lo, DAY_FROM_START <= window_hi) %>%
    group_by(STUDY_ID) %>%
    summarise(
      n_days             = n(),
      RHR_mean           = mean(RESTING_HEART_RATE, na.rm = TRUE),
      HRV_mean           = mean(HRV,                na.rm = TRUE),
      STRAIN_mean        = mean(SCALED_DAY_STRAIN,  na.rm = TRUE),
      SLEEP_DURATION_mean    = mean(MINUTES_OF_SLEEP / 60, na.rm = TRUE),
      SLEEP_CONSISTENCY_mean = mean(SLEEP_CONSISTENCY, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      below_threshold = n_days < min_days,
      across(ends_with("_mean"),
             ~ if_else(below_threshold | is.nan(.x), NA_real_, .x))
    ) %>%
    select(-below_threshold)

  # Rename to canonical: <METRIC>_<LABEL>_MEAN  and  <LABEL>_N_DAYS
  out %>%
    rename_with(~ paste0(sub("_mean$", "", .x), "_", label, "_MEAN"),
                ends_with("_mean")) %>%
    rename(!!paste0(label, "_N_DAYS") := n_days)
}

# --- Apply to each window ---
cat("\nComputing window means...\n")
baseline_means <- compute_window_means(baseline_daily,
                                       window_lo = -30, window_hi = -1,
                                       min_days = 8,    label = "BASELINE")
mo3_means      <- compute_window_means(study_daily,
                                       window_lo = 75,  window_hi = 105,
                                       min_days = 8,    label = "MO3")
mo6_means      <- compute_window_means(study_daily,
                                       window_lo = 165, window_hi = 195,
                                       min_days = 8,    label = "MO6")
study_means    <- compute_window_means(study_daily,
                                       window_lo = 0,   window_hi = 180,
                                       min_days = 45,   label = "STUDY")

# --- Assemble per-user output ---
whoop_means <- roster %>%
  select(STUDY_ID) %>%
  left_join(baseline_means, by = "STUDY_ID") %>%
  left_join(mo3_means,      by = "STUDY_ID") %>%
  left_join(mo6_means,      by = "STUDY_ID") %>%
  left_join(study_means,    by = "STUDY_ID")

stopifnot(nrow(whoop_means) == nrow(roster))
stopifnot(!any(duplicated(whoop_means$STUDY_ID)))

# --- Coverage prints ---
cat("\n=== Per-window coverage (users meeting min-days threshold) ===\n")
report_window <- function(label, min_days, total_window_days) {
  n_days_col <- paste0(label, "_N_DAYS")
  mean_col   <- paste0("RHR_", label, "_MEAN")  # any *_MEAN works
  n_days_vals <- whoop_means[[n_days_col]]
  n_with_data <- sum(!is.na(n_days_vals))
  n_below     <- sum(!is.na(n_days_vals) & n_days_vals < min_days)
  n_qualifies <- sum(!is.na(whoop_means[[mean_col]]))
  n_no_data   <- sum(is.na(n_days_vals))
  cat(sprintf("  %-10s (min %2d of %3d days):  n_pass=%5d  below=%5d  no_data=%5d  total=%d\n",
              label, min_days, total_window_days, n_qualifies, n_below, n_no_data, nrow(whoop_means)))
}
report_window("BASELINE", 8,  30)
report_window("MO3",      8,  30)
report_window("MO6",      8,  30)
report_window("STUDY",    45, 180)

cat("\n=== Window-mean distributions (passing users only) ===\n")
for (label in c("BASELINE", "MO3", "MO6", "STUDY")) {
  cat("\n--- ", label, " ---\n", sep = "")
  for (m in c("RHR", "HRV", "STRAIN", "SLEEP_DURATION", "SLEEP_CONSISTENCY")) {
    col <- paste0(m, "_", label, "_MEAN")
    cat("  ", col, ":\n", sep = "")
    print(summary(whoop_means[[col]]))
  }
}

cat("\nwhoop_means frame:", nrow(whoop_means), "rows,", ncol(whoop_means), "cols\n")

# --- Write ---
write_rds(whoop_means, file.path(derived_dir, "05_whoop_window_means.rds"))
write_csv(whoop_means, file.path(derived_dir, "05_whoop_window_means.csv"))
cat("\nWrote:\n")
cat("  ", file.path(derived_dir, "05_whoop_window_means.rds"), "\n")
cat("  ", file.path(derived_dir, "05_whoop_window_means.csv"), "\n")
