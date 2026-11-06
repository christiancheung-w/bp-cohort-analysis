# 09_glp1_flag.R
# Standalone additive slice (NOT wired into 99_assemble_user_frame.R).
# Builds a user-level GLP-1 self-report flag, analogous to the antihypertensive
# HTN_MED_BASELINE binary, for use as a separate sensitivity-analysis partition.
# See Build/BUILD_TODO.md; join to user_analytic_frame_1X by STUDY_ID at analysis time.
#
# Rationale: unlike HTN meds, there is NO dedicated GLP-1 pull. The only signal is
# the per-user-day GLP1 journal flag (behavior_tracker_id = 203) buried in the daily
# sleep/activity/recovery pulls. This slice aggregates those daily flags up to one
# row per roster user.
#
# Inputs (both daily pulls, so "ever" spans baseline + study windows):
#   Roster CSV                                        — STUDY_ID universe (1X cohort)
#   2026-05-13 Pulls/...DAILY_SLEEP..._RECOVERY_1X    — study-window daily GLP1 flag
#   2026-05-21 Pulls/...RECOVERY_BASELINE_1X          — baseline-window daily GLP1 flag
#
# Window: the study pull extends 180 days past STUDY_END_DATE (the 6-month reading).
# Rows after STUDY_END_DATE are dropped, so GLP1_EVER covers 30 days before baseline
# through the 6-month reading, i.e. it ends at the outcome, like HTN_MED_BASELINE
# (whose look-back is 180 days; the GLP-1 look-back is limited to 30 days by the pull).
#
# Derives (one row per roster user):
#   - GLP1_LOG_COUNT_STUDY     # of study-window days with GLP1 = TRUE
#   - GLP1_LOG_COUNT_BASELINE  # of baseline-window days with GLP1 = TRUE
#   - GLP1_LOG_COUNT           total TRUE-days across both windows
#   - GLP1_EVER                logical; TRUE if GLP1_LOG_COUNT > 0, else FALSE
#                              (users with no daily rows default to FALSE, mirroring
#                              HTN's never_logged default)
#
# NB: journal self-report — a GLP1_EVER == FALSE user is indistinguishable from a
# true non-user who simply never logged. Under-ascertainment is a known limitation
# (same caveat as HTN_MED_BASELINE). Window-split counts are retained to support a
# later timing-aware refinement without re-pulling.

library(tidyverse)

project_dir   <- "./PRIMARY"
data_dir      <- file.path(project_dir, "Query Results")
roster_dir    <- file.path(project_dir, "Roster")
build_dir     <- file.path(project_dir, "Build")
derived_dir   <- file.path(build_dir, "derived")
pull_dir_0513 <- file.path(data_dir, "2026-05-13 Pulls")
pull_dir_0521 <- file.path(data_dir, "2026-05-21 Pulls")

dir.create(derived_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read inputs (only the two columns we need from the 8.8M-row daily pulls) ---
roster <- read_csv(file.path(roster_dir, "ROSTER_1X_2026-05-12.csv"))

read_glp1_daily <- function(path, with_dates = FALSE) {
  cols_keep <- c("STUDY_ID", "GLP1", if (with_dates) c("LOWER_DAYS", "STUDY_END_DATE"))
  read_csv(path, col_select = all_of(cols_keep),
           col_types = cols(STUDY_ID = col_character(), GLP1 = col_character(),
                            .default = col_date()))
}
daily_study    <- read_glp1_daily(file.path(pull_dir_0513,
                    "DAILY_SLEEP_ACTIVITIES_RECOVERY_1X_2026-05-13.csv"),
                    with_dates = TRUE)
daily_baseline <- read_glp1_daily(file.path(pull_dir_0521,
                    "DAILY_SLEEP_ACTIVITIES_RECOVERY_BASELINE_1X_2026-05-21.csv"))

# Truncate the study window at the 6-month reading (drop the +180-day tail of the pull).
n_before <- nrow(daily_study)
daily_study <- daily_study %>% filter(LOWER_DAYS <= STUDY_END_DATE)
cat("Study daily rows after STUDY_END_DATE dropped:", n_before - nrow(daily_study),
    sprintf("(%.1f%% of %s)\n", 100 * (n_before - nrow(daily_study)) / n_before, format(n_before, big.mark = ",")))

# --- Pre-flight ---
cat("Roster rows:          ", nrow(roster),         "  users:", n_distinct(roster$STUDY_ID),         "\n")
cat("Study daily rows:     ", nrow(daily_study),    "  users:", n_distinct(daily_study$STUDY_ID),    "\n")
cat("Baseline daily rows:  ", nrow(daily_baseline), "  users:", n_distinct(daily_baseline$STUDY_ID), "\n")

stopifnot(!any(duplicated(roster$STUDY_ID)))

# Daily-pull users should live within the 1X roster universe.
stray_study    <- setdiff(unique(daily_study$STUDY_ID),    roster$STUDY_ID)
stray_baseline <- setdiff(unique(daily_baseline$STUDY_ID), roster$STUDY_ID)
cat("Daily users NOT in roster — study:", length(stray_study),
    " baseline:", length(stray_baseline), "\n")
stopifnot(length(stray_study) == 0, length(stray_baseline) == 0)

# --- Aggregate daily GLP1 flags to user level ---
# Robust to any casing/quoting: compare trimmed string to "TRUE". Blank -> not-TRUE.
count_true <- function(df, out_name) {
  df %>%
    mutate(is_true = str_trim(GLP1) == "TRUE" & !is.na(GLP1)) %>%
    group_by(STUDY_ID) %>%
    summarise("{out_name}" := sum(is_true), .groups = "drop")
}
glp1_study    <- count_true(daily_study,    "GLP1_LOG_COUNT_STUDY")
glp1_baseline <- count_true(daily_baseline, "GLP1_LOG_COUNT_BASELINE")

# --- Assemble onto roster universe (absent users -> 0 / FALSE) ---
glp1_flag <- roster %>%
  select(STUDY_ID) %>%
  left_join(glp1_study,    by = "STUDY_ID") %>%
  left_join(glp1_baseline, by = "STUDY_ID") %>%
  mutate(
    GLP1_LOG_COUNT_STUDY    = replace_na(GLP1_LOG_COUNT_STUDY,    0L),
    GLP1_LOG_COUNT_BASELINE = replace_na(GLP1_LOG_COUNT_BASELINE, 0L),
    GLP1_LOG_COUNT          = GLP1_LOG_COUNT_STUDY + GLP1_LOG_COUNT_BASELINE,
    GLP1_EVER               = GLP1_LOG_COUNT > 0
  )

# --- Invariants ---
stopifnot(nrow(glp1_flag) == nrow(roster))
stopifnot(!any(duplicated(glp1_flag$STUDY_ID)))
stopifnot(!any(is.na(glp1_flag$GLP1_EVER)))

# --- Summary prints ---
cat("\nglp1_flag frame:", nrow(glp1_flag), "rows,", ncol(glp1_flag), "cols\n")
cat("\nGLP1_EVER distribution:\n")
print(glp1_flag %>% count(GLP1_EVER, name = "n") %>%
      mutate(pct = round(100 * n / sum(n), 2)))
cat("\nDistinct GLP1 ever-users by window (pre-union):\n")
cat("  study-window only :", sum(glp1_flag$GLP1_LOG_COUNT_STUDY    > 0), "\n")
cat("  baseline-window   :", sum(glp1_flag$GLP1_LOG_COUNT_BASELINE > 0), "\n")
cat("  union (GLP1_EVER) :", sum(glp1_flag$GLP1_EVER), "\n")
cat("\nGLP1_LOG_COUNT summary among ever-users:\n")
print(summary(glp1_flag$GLP1_LOG_COUNT[glp1_flag$GLP1_EVER]))

# --- Write ---
write_rds(glp1_flag, file.path(derived_dir, "09_glp1_flag.rds"))
write_csv(glp1_flag, file.path(derived_dir, "09_glp1_flag.csv"))
cat("\nWrote:\n")
cat("  ", file.path(derived_dir, "09_glp1_flag.rds"), "\n")
cat("  ", file.path(derived_dir, "09_glp1_flag.csv"), "\n")
