# 07_engagement.R
# Phase 1 / Slice 7 of the user-level frame build (see Build/BUILD_TODO.md).
# Computes §12.7 engagement predictors for §8.
#
# Inputs:
#   Roster CSV                                                   — STUDY_ID universe
#   Query Results/2026-05-12 Pulls/APP_ENGAGEMENT_1X             — pre-aggregated 6mo behavior summary
#   Query Results/2026-05-21 Pulls/MET_OPTIMAL_STRAIN_1X         — WHOOP-canonical optimal-strain
#   Query Results/2026-05-12 Pulls/BP_MEASURES_1X                — count BP readings in [0, 180]
#   Query Results/2026-06-11 Pulls/DAILYSLEEP_TARGET_1X          — per-user avg sleep-need fulfillment (%)
#
# Derives (one row per roster user, n = 28,290):
#   - PCT_DAYS_JOURNALED              from APP_ENGAGEMENT
#   - AVG_DAILY_TIME_IN_APP_MIN       raw, from APP_ENGAGEMENT (retained for reference)
#   - LOG_AVG_DAILY_TIME_IN_APP_MIN   log(AVG + 1), SAP §12.7 form
#   - PCT_DAYS_ENGAGED_NOTIFICATIONS  from APP_ENGAGEMENT
#   - PCT_DAYS_USED_COACH             from APP_ENGAGEMENT
#   - SLEEP_PCT_OF_TARGET             from APP_ENGAGEMENT
#   - PCT_DAYS_MET_OPTIMAL_STRAIN     from MET_OPTIMAL_STRAIN pull (resolved Re-pull #2)
#   - N_BP_READINGS_STUDY             raw count, BP_MEASURES rows in [0, 180]
#   - LOG_N_BP_READINGS_STUDY         log(N + 1), SAP §12.7 form
#   - AVG_SLEEP_NEED_PCT              from DAILYSLEEP_TARGET pull — mean % of nightly
#                                     sleep-need fulfilled over the study window
#                                     (2026-06-11; replaces SLEEP_PCT_OF_TARGET as the
#                                     §8 sleep predictor — far better coverage)
#   - NIGHTS_WITH_SLEEP_DATA          coverage metadata for AVG_SLEEP_NEED_PCT
#
# Validation: pre-flight (STUDY_ID alignment, no infiltration, no duplicates) +
# row-count invariant after each join.

library(tidyverse)

project_dir   <- "./PRIMARY"
data_dir      <- file.path(project_dir, "Query Results")
roster_dir    <- file.path(project_dir, "Roster")
build_dir     <- file.path(project_dir, "Build")
derived_dir   <- file.path(build_dir, "derived")
pull_dir_0512 <- file.path(data_dir, "2026-05-12 Pulls")
pull_dir_0521 <- file.path(data_dir, "2026-05-21 Pulls")
pull_dir_0611 <- file.path(data_dir, "2026-06-11 Pulls")

dir.create(derived_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read inputs ---
roster <- read_csv(file.path(roster_dir, "ROSTER_1X_2026-05-12.csv"))

engagement <- read_csv(
  file.path(pull_dir_0512, "APP_ENGAGEMENT_1X_2026-05-12.csv"),
  col_select = c("STUDY_ID",
                 "PCT_DAYS_JOURNALED",
                 "AVG_DAILY_TIME_IN_APP_MIN",
                 "PCT_DAYS_ENGAGED_NOTIFICATIONS",
                 "PCT_DAYS_USED_COACH",
                 "SLEEP_PCT_OF_TARGET")
)

optimal_strain <- read_csv(
  file.path(pull_dir_0521, "MET_OPTIMAL_STRAIN_1X_2026-05-21.csv"),
  col_select = c("STUDY_ID", "PCT_DAYS_MET_OPTIMAL_STRAIN")
)

bp <- read_csv(file.path(pull_dir_0512, "BP_MEASURES_1X_2026-05-12.csv"),
               col_select = c("STUDY_ID", "DATE_TIME"))

sleep_need <- read_csv(
  file.path(pull_dir_0611, "DAILYSLEEP_TARGET_1X_2026-06-11.csv"),
  col_select = c("STUDY_ID", "NIGHTS_WITH_SLEEP_DATA", "AVG_SLEEP_NEED_PCT")
)

# --- Pre-flight ---
cat("Roster users:          ", n_distinct(roster$STUDY_ID),         "\n")
cat("Engagement users:      ", n_distinct(engagement$STUDY_ID),     "\n")
cat("Optimal-strain users:  ", n_distinct(optimal_strain$STUDY_ID), "\n")
cat("BP_MEASURES rows:      ", nrow(bp),                            "\n")
cat("BP_MEASURES users:     ", n_distinct(bp$STUDY_ID),              "\n")
cat("Sleep-need users:      ", n_distinct(sleep_need$STUDY_ID),      "\n")

stopifnot(!any(duplicated(engagement$STUDY_ID)))
stopifnot(!any(duplicated(optimal_strain$STUDY_ID)))
stopifnot(!any(duplicated(sleep_need$STUDY_ID)))
stopifnot(all(engagement$STUDY_ID     %in% roster$STUDY_ID))
stopifnot(all(optimal_strain$STUDY_ID %in% roster$STUDY_ID))
stopifnot(all(bp$STUDY_ID             %in% roster$STUDY_ID))
stopifnot(all(sleep_need$STUDY_ID     %in% roster$STUDY_ID))

# --- N_BP_READINGS in per-user study window [STUDY_START_DATE, STUDY_END_DATE] ---
# Per SAP §12.7: "Count of manual-cuff BP readings the user logged in the
# 6-month study window." Each user's study window is bounded by their own
# STUDY_START_DATE / STUDY_END_DATE (typically ~180 days but varies); using
# per-user bounds is SAP-precise. Previous implementation used a fixed
# [0, 180] day cutoff which slightly mis-classified edge-case readings.
bp_counts <- bp %>%
  inner_join(roster %>% select(STUDY_ID, STUDY_START_DATE, STUDY_END_DATE),
             by = "STUDY_ID") %>%
  mutate(READING_DATE = as.Date(DATE_TIME)) %>%
  filter(READING_DATE >= STUDY_START_DATE,
         READING_DATE <= STUDY_END_DATE) %>%
  count(STUDY_ID, name = "N_BP_READINGS_STUDY")

cat("\nBP_MEASURES rows in per-user study window [START, END]:",
    sum(bp_counts$N_BP_READINGS_STUDY),
    "  users with ≥1 reading:", nrow(bp_counts), "\n")

# --- Assemble per-user output ---
engagement_frame <- roster %>%
  select(STUDY_ID) %>%
  left_join(engagement, by = "STUDY_ID") %>%
  left_join(optimal_strain, by = "STUDY_ID") %>%
  left_join(bp_counts, by = "STUDY_ID") %>%
  left_join(sleep_need, by = "STUDY_ID") %>%
  mutate(
    N_BP_READINGS_STUDY            = coalesce(N_BP_READINGS_STUDY, 0L),
    LOG_AVG_DAILY_TIME_IN_APP_MIN  = log(AVG_DAILY_TIME_IN_APP_MIN + 1),
    LOG_N_BP_READINGS_STUDY        = log(N_BP_READINGS_STUDY + 1)
  )

stopifnot(nrow(engagement_frame) == nrow(roster))
stopifnot(!any(duplicated(engagement_frame$STUDY_ID)))

# --- Coverage / summary prints ---
cat("\nengagement_frame:", nrow(engagement_frame), "rows,",
    ncol(engagement_frame), "cols\n")

cat("\nDistributions of §12.7 predictors:\n")
for (col in c("PCT_DAYS_JOURNALED",
              "AVG_DAILY_TIME_IN_APP_MIN",
              "LOG_AVG_DAILY_TIME_IN_APP_MIN",
              "PCT_DAYS_ENGAGED_NOTIFICATIONS",
              "PCT_DAYS_USED_COACH",
              "SLEEP_PCT_OF_TARGET",
              "PCT_DAYS_MET_OPTIMAL_STRAIN",
              "N_BP_READINGS_STUDY",
              "LOG_N_BP_READINGS_STUDY",
              "NIGHTS_WITH_SLEEP_DATA",
              "AVG_SLEEP_NEED_PCT")) {
  cat("\n", col, ":\n", sep = "")
  print(summary(engagement_frame[[col]]))
}

# --- Write ---
write_rds(engagement_frame, file.path(derived_dir, "07_engagement.rds"))
write_csv(engagement_frame, file.path(derived_dir, "07_engagement.csv"))
cat("\nWrote:\n")
cat("  ", file.path(derived_dir, "07_engagement.rds"), "\n")
cat("  ", file.path(derived_dir, "07_engagement.csv"), "\n")
