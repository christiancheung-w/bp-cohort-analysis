# 06_le8_tobacco_nicotine.R
# Phase 1 / Slice 6 of the user-level frame build (see Build/BUILD_TODO.md).
# Computes the §12.4 LE8 component 3 tobacco/nicotine predictors. The other
# §12.4 LE8 predictors (sleep duration mean, strain mean over the 6-month study
# window) are already derived in slice 5 as SLEEP_DURATION_STUDY_MEAN and
# STRAIN_STUDY_MEAN.
#
# Inputs:
#   Roster CSV                                                  — STUDY_ID universe
#   Query Results/2026-05-13 Pulls/DAILY_SLEEP_ACTIVITIES_RECOVERY_1X
#     — TOBACCO, NICOTINE logical (TRUE/FALSE/NULL) per user-day
#
# Derives (one row per roster user, n = 28,290):
#   - N_TOBACCO_NICOTINE_USE_DAYS    Count of study-window days [0, 180]
#                                    with (TOBACCO == TRUE OR NICOTINE == TRUE).
#                                    Both-logged days count once (per SAP §12.4).
#   - N_TOBACCO_NICOTINE_LOGGED_DAYS Count of study-window days with ≥1 non-null
#                                    entry for TOBACCO or NICOTINE (includes
#                                    "no" answers; data-coverage basis only).
#   - LOG_TOBACCO_NICOTINE_USE_DAYS  log(N_USE_DAYS + 1).  Per SAP §12.4 — LE8
#                                    primary predictor, intensity form.
#   - TOBACCO_NICOTINE_EVER_USED     Binary: TRUE if N_USE_DAYS > 0. Per SAP
#                                    §12.4 — LE8 primary predictor, binary
#                                    user-vs-non-user indicator (affirmative
#                                    tobacco/nicotine use on ≥1 study-window day).
#   - TOBACCO_NICOTINE_EVER_LOGGED   Binary: TRUE if N_LOGGED_DAYS > 0. Tracker-
#                                    engagement / data-coverage flag ONLY (counts
#                                    "no" answers); NOT a use indicator and NOT a
#                                    model predictor — retained for completeness
#                                    reporting.
#
# Users with no daily-pull rows OR only NA TOBACCO/NICOTINE entries land as:
#   N_USE_DAYS = 0; N_LOGGED_DAYS = 0; LOG_USE_DAYS = 0;
#   EVER_USED = FALSE; EVER_LOGGED = FALSE.
# Per SAP §12.4: "0 if the user never logs either tracker."

library(tidyverse)

project_dir   <- "./PRIMARY"
data_dir      <- file.path(project_dir, "Query Results")
roster_dir    <- file.path(project_dir, "Roster")
build_dir     <- file.path(project_dir, "Build")
derived_dir   <- file.path(build_dir, "derived")
pull_dir_0513 <- file.path(data_dir, "2026-05-13 Pulls")

dir.create(derived_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read inputs (just the columns we need from the 8.8M-row pull) ---
roster <- read_csv(file.path(roster_dir, "ROSTER_1X_2026-05-12.csv"))

cat("Reading study daily pull (TOBACCO/NICOTINE columns only)...\n")
daily <- read_csv(
  file.path(pull_dir_0513, "DAILY_SLEEP_ACTIVITIES_RECOVERY_1X_2026-05-13.csv"),
  col_select = c("STUDY_ID", "LOWER_DAYS", "STUDY_START_DATE", "TOBACCO", "NICOTINE")
)

# --- Pre-flight ---
cat("\nRoster users:", n_distinct(roster$STUDY_ID), "\n")
cat("Daily pull rows:", nrow(daily), "  users:", n_distinct(daily$STUDY_ID), "\n")
stopifnot(all(daily$STUDY_ID %in% roster$STUDY_ID))

# --- Filter to study window [0, 180] per SAP §12.4 ---
daily <- daily %>%
  mutate(DAY_FROM_START = as.integer(LOWER_DAYS - STUDY_START_DATE))

window_daily <- daily %>%
  filter(DAY_FROM_START >= 0, DAY_FROM_START <= 180)

cat("Daily rows in study window [0, 180]:", nrow(window_daily), "\n")

# --- Defensive dedupe + per-day flags ---
# Use any() to collapse duplicate (STUDY_ID, LOWER_DAYS) rows safely:
#   - has_use: at least one source reported TOBACCO=TRUE or NICOTINE=TRUE that day
#   - has_log: at least one source had a non-NA entry that day (TRUE or FALSE)
per_day <- window_daily %>%
  group_by(STUDY_ID, LOWER_DAYS) %>%
  summarise(
    has_use = any(TOBACCO | NICOTINE,           na.rm = TRUE),
    has_log = any(!is.na(TOBACCO) | !is.na(NICOTINE)),
    .groups = "drop"
  )

cat("Unique (user, day) pairs in window:", nrow(per_day), "\n")

# --- Per-user counts ---
per_user <- per_day %>%
  group_by(STUDY_ID) %>%
  summarise(
    N_TOBACCO_NICOTINE_USE_DAYS    = sum(has_use),
    N_TOBACCO_NICOTINE_LOGGED_DAYS = sum(has_log),
    .groups = "drop"
  )

cat("Users with any window data:", nrow(per_user), "\n")

# --- Assemble per-roster-user output ---
# Roster users without any daily-pull rows in the study window → 0/FALSE per SAP
le8_tobacco_nicotine <- roster %>%
  select(STUDY_ID) %>%
  left_join(per_user, by = "STUDY_ID") %>%
  mutate(
    N_TOBACCO_NICOTINE_USE_DAYS    = coalesce(N_TOBACCO_NICOTINE_USE_DAYS,    0L),
    N_TOBACCO_NICOTINE_LOGGED_DAYS = coalesce(N_TOBACCO_NICOTINE_LOGGED_DAYS, 0L),
    LOG_TOBACCO_NICOTINE_USE_DAYS  = log(N_TOBACCO_NICOTINE_USE_DAYS + 1),
    TOBACCO_NICOTINE_EVER_USED     = N_TOBACCO_NICOTINE_USE_DAYS > 0,
    TOBACCO_NICOTINE_EVER_LOGGED   = N_TOBACCO_NICOTINE_LOGGED_DAYS > 0
  )

stopifnot(nrow(le8_tobacco_nicotine) == nrow(roster))
stopifnot(!any(duplicated(le8_tobacco_nicotine$STUDY_ID)))

# --- Coverage / summary prints ---
cat("\n=== Tobacco/Nicotine — use-day distribution (study window [0, 180]) ===\n")
print(summary(le8_tobacco_nicotine$N_TOBACCO_NICOTINE_USE_DAYS))

cat("\n=== Tobacco/Nicotine — logged-day distribution ===\n")
print(summary(le8_tobacco_nicotine$N_TOBACCO_NICOTINE_LOGGED_DAYS))

cat("\n=== EVER_USED binary distribution (§6.1 model predictor; N_USE_DAYS > 0) ===\n")
print(le8_tobacco_nicotine %>% count(TOBACCO_NICOTINE_EVER_USED, name = "n"))

cat("\n=== EVER_LOGGED binary distribution (coverage flag; NOT a predictor) ===\n")
print(le8_tobacco_nicotine %>% count(TOBACCO_NICOTINE_EVER_LOGGED, name = "n"))

cat("\n=== Users with N_USE_DAYS > 0 (any reported tobacco/nicotine use day) ===\n")
cat("  n =", sum(le8_tobacco_nicotine$N_TOBACCO_NICOTINE_USE_DAYS > 0),
    "of", nrow(le8_tobacco_nicotine),
    sprintf("(%.2f%%)", 100 * sum(le8_tobacco_nicotine$N_TOBACCO_NICOTINE_USE_DAYS > 0) / nrow(le8_tobacco_nicotine)),
    "\n")

cat("\n=== LOG_USE_DAYS distribution ===\n")
print(summary(le8_tobacco_nicotine$LOG_TOBACCO_NICOTINE_USE_DAYS))

cat("\nle8_tobacco_nicotine frame:", nrow(le8_tobacco_nicotine), "rows,",
    ncol(le8_tobacco_nicotine), "cols\n")

# --- Write ---
write_rds(le8_tobacco_nicotine, file.path(derived_dir, "06_le8_tobacco_nicotine.rds"))
write_csv(le8_tobacco_nicotine, file.path(derived_dir, "06_le8_tobacco_nicotine.csv"))
cat("\nWrote:\n")
cat("  ", file.path(derived_dir, "06_le8_tobacco_nicotine.rds"), "\n")
cat("  ", file.path(derived_dir, "06_le8_tobacco_nicotine.csv"), "\n")
