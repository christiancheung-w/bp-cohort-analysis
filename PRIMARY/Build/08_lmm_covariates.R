# 08_lmm_covariates.R
# Phase 1 / Slice 8 of the user-level frame build (see Build/BUILD_TODO.md).
# LAST Phase 1 substrate before QC (slice 98) + assembly (slice 99).
#
# REWRITTEN 2026-05-22 — now uses the user-location pull (Re-pull #3, resolved
# 2026-05-22). The parked Northern-hemisphere-only scaffold is replaced with a
# proper timezone-aware + hemisphere-aware derivation.
#
# Inputs:
#   Build/derived/01_bp_timepoints.rds                              — *_READING_TS (UTC POSIXct)
#   Query Results/2026-05-22 Pulls/LOCATION_1X       — STUDY_ID, TIMEZONE, COUNTRY
#
# Derives (one row per roster user, n = 28,290):
#   Passthrough/derived metadata:
#     - TIMEZONE       IANA tz database name (e.g., America/New_York)
#     - COUNTRY        ISO-3166-1 alpha-2 (e.g., US, AU)
#     - HEMISPHERE     "Northern" / "Southern" (derived from COUNTRY via static list)
#
#   Per-timepoint covariates (the §5.1 / §5.5 adjusted-LMM covariates):
#     - {BASELINE, MO3, MO6, MO12}_TIME_OF_DAY    decimal hours in user-local clock (0–24)
#     - {BASELINE, MO3, MO6, MO12}_SEASON         "winter"/"spring"/"summer"/"fall",
#                                                  meteorological convention, hemisphere-corrected
#
# Time conversion: UTC *_READING_TS → user-local POSIXct via lubridate::with_tz().
# Done via group-split by TIMEZONE (~50 unique tzs across 28k users → very fast).
#
# Season convention (meteorological, by user's LOCAL calendar month):
#   Northern hemisphere:    Dec/Jan/Feb=winter, Mar/Apr/May=spring,
#                           Jun/Jul/Aug=summer, Sep/Oct/Nov=fall
#   Southern hemisphere:    Dec/Jan/Feb=summer, Mar/Apr/May=fall,
#                           Jun/Jul/Aug=winter, Sep/Oct/Nov=spring

library(tidyverse)
library(lubridate)

project_dir   <- "./PRIMARY"
data_dir      <- file.path(project_dir, "Query Results")
build_dir     <- file.path(project_dir, "Build")
derived_dir   <- file.path(build_dir, "derived")
pull_dir_0522 <- file.path(data_dir, "2026-05-22 Pulls")

dir.create(derived_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read inputs ---
bp_timepoints <- read_rds(file.path(derived_dir, "01_bp_timepoints.rds"))
location <- read_csv(file.path(pull_dir_0522,
                               "LOCATION_1X_2026-05-22.csv"))

# --- Pre-flight ---
cat("bp_timepoints users:", n_distinct(bp_timepoints$STUDY_ID), "\n")
cat("Location pull users:", n_distinct(location$STUDY_ID), "\n")
stopifnot(!any(duplicated(location$STUDY_ID)))
stopifnot(setequal(location$STUDY_ID, bp_timepoints$STUDY_ID))
cat("Cohort match: identical STUDY_ID sets in bp_timepoints and location pull.\n")

# --- Hemisphere lookup ---
# Static list of ISO-3166-1 alpha-2 codes for Southern-hemisphere-majority
# countries (capital / population-center south of equator). Equator-line countries
# assigned by capital location. Anything not in this list defaults to Northern.
sh_countries <- c(
  "AO", "AR", "AU", "BO", "BR", "BW", "CD", "CG", "CL", "EC", "FJ", "ID",
  "KE", "LS", "MG", "MU", "MW", "MZ", "NA", "NZ", "PE", "PG", "PY", "RE",
  "SB", "SZ", "TL", "TO", "TZ", "UY", "VU", "WS", "ZA", "ZM", "ZW"
)

location <- location %>%
  mutate(HEMISPHERE = if_else(COUNTRY %in% sh_countries, "Southern", "Northern"))

cat("\nHemisphere distribution (n users):\n")
print(location %>% count(HEMISPHERE, name = "n"))

# --- Join location onto BP timepoints ---
df <- bp_timepoints %>%
  select(STUDY_ID, BASELINE_READING_TS, MO3_READING_TS, MO6_READING_TS, MO12_READING_TS) %>%
  left_join(location, by = "STUDY_ID")

stopifnot(nrow(df) == nrow(bp_timepoints))

# --- UTC → local-clock conversion (split-apply-combine by TIMEZONE) ---
# with_tz() requires a scalar tz argument, so we group by TIMEZONE and apply
# per group. Within each group all rows share a single tz, so with_tz operates
# on the POSIXct vector in one call.
#
# We extract local decimal-hour and local calendar month into numeric/integer
# columns immediately — POSIXct vectors with mixed tz attributes don't combine
# cleanly across groups.

local_hour <- function(utc_ts, tz_scalar) {
  local <- with_tz(utc_ts, tzone = tz_scalar)
  hour(local) + minute(local) / 60 + second(local) / 3600
}

local_month <- function(utc_ts, tz_scalar) {
  local <- with_tz(utc_ts, tzone = tz_scalar)
  as.integer(month(local))
}

cat("\nConverting UTC reading_ts to user-local clock by timezone group...\n")
df_local <- df %>%
  group_by(TIMEZONE) %>%
  group_modify(~ {
    tz_scalar <- .y$TIMEZONE
    .x %>% mutate(
      BASELINE_TIME_OF_DAY = local_hour(BASELINE_READING_TS, tz_scalar),
      BASELINE_LOCAL_MONTH = local_month(BASELINE_READING_TS, tz_scalar),
      MO3_TIME_OF_DAY      = local_hour(MO3_READING_TS,      tz_scalar),
      MO3_LOCAL_MONTH      = local_month(MO3_READING_TS,     tz_scalar),
      MO6_TIME_OF_DAY      = local_hour(MO6_READING_TS,      tz_scalar),
      MO6_LOCAL_MONTH      = local_month(MO6_READING_TS,     tz_scalar),
      MO12_TIME_OF_DAY     = local_hour(MO12_READING_TS,     tz_scalar),
      MO12_LOCAL_MONTH     = local_month(MO12_READING_TS,    tz_scalar)
    )
  }) %>%
  ungroup()

cat("Conversion complete.\n")

# --- Season derivation (local month + hemisphere) ---
month_to_season <- function(m, hemi) {
  case_when(
    is.na(m) | is.na(hemi)                          ~ NA_character_,
    hemi == "Northern" & m %in% c(12, 1, 2)         ~ "winter",
    hemi == "Northern" & m %in% c(3, 4, 5)          ~ "spring",
    hemi == "Northern" & m %in% c(6, 7, 8)          ~ "summer",
    hemi == "Northern" & m %in% c(9, 10, 11)        ~ "fall",
    hemi == "Southern" & m %in% c(12, 1, 2)         ~ "summer",
    hemi == "Southern" & m %in% c(3, 4, 5)          ~ "fall",
    hemi == "Southern" & m %in% c(6, 7, 8)          ~ "winter",
    hemi == "Southern" & m %in% c(9, 10, 11)        ~ "spring"
  )
}

lmm_covariates <- df_local %>%
  mutate(
    BASELINE_SEASON = month_to_season(BASELINE_LOCAL_MONTH, HEMISPHERE),
    MO3_SEASON      = month_to_season(MO3_LOCAL_MONTH,      HEMISPHERE),
    MO6_SEASON      = month_to_season(MO6_LOCAL_MONTH,      HEMISPHERE),
    MO12_SEASON     = month_to_season(MO12_LOCAL_MONTH,     HEMISPHERE)
  ) %>%
  select(STUDY_ID, TIMEZONE, COUNTRY, HEMISPHERE,
         BASELINE_TIME_OF_DAY, BASELINE_SEASON,
         MO3_TIME_OF_DAY,      MO3_SEASON,
         MO6_TIME_OF_DAY,      MO6_SEASON,
         MO12_TIME_OF_DAY,     MO12_SEASON)

stopifnot(nrow(lmm_covariates) == nrow(bp_timepoints))
stopifnot(!any(duplicated(lmm_covariates$STUDY_ID)))

# --- Coverage / summary prints ---
cat("\nlmm_covariates frame:", nrow(lmm_covariates), "rows,",
    ncol(lmm_covariates), "cols\n")

cat("\nNon-NA coverage per timepoint:\n")
for (tp in c("BASELINE", "MO3", "MO6", "MO12")) {
  ttd <- paste0(tp, "_TIME_OF_DAY")
  ssn <- paste0(tp, "_SEASON")
  cat(sprintf("  %-8s  TIME_OF_DAY non-NA = %5d   SEASON non-NA = %5d\n",
              tp,
              sum(!is.na(lmm_covariates[[ttd]])),
              sum(!is.na(lmm_covariates[[ssn]]))))
}

cat("\nTIME_OF_DAY distributions (decimal hours, 0–24, user-local clock):\n")
for (tp in c("BASELINE", "MO3", "MO6", "MO12")) {
  col <- paste0(tp, "_TIME_OF_DAY")
  cat("\n  ", col, ":\n", sep = "")
  print(summary(lmm_covariates[[col]]))
}

cat("\nSEASON distributions (hemisphere-corrected):\n")
for (tp in c("BASELINE", "MO3", "MO6", "MO12")) {
  col <- paste0(tp, "_SEASON")
  cat("\n  ", col, ":\n", sep = "")
  print(lmm_covariates %>% count(.data[[col]], name = "n"))
}

cat("\nCross-check: hemisphere × baseline season distribution\n")
print(lmm_covariates %>% count(HEMISPHERE, BASELINE_SEASON, name = "n"))

# --- Write ---
write_rds(lmm_covariates, file.path(derived_dir, "08_lmm_covariates.rds"))
write_csv(lmm_covariates, file.path(derived_dir, "08_lmm_covariates.csv"))
cat("\nWrote:\n")
cat("  ", file.path(derived_dir, "08_lmm_covariates.rds"), "\n")
cat("  ", file.path(derived_dir, "08_lmm_covariates.csv"), "\n")
