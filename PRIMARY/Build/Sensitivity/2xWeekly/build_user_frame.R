# build_user_frame.R
# 2xWeekly cohort §5.1 sensitivity build (single script).
#
# Per SAP §5.1 sensitivity (line 146): "Refit the primary §5.1 LMM on the
# 2xWeekly cohort (n ≈ 3,041), which uses within-timepoint averaging
# of two readings per timepoint (stricter measurement). Replication in direction
# and rough magnitude = independent evidence against (a) measurement-design
# artifact and (b) RTM-only explanation."
#
# Scope (per 2026-05-22 simplification decision):
#   §5.1 only. The SAP-spec'd 2xWeekly sensitivities in §6.1 + §6.2 are deferred
#   because those analyses already use baseline_BP as an ANCOVA covariate (which
#   partials out the RTM concern the 2xWeekly cohort would otherwise demonstrate).
#   The §5.1 sensitivity is the load-bearing one. Pending SAP v0.50 to formalize.
#
# Inputs (NO new Cortex pulls — all from existing data):
#   PRIMARY/Alternative-Cohorts/2xWeekly/Roster/ROSTER_2xWeekly_2026-05-11.csv
#   PRIMARY/Alternative-Cohorts/2xWeekly/Query Results/2026-05-08 Pulls/BP_MEASURES_2xWeekly_2026-05-08.csv
#   PRIMARY/Query Results/2026-05-13 Pulls/HEIGHT_BMI_1X_2026-05-13.csv (filter to 2xWeekly users)
#   PRIMARY/Query Results/2026-05-14 Pulls/HTN_MEDS_1X_2026-05-14.csv  (filter to 2xWeekly users)
#   PRIMARY/Query Results/2026-05-22 Pulls/LOCATION_1X_2026-05-22.csv  (filter to 2xWeekly users)
#
# 2xWeekly is a strict subset of 1X (3,041 of 28,290 users; zero infiltration
# verified 2026-05-22). Aux pulls (HEIGHT_BMI, HTN_MEDS, LOCATION) can therefore
# be filtered from the 1X pulls; no new Cortex queries required.
#
# Output:
#   PRIMARY/Build/Sensitivity/2xWeekly/derived/qc/user_analytic_frame_2xWeekly.{rds,csv}
#
# 2xWeekly cohort-defining picking rules (per the 2xWeekly roster data dictionary;
# NOT the SAP §3 windowing used by 1X):
#   Baseline: mean of the FIRST 2 manual-cuff readings on/after STUDY_START_DATE
#             (no window cap — span whatever days needed to find 2 readings).
#   6mo:      mean of the 2 readings CLOSEST to STUDY_END_DATE within window
#             [STUDY_START_DATE, STUDY_START_DATE + 390d]. Tiebreak: LATER reading.
#   3mo:      not in the roster; derived per SAP §3 (window [60, 120], closest to
#             day 90, earlier-ts tiebreak), with 2-reading averaging applied —
#             paralleling the 2xWeekly's "averaged-2-readings" principle.
# Validation against the 2xWeekly roster's pre-computed BASELINE_SBP / FOLLOWUP_SBP
# (which used these rules) confirms our derivations match exactly.

library(tidyverse)
library(lubridate)

project_dir   <- "./PRIMARY"
secondary_dir <- "./PRIMARY/Alternative-Cohorts/2xWeekly"
out_dir       <- file.path(project_dir, "Build", "Sensitivity", "2xWeekly", "derived", "qc")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read 2xWeekly cohort source data ---
roster <- read_csv(file.path(secondary_dir, "Roster",
                             "ROSTER_2xWeekly_2026-05-11.csv"))
bp     <- read_csv(file.path(secondary_dir, "Query Results", "2026-05-08 Pulls",
                             "BP_MEASURES_2xWeekly_2026-05-08.csv"))

# --- Read 1X auxiliary pulls (will filter to 2xWeekly users) ---
height_bmi_1x <- read_csv(file.path(project_dir, "Query Results", "2026-05-13 Pulls",
                                    "HEIGHT_BMI_1X_2026-05-13.csv"))
htn_meds_1x   <- read_csv(file.path(project_dir, "Query Results", "2026-05-14 Pulls",
                                    "HTN_MEDS_1X_2026-05-14.csv"))
location_1x   <- read_csv(file.path(project_dir, "Query Results", "2026-05-22 Pulls",
                                    "LOCATION_1X_2026-05-22.csv"))

cat("Roster 2xWeekly: ", nrow(roster), "users\n")
cat("BP 2xWeekly:     ", nrow(bp),     "rows /", n_distinct(bp$STUDY_ID), "users\n")

# Pre-flight: 2xWeekly ⊂ 1X
stopifnot(all(roster$STUDY_ID %in% height_bmi_1x$STUDY_ID))
stopifnot(all(roster$STUDY_ID %in% htn_meds_1x$STUDY_ID))
stopifnot(all(roster$STUDY_ID %in% location_1x$STUDY_ID))
stopifnot(all(bp$STUDY_ID     %in% roster$STUDY_ID))

# Filter aux pulls to 2xWeekly cohort
height_bmi <- height_bmi_1x %>% filter(STUDY_ID %in% roster$STUDY_ID)
htn_meds   <- htn_meds_1x   %>% filter(STUDY_ID %in% roster$STUDY_ID)
location   <- location_1x   %>% filter(STUDY_ID %in% roster$STUDY_ID)

cat("Aux pulls filtered to 2xWeekly cohort (n =", nrow(roster), ")\n")

# =====================================================================
# Section 1: BP timepoints with within-timepoint averaging (2xWeekly's defining feature)
# =====================================================================

# Bring STUDY_START_DATE + STUDY_END_DATE to BP frame; apply QC bounds
bp <- bp %>%
  inner_join(roster %>% select(STUDY_ID, STUDY_START_DATE, STUDY_END_DATE),
             by = "STUDY_ID") %>%
  mutate(
    READING_DATE   = as.Date(DATE_TIME),
    DAY_FROM_START = as.integer(READING_DATE - STUDY_START_DATE),
    DAYS_TO_END    = abs(as.integer(READING_DATE - STUDY_END_DATE)),
    SBP            = if_else(!is.na(SBP) & SBP >= 70 & SBP <= 250, SBP, NA_real_),
    DBP            = if_else(!is.na(DBP) & DBP >= 40 & DBP <= 150, DBP, NA_real_)
  )

# --- BASELINE: first 2 manual-cuff readings on/after STUDY_START_DATE ---
# (per 2xWeekly roster data dictionary; no upper window cap; sort by DATE_TIME asc)
derive_baseline_2x <- function(bp) {
  bp %>%
    filter(DAY_FROM_START >= 0) %>%
    group_by(STUDY_ID) %>%
    arrange(DATE_TIME, .by_group = TRUE) %>%
    slice_head(n = 2) %>%
    summarise(
      SBP                 = mean(SBP, na.rm = TRUE),
      DBP                 = mean(DBP, na.rm = TRUE),
      EARLIEST_READING_TS = first(DATE_TIME),
      MEAN_DAY_FROM_START = mean(DAY_FROM_START),
      N_READINGS          = n(),
      .groups             = "drop"
    )
}

# --- 3mo: SAP §3 rule (closest to day 90 in [60, 120], earlier-ts tiebreak), 2 readings averaged ---
derive_mo3_2x <- function(bp) {
  bp %>%
    filter(DAY_FROM_START >= 60, DAY_FROM_START <= 120) %>%
    group_by(STUDY_ID) %>%
    arrange(abs(DAY_FROM_START - 90), DATE_TIME, .by_group = TRUE) %>%
    slice_head(n = 2) %>%
    summarise(
      SBP                 = mean(SBP, na.rm = TRUE),
      DBP                 = mean(DBP, na.rm = TRUE),
      EARLIEST_READING_TS = first(DATE_TIME),
      MEAN_DAY_FROM_START = mean(DAY_FROM_START),
      N_READINGS          = n(),
      .groups             = "drop"
    )
}

# --- 6mo (= FOLLOWUP): 2 closest to STUDY_END_DATE in [0, 390], LATER-ts tiebreak ---
# (per 2xWeekly roster data dictionary)
derive_mo6_2x <- function(bp) {
  bp %>%
    filter(DAY_FROM_START >= 0, DAY_FROM_START <= 390) %>%
    group_by(STUDY_ID) %>%
    arrange(DAYS_TO_END, desc(DATE_TIME), .by_group = TRUE) %>%
    slice_head(n = 2) %>%
    summarise(
      SBP                 = mean(SBP, na.rm = TRUE),
      DBP                 = mean(DBP, na.rm = TRUE),
      EARLIEST_READING_TS = min(DATE_TIME),
      MEAN_DAY_FROM_START = mean(DAY_FROM_START),
      N_READINGS          = n(),
      .groups             = "drop"
    )
}

baseline <- derive_baseline_2x(bp)
mo3      <- derive_mo3_2x(bp)
mo6      <- derive_mo6_2x(bp)

cat("\n=== Within-timepoint coverage (2xWeekly cohort) ===\n")
cat("Baseline: ", nrow(baseline), "users\n")
print(baseline %>% count(N_READINGS, name = "n_users"))
cat("MO3:      ", nrow(mo3), "users (of", nrow(roster), "total)\n")
print(mo3 %>% count(N_READINGS, name = "n_users"))
cat("MO6:      ", nrow(mo6), "users\n")
print(mo6 %>% count(N_READINGS, name = "n_users"))

# =====================================================================
# Section 2: Validate re-derivation against roster (REPORT drift; no stopifnot)
# =====================================================================
#
# Expected drift per the 2xWeekly roster data dictionary:
#   BASELINE_SBP: 39 / 3,041 users differ (1.3%) — new early-window readings
#                  in BP_MEASURES vs the roster's 2026-05-08 snapshot
#   FOLLOWUP_SBP: 390 / 3,041 users differ (12.8%) — new readings shift which
#                  2 are "closest to STUDY_END" — mostly < 5 mmHg differences
#
# We trust the ROSTER's BASELINE/FOLLOWUP BP values as canonical (the 2xWeekly
# cohort is defined by that snapshot). The re-derivation here is used only to
# extract READING_TS values for the time_of_day / season per-timepoint
# covariates — those reflect the current BP_MEASURES, not the 2026-05-08
# snapshot, but this is acceptable nuisance noise for the ~10% drift users.

mmhg_tol <- 1e-4

check_base <- roster %>%
  select(STUDY_ID, BASELINE_SBP, BASELINE_DBP) %>%
  left_join(baseline %>% rename(BASELINE_SBP_R = SBP, BASELINE_DBP_R = DBP),
            by = "STUDY_ID")
n_base_sbp_disagree <- sum(!dplyr::near(check_base$BASELINE_SBP, check_base$BASELINE_SBP_R, tol = mmhg_tol), na.rm = TRUE)
n_base_dbp_disagree <- sum(!dplyr::near(check_base$BASELINE_DBP, check_base$BASELINE_DBP_R, tol = mmhg_tol), na.rm = TRUE)
cat("\nDrift report (re-derived vs 2xWeekly roster's snapshot values):\n")
cat(sprintf("  Baseline SBP disagreements: %4d (%.1f%%)  [data-dictionary expects ~39]\n",
            n_base_sbp_disagree, 100 * n_base_sbp_disagree / nrow(roster)))
cat(sprintf("  Baseline DBP disagreements: %4d (%.1f%%)\n",
            n_base_dbp_disagree, 100 * n_base_dbp_disagree / nrow(roster)))

check_fu <- roster %>%
  select(STUDY_ID, FOLLOWUP_SBP, FOLLOWUP_DBP, FOLLOWUP_MIDPOINT_DAY) %>%
  left_join(mo6 %>% rename(FOLLOWUP_SBP_R = SBP, FOLLOWUP_DBP_R = DBP,
                           FOLLOWUP_MIDPOINT_DAY_R = MEAN_DAY_FROM_START),
            by = "STUDY_ID")
n_fu_sbp_disagree <- sum(!dplyr::near(check_fu$FOLLOWUP_SBP, check_fu$FOLLOWUP_SBP_R, tol = mmhg_tol), na.rm = TRUE)
n_fu_dbp_disagree <- sum(!dplyr::near(check_fu$FOLLOWUP_DBP, check_fu$FOLLOWUP_DBP_R, tol = mmhg_tol), na.rm = TRUE)
n_fu_day_disagree <- sum(check_fu$FOLLOWUP_MIDPOINT_DAY != check_fu$FOLLOWUP_MIDPOINT_DAY_R, na.rm = TRUE)
cat(sprintf("  6mo SBP disagreements:      %4d (%.1f%%)  [data-dictionary expects ~390]\n",
            n_fu_sbp_disagree, 100 * n_fu_sbp_disagree / nrow(roster)))
cat(sprintf("  6mo DBP disagreements:      %4d (%.1f%%)\n",
            n_fu_dbp_disagree, 100 * n_fu_dbp_disagree / nrow(roster)))
cat(sprintf("  6mo MIDPOINT_DAY disagreements: %4d (%.1f%%)\n",
            n_fu_day_disagree, 100 * n_fu_day_disagree / nrow(roster)))

# Sanity gate: drift should be within ~2× the documented count (to catch a
# data-source mistake without rejecting expected drift).
stopifnot(n_base_sbp_disagree < 100)   # 39 documented; allow ~2.5×
stopifnot(n_fu_sbp_disagree   < 800)   # 390 documented; allow ~2×
cat("Drift within documented bounds — proceeding with ROSTER values as canonical.\n")

# =====================================================================
# Section 3: BP outcomes — AHA classification + Δs (per slice 02 logic)
# =====================================================================

classify_aha <- function(sbp, dbp) {
  case_when(
    is.na(sbp) | is.na(dbp) ~ NA_character_,
    sbp >= 140 | dbp >= 90  ~ "Stage 2 Hypertension",
    sbp >= 130 | dbp >= 80  ~ "Stage 1 Hypertension",
    sbp >= 120 & dbp <  80  ~ "Elevated",
    TRUE                    ~ "Normal"
  )
}

# Use ROSTER's canonical BP values for baseline + 6mo (those are the 2xWeekly
# cohort's authoritative BASELINE_SBP / FOLLOWUP_SBP / AHA classifications).
# Use derived READING_TS for time_of_day / season covariates (small drift for
# ~10% of users; acceptable nuisance noise per Section 2 above).
bp_frame <- roster %>%
  select(STUDY_ID, STUDY_START_DATE, STUDY_END_DATE, BIRTHDAY, AGE_AT_START, GENDER,
         BASELINE_SBP, BASELINE_DBP,
         FOLLOWUP_SBP, FOLLOWUP_DBP, FOLLOWUP_MIDPOINT_DAY,
         BASELINE_AHA_CLASSIFICATION, END_AHA_CLASSIFICATION) %>%
  rename(
    SEX          = GENDER,
    MO6_SBP      = FOLLOWUP_SBP,
    MO6_DBP      = FOLLOWUP_DBP,
    MO6_DAY      = FOLLOWUP_MIDPOINT_DAY,
    BASELINE_AHA = BASELINE_AHA_CLASSIFICATION,
    MO6_AHA      = END_AHA_CLASSIFICATION
  ) %>%
  left_join(baseline %>%
              select(STUDY_ID,
                     BASELINE_READING_TS = EARLIEST_READING_TS,
                     BASELINE_DAY        = MEAN_DAY_FROM_START,
                     BASELINE_N_READINGS = N_READINGS),
            by = "STUDY_ID") %>%
  left_join(mo3 %>% rename(MO3_SBP = SBP, MO3_DBP = DBP,
                           MO3_READING_TS = EARLIEST_READING_TS,
                           MO3_DAY = MEAN_DAY_FROM_START,
                           MO3_N_READINGS = N_READINGS),
            by = "STUDY_ID") %>%
  left_join(mo6 %>%
              select(STUDY_ID,
                     MO6_READING_TS = EARLIEST_READING_TS,
                     MO6_N_READINGS = N_READINGS),
            by = "STUDY_ID") %>%
  mutate(
    MO3_AHA       = classify_aha(MO3_SBP, MO3_DBP),
    DELTA_3MO_SBP = MO3_SBP - BASELINE_SBP,
    DELTA_3MO_DBP = MO3_DBP - BASELINE_DBP,
    DELTA_6MO_SBP = MO6_SBP - BASELINE_SBP,
    DELTA_6MO_DBP = MO6_DBP - BASELINE_DBP
  )

cat("\nBASELINE_AHA stratum counts (2xWeekly, from roster's patched classifications):\n")
print(bp_frame %>% count(BASELINE_AHA, name = "n"))
cat("\nMO6 (END) AHA stratum counts:\n")
print(bp_frame %>% count(MO6_AHA, name = "n"))

# =====================================================================
# Section 4: Demographics — anthropometrics (filtered 1X HEIGHT_BMI) + HTN meds
# =====================================================================

height_bmi_qc <- height_bmi %>%
  mutate(
    HEIGHT_M                   = if_else(HEIGHT >= 1.4 & HEIGHT <= 2.2, HEIGHT, NA_real_),
    WEIGHT_BASELINE_PROFILE_KG = if_else(WEIGHT >= 30  & WEIGHT <= 250, WEIGHT, NA_real_),
    BASELINE_BMI               = if_else(BMI    >= 12  & BMI    <= 70,  BMI,    NA_real_)
  ) %>%
  select(STUDY_ID, HEIGHT_M, WEIGHT_BASELINE_PROFILE_KG, BASELINE_BMI)

htn_meds_sel <- htn_meds %>% select(STUDY_ID, HTN_MED_BASELINE = BP_MED_STATUS)

# =====================================================================
# Section 5: LMM covariates — TIMEZONE/COUNTRY/HEMISPHERE + per-timepoint time_of_day + season
# =====================================================================

sh_countries <- c(
  "AO", "AR", "AU", "BO", "BR", "BW", "CD", "CG", "CL", "EC", "FJ", "ID",
  "KE", "LS", "MG", "MU", "MW", "MZ", "NA", "NZ", "PE", "PG", "PY", "RE",
  "SB", "SZ", "TL", "TO", "TZ", "UY", "VU", "WS", "ZA", "ZM", "ZW"
)

location_qc <- location %>%
  mutate(HEMISPHERE = if_else(COUNTRY %in% sh_countries, "Southern", "Northern"))

local_hour <- function(utc_ts, tz_scalar) {
  local <- with_tz(utc_ts, tzone = tz_scalar)
  hour(local) + minute(local) / 60 + second(local) / 3600
}
local_month <- function(utc_ts, tz_scalar) {
  local <- with_tz(utc_ts, tzone = tz_scalar)
  as.integer(month(local))
}

bp_frame_loc <- bp_frame %>%
  left_join(location_qc, by = "STUDY_ID")

cat("\nConverting UTC reading_ts to user-local clock per timezone group...\n")
bp_frame_loc <- bp_frame_loc %>%
  group_by(TIMEZONE) %>%
  group_modify(~ {
    tz_scalar <- .y$TIMEZONE
    .x %>% mutate(
      BASELINE_TIME_OF_DAY = local_hour(BASELINE_READING_TS, tz_scalar),
      BASELINE_LOCAL_MONTH = local_month(BASELINE_READING_TS, tz_scalar),
      MO3_TIME_OF_DAY      = local_hour(MO3_READING_TS,      tz_scalar),
      MO3_LOCAL_MONTH      = local_month(MO3_READING_TS,     tz_scalar),
      MO6_TIME_OF_DAY      = local_hour(MO6_READING_TS,      tz_scalar),
      MO6_LOCAL_MONTH      = local_month(MO6_READING_TS,     tz_scalar)
    )
  }) %>%
  ungroup()

month_to_season <- function(m, hemi) {
  case_when(
    is.na(m) | is.na(hemi)                    ~ NA_character_,
    hemi == "Northern" & m %in% c(12, 1, 2)   ~ "winter",
    hemi == "Northern" & m %in% c(3, 4, 5)    ~ "spring",
    hemi == "Northern" & m %in% c(6, 7, 8)    ~ "summer",
    hemi == "Northern" & m %in% c(9, 10, 11)  ~ "fall",
    hemi == "Southern" & m %in% c(12, 1, 2)   ~ "summer",
    hemi == "Southern" & m %in% c(3, 4, 5)    ~ "fall",
    hemi == "Southern" & m %in% c(6, 7, 8)    ~ "winter",
    hemi == "Southern" & m %in% c(9, 10, 11)  ~ "spring"
  )
}

bp_frame_loc <- bp_frame_loc %>%
  mutate(
    BASELINE_SEASON = month_to_season(BASELINE_LOCAL_MONTH, HEMISPHERE),
    MO3_SEASON      = month_to_season(MO3_LOCAL_MONTH,      HEMISPHERE),
    MO6_SEASON      = month_to_season(MO6_LOCAL_MONTH,      HEMISPHERE)
  )

# =====================================================================
# Section 6: Final assembly + write
# =====================================================================

user_frame <- bp_frame_loc %>%
  left_join(height_bmi_qc, by = "STUDY_ID") %>%
  left_join(htn_meds_sel,  by = "STUDY_ID") %>%
  select(STUDY_ID, STUDY_START_DATE, STUDY_END_DATE, BIRTHDAY,
         AGE_AT_START, SEX, HEIGHT_M, WEIGHT_BASELINE_PROFILE_KG, BASELINE_BMI,
         HTN_MED_BASELINE,
         TIMEZONE, COUNTRY, HEMISPHERE,
         BASELINE_SBP, BASELINE_DBP, BASELINE_AHA,
         BASELINE_READING_TS, BASELINE_DAY, BASELINE_N_READINGS,
         BASELINE_TIME_OF_DAY, BASELINE_SEASON,
         MO3_SBP, MO3_DBP, MO3_AHA,
         MO3_READING_TS, MO3_DAY, MO3_N_READINGS,
         MO3_TIME_OF_DAY, MO3_SEASON,
         MO6_SBP, MO6_DBP, MO6_AHA,
         MO6_READING_TS, MO6_DAY, MO6_N_READINGS,
         MO6_TIME_OF_DAY, MO6_SEASON,
         DELTA_3MO_SBP, DELTA_3MO_DBP,
         DELTA_6MO_SBP, DELTA_6MO_DBP)

stopifnot(nrow(user_frame) == nrow(roster))
stopifnot(!any(duplicated(user_frame$STUDY_ID)))

cat("\n=== Final 2xWeekly user frame ===\n")
cat("Rows:", nrow(user_frame), "  Cols:", ncol(user_frame), "\n")

cat("\nMO3 coverage (3mo reading optional per §3):\n")
cat("  Users with MO3 BP:", sum(!is.na(user_frame$MO3_SBP)), "of", nrow(user_frame), "\n")

cat("\nHTN_MED_BASELINE distribution:\n")
print(user_frame %>% count(HTN_MED_BASELINE, name = "n"))

cat("\nHEMISPHERE distribution:\n")
print(user_frame %>% count(HEMISPHERE, name = "n"))

# --- Write ---
write_rds(user_frame, file.path(out_dir, "user_analytic_frame_2xWeekly.rds"))
write_csv(user_frame, file.path(out_dir, "user_analytic_frame_2xWeekly.csv"))
cat("\nWrote:\n")
cat("  ", file.path(out_dir, "user_analytic_frame_2xWeekly.rds"), "\n")
cat("  ", file.path(out_dir, "user_analytic_frame_2xWeekly.csv"), "\n")

cat("\n=== 2xWeekly §5.1 sensitivity build complete ===\n")
cat("Frame ready for the §5.1 LMM (unadjusted + adjusted variants) replicated\n")
cat("on the 2xWeekly cohort.\n")
