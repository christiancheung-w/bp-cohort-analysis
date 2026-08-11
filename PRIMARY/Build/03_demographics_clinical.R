# 03_demographics_clinical.R
# Phase 1 / Slice 3 of the user-level frame build (see Build/BUILD_TODO.md).
#
# Inputs:
#   Roster CSV                                  — AGE_AT_START, GENDER
#   Query Results/2026-05-13 Pulls/HEIGHT_BMI   — HEIGHT, WEIGHT, BMI snapshot at STUDY_START_DATE
#   Query Results/2026-05-14 Pulls/HTN_MEDS     — BP-medication journal-tracker fields + BP_MED_STATUS
#
# Derives (one row per roster user, n = 28,290):
#   - AGE_AT_START         passthrough from roster
#   - SEX                  passthrough from roster.GENDER (lowercase character, 4 levels in data:
#                          male / female / non-binary / "i prefer not to choose").
#                          Downstream binarization (if needed) is an analysis-stage decision.
#   - HEIGHT_M             from HEIGHT_BMI.HEIGHT, NA-out if outside [1.4, 2.2] m
#   - WEIGHT_BASELINE_PROFILE_KG   from HEIGHT_BMI.WEIGHT, NA-out if outside [30, 250] kg
#                          (snapshot at STUDY_START_DATE; distinct from the windowed PRE_WEIGHT_KG
#                          aggregate that comes from the WEIGHT pull for Δweight)
#   - BASELINE_BMI         from HEIGHT_BMI.BMI, NA-out if outside [12, 70] kg/m²
#   - HTN_MED_BASELINE     from HTN_MEDS.BP_MED_STATUS (4 levels: stable_on_meds /
#                          recent_initiator / mid_study_initiator / never_logged)
#
# Validation (gold-standard, fail-fast):
#   BP_MED_STATUS_R recomputed from raw HTN_MEDS fields per SAP convention
#   (boundary at 56 days pre-study) exactly reproduces source BP_MED_STATUS.
#   Pattern mirrors build_analysis_frame_1X.R:241-253.

library(tidyverse)

project_dir   <- "./PRIMARY"
data_dir      <- file.path(project_dir, "Query Results")
roster_dir    <- file.path(project_dir, "Roster")
build_dir     <- file.path(project_dir, "Build")
derived_dir   <- file.path(build_dir, "derived")
pull_dir_0513 <- file.path(data_dir, "2026-05-13 Pulls")
pull_dir_0514 <- file.path(data_dir, "2026-05-14 Pulls")

dir.create(derived_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read inputs ---
roster     <- read_csv(file.path(roster_dir,    "ROSTER_1X_2026-05-12.csv"))
height_bmi <- read_csv(file.path(pull_dir_0513, "HEIGHT_BMI_1X_2026-05-13.csv"))
htn_meds   <- read_csv(file.path(pull_dir_0514, "HTN_MEDS_1X_2026-05-14.csv"))

# --- Pre-flight ---
cat("Roster rows:     ", nrow(roster),     "  users:", n_distinct(roster$STUDY_ID),     "\n")
cat("HEIGHT_BMI rows: ", nrow(height_bmi), "  users:", n_distinct(height_bmi$STUDY_ID), "\n")
cat("HTN_MEDS rows:   ", nrow(htn_meds),   "  users:", n_distinct(htn_meds$STUDY_ID),   "\n")

stopifnot(!any(duplicated(roster$STUDY_ID)))
stopifnot(!any(duplicated(height_bmi$STUDY_ID)))
stopifnot(!any(duplicated(htn_meds$STUDY_ID)))

# 1X cohort is the universe — both substrates should align exactly with roster.
stopifnot(setequal(height_bmi$STUDY_ID, roster$STUDY_ID))
stopifnot(setequal(htn_meds$STUDY_ID,   roster$STUDY_ID))

# --- Validation: BP_MED_STATUS_R recomputation vs. source BP_MED_STATUS ---
# Convention (SAP): stable_on_meds if first log was at least 56 days before
# STUDY_START_DATE; recent_initiator if within 56 days pre-study;
# mid_study_initiator if on/after STUDY_START_DATE; never_logged if no log
# in (STUDY_START_DATE - 180d, STUDY_END_DATE).
htn_meds <- htn_meds %>%
  mutate(
    BP_MED_STATUS_R = case_when(
      is.na(BP_MED_FIRST_LOG_DATE)                  ~ "never_logged",
      BP_MED_DAYS_FIRST_LOG_VS_STUDY_START <= -56   ~ "stable_on_meds",
      BP_MED_DAYS_FIRST_LOG_VS_STUDY_START <   0    ~ "recent_initiator",
      BP_MED_DAYS_FIRST_LOG_VS_STUDY_START >=  0    ~ "mid_study_initiator"
    )
  )

n_bp_med_disagree <- sum(htn_meds$BP_MED_STATUS != htn_meds$BP_MED_STATUS_R, na.rm = TRUE)
n_bp_med_na_r     <- sum(is.na(htn_meds$BP_MED_STATUS_R))
n_bp_med_na_src   <- sum(is.na(htn_meds$BP_MED_STATUS))

cat("\nBP_MED_STATUS validation (R recomputation vs. source):\n")
cat("  Disagreements (row-by-row):", n_bp_med_disagree, "\n")
cat("  NA in R-derived:", n_bp_med_na_r, "\n")
cat("  NA in source:   ", n_bp_med_na_src, "\n")
cat("  Per-stratum n's — source (BP_MED_STATUS):\n")
print(htn_meds %>% count(BP_MED_STATUS, name = "n_source"))
cat("  Per-stratum n's — R-derived (BP_MED_STATUS_R):\n")
print(htn_meds %>% count(BP_MED_STATUS_R, name = "n_R"))
stopifnot(identical(htn_meds$BP_MED_STATUS, htn_meds$BP_MED_STATUS_R))
cat("BP_MED_STATUS reproduced exactly from raw HTN_MEDS fields.\n")

# --- Anthropometric NA-outlier handling ---
# Bounds chosen as defensible adult ranges; mirror build_analysis_frame_1X.R:200-211.
n_height_bad <- sum(!is.na(height_bmi$HEIGHT) & (height_bmi$HEIGHT < 1.4 | height_bmi$HEIGHT > 2.2))
n_weight_bad <- sum(!is.na(height_bmi$WEIGHT) & (height_bmi$WEIGHT < 30  | height_bmi$WEIGHT > 250))
n_bmi_bad    <- sum(!is.na(height_bmi$BMI)    & (height_bmi$BMI    < 12  | height_bmi$BMI    > 70))
cat("\nAnthropometric out-of-range values NA'd — HEIGHT:", n_height_bad,
    " WEIGHT:", n_weight_bad, " BMI:", n_bmi_bad, "\n")

height_bmi <- height_bmi %>%
  mutate(
    HEIGHT_M           = if_else(HEIGHT >= 1.4 & HEIGHT <= 2.2, HEIGHT, NA_real_),
    WEIGHT_BASELINE_PROFILE_KG = if_else(WEIGHT >= 30  & WEIGHT <= 250, WEIGHT, NA_real_),
    BASELINE_BMI       = if_else(BMI    >= 12  & BMI    <= 70,  BMI,    NA_real_)
  )

# --- Assemble §12.2 frame ---
demo_clinical <- roster %>%
  select(STUDY_ID, AGE_AT_START, GENDER) %>%
  rename(SEX = GENDER) %>%
  left_join(
    height_bmi %>% select(STUDY_ID, HEIGHT_M, WEIGHT_BASELINE_PROFILE_KG, BASELINE_BMI),
    by = "STUDY_ID"
  ) %>%
  left_join(
    htn_meds %>% select(STUDY_ID, HTN_MED_BASELINE = BP_MED_STATUS),
    by = "STUDY_ID"
  )

stopifnot(nrow(demo_clinical) == nrow(roster))
stopifnot(!any(duplicated(demo_clinical$STUDY_ID)))

# --- Coverage / summary prints ---
cat("\ndemo_clinical frame:", nrow(demo_clinical), "rows,", ncol(demo_clinical), "cols\n")

cat("\nAge distribution (AGE_AT_START):\n")
print(summary(demo_clinical$AGE_AT_START))

cat("\nSex distribution:\n")
print(demo_clinical %>% count(SEX, name = "n"))

cat("\nAnthropometric coverage (after NA-out):\n")
cat("  Non-NA HEIGHT_M:           ", sum(!is.na(demo_clinical$HEIGHT_M)),
    " / ", nrow(demo_clinical), "\n", sep = "")
cat("  Non-NA WEIGHT_BASELINE_PROFILE_KG: ", sum(!is.na(demo_clinical$WEIGHT_BASELINE_PROFILE_KG)),
    " / ", nrow(demo_clinical), "\n", sep = "")
cat("  Non-NA BASELINE_BMI:       ", sum(!is.na(demo_clinical$BASELINE_BMI)),
    " / ", nrow(demo_clinical), "\n", sep = "")

cat("\nBASELINE_BMI distribution (after NA-out):\n")
print(summary(demo_clinical$BASELINE_BMI))

cat("\nHTN_MED_BASELINE distribution:\n")
print(demo_clinical %>% count(HTN_MED_BASELINE, name = "n"))

# --- Write ---
write_rds(demo_clinical, file.path(derived_dir, "03_demographics_clinical.rds"))
write_csv(demo_clinical, file.path(derived_dir, "03_demographics_clinical.csv"))
cat("\nWrote:\n")
cat("  ", file.path(derived_dir, "03_demographics_clinical.rds"), "\n")
cat("  ", file.path(derived_dir, "03_demographics_clinical.csv"), "\n")
