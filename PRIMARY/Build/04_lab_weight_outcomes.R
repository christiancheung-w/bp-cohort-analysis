# 04_lab_weight_outcomes.R
# Phase 1 / Slice 4 of the user-level frame build (see Build/BUILD_TODO.md).
#
# Inputs:
#   Roster CSV                                         — STUDY_ID universe
#   Query Results/2026-05-12 Pulls/WEIGHT              — windowed pre/post mean MEASURED weight + counts
#   Query Results/2026-05-12 Pulls/AL (labs)           — 10 analytes × {pre, post} × {value, date}
#
# Derives (one row per roster user, n = 28,290):
#
#   Weight (SAP §4.2 — uses MEASURED, not the profile snapshot from slice 3):
#     - WEIGHT_BASELINE_MEASURED_KG     (from WEIGHT.PRE_WEIGHT_KG;  bounds [30, 250] applied)
#     - WEIGHT_FOLLOWUP_MEASURED_KG     (from WEIGHT.POST_WEIGHT_KG; bounds [30, 250] applied)
#     - WEIGHT_BASELINE_MEASURED_COUNT  (from WEIGHT.PRE_WEIGHT_COUNT)
#     - WEIGHT_FOLLOWUP_MEASURED_COUNT  (from WEIGHT.POST_WEIGHT_COUNT)
#     - DELTA_WEIGHT_KG                 (FOLLOWUP - BASELINE on the measured values)
#
#   Labs (per SAP §4.2 — Δ for the 6 specified analytes only):
#     For ALL 10 analytes (TC, HDL, LDL, triglycerides, ApoB, glucose, HbA1c,
#     vitamin D, creatinine, eGFR): preserve PRE_<X>, POST_<X>, PRE_<X>_DATE, POST_<X>_DATE
#     from the source pull verbatim (no PRE/POST naming ambiguity — labs are
#     always measured).
#     Compute DELTA_<X> = POST - PRE for the 6 SAP-specified analytes:
#       DELTA_TOTAL_CHOLESTEROL, DELTA_HDL_CHOLESTEROL, DELTA_LDL_CHOLESTEROL,
#       DELTA_TRIGLYCERIDES, DELTA_GLUCOSE, DELTA_HBA1C.
#     ApoB / vitamin D / creatinine / eGFR: values + dates stored; Δ skipped per
#     strict SAP scope (can be added in a future slice if scope expands).
#
# Validation:
#   No roster gold-standard (these are new derivations). Pre-flight: row counts,
#   STUDY_ID alignment, no dups. Coverage: paired-completeness counts per analyte
#   (n with both pre AND post non-NA). Weight outliers flagged via bounds.

library(tidyverse)

project_dir <- "./PRIMARY"
data_dir    <- file.path(project_dir, "Query Results")
roster_dir  <- file.path(project_dir, "Roster")
build_dir   <- file.path(project_dir, "Build")
derived_dir <- file.path(build_dir, "derived")
pull_dir    <- file.path(data_dir, "2026-05-12 Pulls")

dir.create(derived_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read inputs ---
roster <- read_csv(file.path(roster_dir, "ROSTER_1X_2026-05-12.csv"))
weight <- read_csv(file.path(pull_dir,   "WEIGHT_1X_2026-05-12.csv"))
labs   <- read_csv(file.path(pull_dir,   "AL_1X_2026-05-12.csv"))

# --- Pre-flight ---
cat("Roster rows:", nrow(roster), "  users:", n_distinct(roster$STUDY_ID), "\n")
cat("Weight rows:", nrow(weight), "  users:", n_distinct(weight$STUDY_ID), "\n")
cat("Labs rows:  ", nrow(labs),   "  users:", n_distinct(labs$STUDY_ID),   "\n")

stopifnot(!any(duplicated(weight$STUDY_ID)))
stopifnot(!any(duplicated(labs$STUDY_ID)))

# WEIGHT and AL (labs) are conditional pulls — users only appear if they have
# device-measured weights or blood-panel results. Roster users without these
# measurements are absent from the pull and will land as NA after the left-join.
# Verify no infiltration (pull users ⊂ roster); substrate-completeness coverage
# reported below and per-analyte in the coverage block.
stopifnot(all(weight$STUDY_ID %in% roster$STUDY_ID))
stopifnot(all(labs$STUDY_ID   %in% roster$STUDY_ID))

cat("Weight pull coverage:", n_distinct(weight$STUDY_ID), "/", n_distinct(roster$STUDY_ID),
    sprintf("(%.1f%%)", 100 * n_distinct(weight$STUDY_ID) / n_distinct(roster$STUDY_ID)), "\n")
cat("Labs   pull coverage:", n_distinct(labs$STUDY_ID),   "/", n_distinct(roster$STUDY_ID),
    sprintf("(%.1f%%)", 100 * n_distinct(labs$STUDY_ID)   / n_distinct(roster$STUDY_ID)), "\n")

# --- Weight: outlier bounds + derive ---
n_pre_wt_bad  <- sum(!is.na(weight$PRE_WEIGHT_KG)  & (weight$PRE_WEIGHT_KG  < 30 | weight$PRE_WEIGHT_KG  > 250))
n_post_wt_bad <- sum(!is.na(weight$POST_WEIGHT_KG) & (weight$POST_WEIGHT_KG < 30 | weight$POST_WEIGHT_KG > 250))
cat("\nWeight outlier check ([30, 250] kg bounds):\n")
cat("  PRE_WEIGHT_KG  out-of-bounds (NA'd):", n_pre_wt_bad,  "\n")
cat("  POST_WEIGHT_KG out-of-bounds (NA'd):", n_post_wt_bad, "\n")

weight <- weight %>%
  mutate(
    WEIGHT_BASELINE_MEASURED_KG    = if_else(PRE_WEIGHT_KG  >= 30 & PRE_WEIGHT_KG  <= 250, PRE_WEIGHT_KG,  NA_real_),
    WEIGHT_FOLLOWUP_MEASURED_KG    = if_else(POST_WEIGHT_KG >= 30 & POST_WEIGHT_KG <= 250, POST_WEIGHT_KG, NA_real_),
    WEIGHT_BASELINE_MEASURED_COUNT = PRE_WEIGHT_COUNT,
    WEIGHT_FOLLOWUP_MEASURED_COUNT = POST_WEIGHT_COUNT,
    DELTA_WEIGHT_KG                = WEIGHT_FOLLOWUP_MEASURED_KG - WEIGHT_BASELINE_MEASURED_KG
  )

# --- Labs: compute Δ for 6 SAP-specified analytes + derive non-HDL cholesterol ---
# Non-HDL = TC − HDL (per SAP §6.1 +Δlabs variant predictor set). Slice 98 will
# re-derive non-HDL from QC'd PRE/POST values after unit conversion + bounds.
labs <- labs %>%
  mutate(
    DELTA_TOTAL_CHOLESTEROL = POST_TOTAL_CHOLESTEROL - PRE_TOTAL_CHOLESTEROL,
    DELTA_HDL_CHOLESTEROL   = POST_HDL_CHOLESTEROL   - PRE_HDL_CHOLESTEROL,
    DELTA_LDL_CHOLESTEROL   = POST_LDL_CHOLESTEROL   - PRE_LDL_CHOLESTEROL,
    DELTA_TRIGLYCERIDES     = POST_TRIGLYCERIDES     - PRE_TRIGLYCERIDES,
    DELTA_GLUCOSE           = POST_GLUCOSE           - PRE_GLUCOSE,
    DELTA_HBA1C             = POST_HBA1C             - PRE_HBA1C,
    PRE_NON_HDL_CHOLESTEROL   = PRE_TOTAL_CHOLESTEROL  - PRE_HDL_CHOLESTEROL,
    POST_NON_HDL_CHOLESTEROL  = POST_TOTAL_CHOLESTEROL - POST_HDL_CHOLESTEROL,
    DELTA_NON_HDL_CHOLESTEROL = POST_NON_HDL_CHOLESTEROL - PRE_NON_HDL_CHOLESTEROL
  )

# --- Assemble per-user output ---
lab_weight <- roster %>%
  select(STUDY_ID) %>%
  left_join(
    weight %>% select(STUDY_ID,
                      WEIGHT_BASELINE_MEASURED_KG, WEIGHT_BASELINE_MEASURED_COUNT,
                      WEIGHT_FOLLOWUP_MEASURED_KG, WEIGHT_FOLLOWUP_MEASURED_COUNT,
                      DELTA_WEIGHT_KG),
    by = "STUDY_ID"
  ) %>%
  left_join(
    labs %>% select(STUDY_ID,
                    starts_with("PRE_"), starts_with("POST_"),
                    starts_with("DELTA_")),
    by = "STUDY_ID"
  )

stopifnot(nrow(lab_weight) == nrow(roster))
stopifnot(!any(duplicated(lab_weight$STUDY_ID)))

# --- Coverage / summary prints ---

cat("\n=== Weight (SAP §4.2) ===\n")
cat("  Users with WEIGHT_BASELINE_MEASURED_KG (post-bounds):", sum(!is.na(lab_weight$WEIGHT_BASELINE_MEASURED_KG)), "\n")
cat("  Users with WEIGHT_FOLLOWUP_MEASURED_KG (post-bounds):", sum(!is.na(lab_weight$WEIGHT_FOLLOWUP_MEASURED_KG)), "\n")
cat("  Users with paired weights (DELTA_WEIGHT_KG computable):", sum(!is.na(lab_weight$DELTA_WEIGHT_KG)), "\n")
cat("  DELTA_WEIGHT_KG distribution:\n")
print(summary(lab_weight$DELTA_WEIGHT_KG))

cat("\n=== Labs: paired-completeness per analyte ===\n")
analytes <- c("TOTAL_CHOLESTEROL", "HDL_CHOLESTEROL", "LDL_CHOLESTEROL",
              "NON_HDL_CHOLESTEROL",
              "TRIGLYCERIDES", "APOLIPOPROTEIN_B", "GLUCOSE", "HBA1C",
              "VITAMIN_D", "CREATININE", "EGFR")
for (a in analytes) {
  pre_col   <- paste0("PRE_",  a)
  post_col  <- paste0("POST_", a)
  delta_col <- paste0("DELTA_", a)
  n_pre     <- sum(!is.na(lab_weight[[pre_col]]))
  n_post    <- sum(!is.na(lab_weight[[post_col]]))
  n_paired  <- sum(!is.na(lab_weight[[pre_col]]) & !is.na(lab_weight[[post_col]]))
  has_delta <- delta_col %in% names(lab_weight)
  tag       <- if (has_delta) "  (Δ computed)" else "  (no Δ — out of SAP scope)"
  cat(sprintf("  %-22s  pre = %5d  post = %5d  paired = %5d%s\n",
              a, n_pre, n_post, n_paired, tag))
}

cat("\n=== Δ-analyte distributions (SAP §4.2 — 6 analytes) ===\n")
for (a in c("TOTAL_CHOLESTEROL", "HDL_CHOLESTEROL", "LDL_CHOLESTEROL",
            "NON_HDL_CHOLESTEROL",
            "TRIGLYCERIDES", "GLUCOSE", "HBA1C")) {
  delta_col <- paste0("DELTA_", a)
  cat("  ", delta_col, ":\n", sep = "")
  print(summary(lab_weight[[delta_col]]))
}

cat("\nlab_weight frame:", nrow(lab_weight), "rows,", ncol(lab_weight), "cols\n")

# --- Write ---
write_rds(lab_weight, file.path(derived_dir, "04_lab_weight_outcomes.rds"))
write_csv(lab_weight, file.path(derived_dir, "04_lab_weight_outcomes.csv"))
cat("\nWrote:\n")
cat("  ", file.path(derived_dir, "04_lab_weight_outcomes.rds"), "\n")
cat("  ", file.path(derived_dir, "04_lab_weight_outcomes.csv"), "\n")
