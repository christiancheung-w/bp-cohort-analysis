# 99_assemble_user_frame.R
# Phase 1 / Slice 99: final assembly of the canonical user-level analytic frame.
#
# Reads all 8 QC'd substrates from Build/derived/qc/ and left-joins them onto a
# roster-anchor base (STUDY_ID, STUDY_START_DATE, STUDY_END_DATE, BIRTHDAY).
# Mirrors the build_analysis_frame_1X.R pattern:
#   read → pre-flight → verify-shared → drop → left-join → invariant.
#
# Output: Build/derived/qc/user_analytic_frame_1X.{rds, csv} — the user-level
# (one row per user) substrate for all §5+ SAP analyses except §7.2.
# (A sibling BP-reading-level frame for §7.2 will come from a separate Phase 2 build.)

library(tidyverse)

project_dir <- "./PRIMARY"
roster_dir  <- file.path(project_dir, "Roster")
derived_dir <- file.path(project_dir, "Build", "derived")
qc_dir      <- file.path(derived_dir, "qc")

# --- Read roster (for anchors not in any substrate) + QC'd substrates ---
roster <- read_csv(file.path(roster_dir, "ROSTER_1X_2026-05-12.csv"))

bp_tp  <- read_rds(file.path(qc_dir, "01_bp_timepoints.rds"))
bp_out <- read_rds(file.path(qc_dir, "02_bp_outcomes.rds"))
demo   <- read_rds(file.path(qc_dir, "03_demographics_clinical.rds"))
lab_wt <- read_rds(file.path(qc_dir, "04_lab_weight_outcomes.rds"))
whoop  <- read_rds(file.path(qc_dir, "05_whoop_window_means.rds"))
le8tn  <- read_rds(file.path(qc_dir, "06_le8_tobacco_nicotine.rds"))
eng    <- read_rds(file.path(qc_dir, "07_engagement.rds"))
lmm    <- read_rds(file.path(qc_dir, "08_lmm_covariates.rds"))

# --- Pre-flight ---
substrates <- list(
  "01_bp_timepoints"          = bp_tp,
  "02_bp_outcomes"            = bp_out,
  "03_demographics_clinical"  = demo,
  "04_lab_weight_outcomes"    = lab_wt,
  "05_whoop_window_means"     = whoop,
  "06_le8_tobacco_nicotine"   = le8tn,
  "07_engagement"             = eng,
  "08_lmm_covariates"         = lmm
)

cat("Roster:", nrow(roster), "users\n")
cat("\nSubstrate pre-flight:\n")
for (nm in names(substrates)) {
  df <- substrates[[nm]]
  stopifnot(!any(duplicated(df$STUDY_ID)))
  stopifnot(setequal(df$STUDY_ID, roster$STUDY_ID))
  cat(sprintf("  %-30s  %5d rows × %3d cols  ✓ (cohort matches roster)\n",
              nm, nrow(df), ncol(df)))
}

# --- Build anchor frame (roster cols not in substrates) ---
# Slice 1 already carries STUDY_START_DATE; we add STUDY_END_DATE and BIRTHDAY.
# BASELINE_SBP/DBP, FOLLOWUP_SBP/DBP, DELTA_SBP/DBP, AHA_CLASSIFICATION etc. are
# the legacy roster derivations — we pull the QC'd substrate versions instead.
anchor <- roster %>%
  select(STUDY_ID, STUDY_END_DATE, BIRTHDAY)

# --- Sequential left-join, dropping shared columns from the right side ---
join_with_dedupe <- function(left, right, right_name) {
  shared <- setdiff(intersect(names(left), names(right)), "STUDY_ID")
  if (length(shared) > 0) {
    cat(sprintf("  Dropping %d already-present cols from %s: %s\n",
                length(shared), right_name, paste(shared, collapse = ", ")))
    right <- right %>% select(-all_of(shared))
  }
  result <- left %>% left_join(right, by = "STUDY_ID")
  stopifnot(nrow(result) == nrow(left))
  result
}

cat("\n=== Assembling analytic frame ===\n")
frame <- anchor
cat(sprintf("  Start: anchor frame  %5d rows × %3d cols\n", nrow(frame), ncol(frame)))

for (nm in names(substrates)) {
  frame <- join_with_dedupe(frame, substrates[[nm]], nm)
  cat(sprintf("  + %-30s  →   %5d rows × %3d cols\n", nm, nrow(frame), ncol(frame)))
}

stopifnot(nrow(frame) == nrow(roster))
stopifnot(!any(duplicated(frame$STUDY_ID)))

# --- Summary prints ---
cat(sprintf("\nFinal analytic frame: %d rows × %d cols\n", nrow(frame), ncol(frame)))

cat("\nColumn class summary:\n")
print(tibble(column = names(frame), class = sapply(frame, function(x) class(x)[1])) %>%
        count(class, name = "n_cols") %>% arrange(desc(n_cols)))

cat("\nSpot-check: first 5 users, key columns:\n")
print(frame %>%
        select(STUDY_ID, STUDY_START_DATE, STUDY_END_DATE,
               AGE_AT_START, SEX,
               BASELINE_AHA, MO6_AHA,
               DELTA_6MO_SBP, SBP_RESPONDER_6MO,
               TIMEZONE, COUNTRY, HEMISPHERE) %>%
        head(5))

# --- Write ---
write_rds(frame, file.path(qc_dir, "user_analytic_frame_1X.rds"))
write_csv(frame, file.path(qc_dir, "user_analytic_frame_1X.csv"))
cat("\nWrote:\n")
cat("  ", file.path(qc_dir, "user_analytic_frame_1X.rds"), "\n")
cat("  ", file.path(qc_dir, "user_analytic_frame_1X.csv"), "\n")

cat("\n=== Phase 1 build complete ===\n")
cat("The canonical user-level analytic frame is ready for §5+ SAP analyses (except §7.2 — that\n")
cat("requires the BP-reading-level long frame from the separate Phase 2 build).\n")
