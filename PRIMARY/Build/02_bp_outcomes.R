# 02_bp_outcomes.R
# Derive BP-based outcome columns from the per-user 4-timepoint frame.
# Phase 1 / Slice 2 of the user-level frame build (see Build/BUILD_TODO.md).
#
# Inputs:
#   Build/derived/01_bp_timepoints.rds  — per-user BASELINE / MO3 / MO6 / MO12 readings
#   Roster CSV                          — for AHA gold-standard validation
#
# Derives (one row per roster user, n = 28,290):
#   - Per-timepoint AHA classification: BASELINE_AHA, MO3_AHA, MO6_AHA, MO12_AHA
#       using the 2026-05-14 v2 continuous-comparison patch logic:
#         SBP >= 140 OR DBP >= 90    -> "Stage 2 Hypertension"
#         SBP >= 130 OR DBP >= 80    -> "Stage 1 Hypertension"
#         SBP >= 120 AND DBP <  80   -> "Elevated"
#         else                       -> "Normal"
#   - Per-timepoint Δ from baseline (SBP and DBP):
#         DELTA_3MO_SBP, DELTA_3MO_DBP,
#         DELTA_6MO_SBP, DELTA_6MO_DBP,
#         DELTA_12MO_SBP, DELTA_12MO_DBP
#   - SBP responder at 6mo (binary, SAP §6.2):
#         SBP_RESPONDER_6MO = (DELTA_6MO_SBP <= -5)
#   - 3-bin BP-change category at 6mo (§5.2) and 12mo (§5.5):
#         improved (Δ <= -5) / stable (-5 < Δ < +5) / worsened (Δ >= +5)
#         separately for SBP and DBP:
#         BPCAT_6MO_SBP,  BPCAT_6MO_DBP,
#         BPCAT_12MO_SBP, BPCAT_12MO_DBP   (12mo on the 12mo subset only)
#
# Validation (gold-standard against patched roster):
#   - BASELINE_AHA exactly matches roster.BASELINE_AHA_CLASSIFICATION  (stopifnot 0 disagreements)
#   - MO6_AHA      exactly matches roster.END_AHA_CLASSIFICATION       (stopifnot 0 disagreements)
#   - Per-stratum n's printed side-by-side (R-derived vs. roster) at both timepoints.

library(tidyverse)

project_dir <- "./PRIMARY"
roster_dir  <- file.path(project_dir, "Roster")
build_dir   <- file.path(project_dir, "Build")
derived_dir <- file.path(build_dir, "derived")

bpcat_levels <- c("improved", "stable", "worsened")

# --- Read inputs ---
bp_timepoints <- read_rds(file.path(derived_dir, "01_bp_timepoints.rds"))
roster        <- read_csv(file.path(roster_dir, "ROSTER_1X_2026-05-12.csv"))

cat("bp_timepoints rows:", nrow(bp_timepoints), " cols:", ncol(bp_timepoints), "\n")
cat("roster rows:", nrow(roster), " cols:", ncol(roster), "\n")
stopifnot(nrow(bp_timepoints) == nrow(roster))
stopifnot(setequal(bp_timepoints$STUDY_ID, roster$STUDY_ID))

# --- AHA classification (2026-05-14 v2 continuous-comparison patch) ---
# Order of case_when branches is load-bearing: Stage 2 catches first (>=140 SBP
# OR >=90 DBP), then Stage 1 (>=130 SBP OR >=80 DBP), then Elevated (>=120 SBP
# AND <80 DBP — the DBP < 80 clause is redundant given the prior branches'
# negation but kept explicit to mirror the SQL patch). Returns NA when either
# SBP or DBP is NA (MO3 / MO12 will be NA for users without that timepoint).
classify_aha <- function(sbp, dbp) {
  case_when(
    is.na(sbp) | is.na(dbp)    ~ NA_character_,
    sbp >= 140 | dbp >= 90     ~ "Stage 2 Hypertension",
    sbp >= 130 | dbp >= 80     ~ "Stage 1 Hypertension",
    sbp >= 120 & dbp <  80     ~ "Elevated",
    TRUE                       ~ "Normal"
  )
}

bp_outcomes <- bp_timepoints %>%
  mutate(
    BASELINE_AHA = classify_aha(BASELINE_SBP, BASELINE_DBP),
    MO3_AHA      = classify_aha(MO3_SBP,      MO3_DBP),
    MO6_AHA      = classify_aha(MO6_SBP,      MO6_DBP),
    MO12_AHA     = classify_aha(MO12_SBP,     MO12_DBP)
  )

# --- Per-timepoint deltas from baseline ---
bp_outcomes <- bp_outcomes %>%
  mutate(
    DELTA_3MO_SBP  = MO3_SBP  - BASELINE_SBP,
    DELTA_3MO_DBP  = MO3_DBP  - BASELINE_DBP,
    DELTA_6MO_SBP  = MO6_SBP  - BASELINE_SBP,
    DELTA_6MO_DBP  = MO6_DBP  - BASELINE_DBP,
    DELTA_12MO_SBP = MO12_SBP - BASELINE_SBP,
    DELTA_12MO_DBP = MO12_DBP - BASELINE_DBP
  )

# --- SBP responder at 6mo (SAP §6.2) ---
bp_outcomes <- bp_outcomes %>%
  mutate(SBP_RESPONDER_6MO = DELTA_6MO_SBP <= -5)

# --- 3-bin BP-change category (SAP §5.2 / §5.5) ---
bin_bp_change <- function(delta) {
  factor(
    case_when(
      is.na(delta) ~ NA_character_,
      delta <= -5  ~ "improved",
      delta >=  5  ~ "worsened",
      TRUE         ~ "stable"
    ),
    levels = bpcat_levels
  )
}

bp_outcomes <- bp_outcomes %>%
  mutate(
    BPCAT_6MO_SBP  = bin_bp_change(DELTA_6MO_SBP),
    BPCAT_6MO_DBP  = bin_bp_change(DELTA_6MO_DBP),
    BPCAT_12MO_SBP = bin_bp_change(DELTA_12MO_SBP),
    BPCAT_12MO_DBP = bin_bp_change(DELTA_12MO_DBP)
  )

# --- Validation: BASELINE_AHA vs. roster.BASELINE_AHA_CLASSIFICATION ---
check_base_aha <- bp_outcomes %>%
  select(STUDY_ID, BASELINE_AHA, BASELINE_SBP, BASELINE_DBP) %>%
  inner_join(
    roster %>% select(STUDY_ID, BASELINE_AHA_ROSTER = BASELINE_AHA_CLASSIFICATION),
    by = "STUDY_ID"
  )

n_base_aha_disagree  <- sum(check_base_aha$BASELINE_AHA != check_base_aha$BASELINE_AHA_ROSTER, na.rm = TRUE)
n_base_aha_na_r      <- sum(is.na(check_base_aha$BASELINE_AHA))
n_base_aha_na_roster <- sum(is.na(check_base_aha$BASELINE_AHA_ROSTER))

cat("\nBaseline AHA validation vs. roster:\n")
cat("  Disagreements (row-by-row):", n_base_aha_disagree, "\n")
cat("  NA in R-derived BASELINE_AHA:", n_base_aha_na_r, "\n")
cat("  NA in roster BASELINE_AHA_CLASSIFICATION:", n_base_aha_na_roster, "\n")
cat("  Per-stratum n's — roster (BASELINE_AHA_CLASSIFICATION):\n")
print(roster %>% count(BASELINE_AHA_CLASSIFICATION, name = "n_roster"))
cat("  Per-stratum n's — R-derived (BASELINE_AHA):\n")
print(bp_outcomes %>% count(BASELINE_AHA, name = "n_R"))
if (n_base_aha_disagree > 0) {
  cat("  First few disagreements (SBP, DBP shown for inspection):\n")
  print(head(check_base_aha %>%
               filter(BASELINE_AHA != BASELINE_AHA_ROSTER), 15))
}
stopifnot(n_base_aha_disagree == 0)
stopifnot(n_base_aha_na_r == 0)
stopifnot(n_base_aha_na_roster == 0)

# --- Validation: MO6_AHA vs. roster.END_AHA_CLASSIFICATION ---
check_end_aha <- bp_outcomes %>%
  select(STUDY_ID, MO6_AHA, MO6_SBP, MO6_DBP) %>%
  inner_join(
    roster %>% select(STUDY_ID, END_AHA_ROSTER = END_AHA_CLASSIFICATION),
    by = "STUDY_ID"
  )

n_end_aha_disagree  <- sum(check_end_aha$MO6_AHA != check_end_aha$END_AHA_ROSTER, na.rm = TRUE)
n_end_aha_na_r      <- sum(is.na(check_end_aha$MO6_AHA))
n_end_aha_na_roster <- sum(is.na(check_end_aha$END_AHA_ROSTER))

cat("\n6mo AHA validation vs. roster (END_AHA_CLASSIFICATION):\n")
cat("  Disagreements (row-by-row):", n_end_aha_disagree, "\n")
cat("  NA in R-derived MO6_AHA:", n_end_aha_na_r, "\n")
cat("  NA in roster END_AHA_CLASSIFICATION:", n_end_aha_na_roster, "\n")
cat("  Per-stratum n's — roster (END_AHA_CLASSIFICATION):\n")
print(roster %>% count(END_AHA_CLASSIFICATION, name = "n_roster"))
cat("  Per-stratum n's — R-derived (MO6_AHA):\n")
print(bp_outcomes %>% count(MO6_AHA, name = "n_R"))
if (n_end_aha_disagree > 0) {
  cat("  First few disagreements (MO6_SBP, MO6_DBP shown for inspection):\n")
  print(head(check_end_aha %>%
               filter(MO6_AHA != END_AHA_ROSTER), 15))
}
stopifnot(n_end_aha_disagree == 0)
stopifnot(n_end_aha_na_r == 0)
stopifnot(n_end_aha_na_roster == 0)

# --- Coverage / summary prints (new derivations, no roster gold standard) ---
cat("\n3mo AHA coverage (new — not in roster):\n")
print(bp_outcomes %>% count(MO3_AHA, name = "n"))
cat("\n12mo AHA coverage (new — not in roster):\n")
print(bp_outcomes %>% count(MO12_AHA, name = "n"))

cat("\nSBP responder at 6mo (ΔSBP <= -5):\n")
print(bp_outcomes %>% count(SBP_RESPONDER_6MO, name = "n"))

cat("\n6mo 3-bin BP-change category:\n")
cat("  SBP (BPCAT_6MO_SBP):\n")
print(bp_outcomes %>% count(BPCAT_6MO_SBP, name = "n"))
cat("  DBP (BPCAT_6MO_DBP):\n")
print(bp_outcomes %>% count(BPCAT_6MO_DBP, name = "n"))

cat("\n12mo 3-bin BP-change category (12mo subset):\n")
cat("  SBP (BPCAT_12MO_SBP):\n")
print(bp_outcomes %>% count(BPCAT_12MO_SBP, name = "n"))
cat("  DBP (BPCAT_12MO_DBP):\n")
print(bp_outcomes %>% count(BPCAT_12MO_DBP, name = "n"))

# --- Assemble outcomes-only frame and write ---
bp_outcomes_out <- bp_outcomes %>%
  select(STUDY_ID,
         BASELINE_AHA, MO3_AHA, MO6_AHA, MO12_AHA,
         DELTA_3MO_SBP,  DELTA_3MO_DBP,
         DELTA_6MO_SBP,  DELTA_6MO_DBP,
         DELTA_12MO_SBP, DELTA_12MO_DBP,
         SBP_RESPONDER_6MO,
         BPCAT_6MO_SBP,  BPCAT_6MO_DBP,
         BPCAT_12MO_SBP, BPCAT_12MO_DBP)

stopifnot(nrow(bp_outcomes_out) == nrow(roster))
stopifnot(!any(duplicated(bp_outcomes_out$STUDY_ID)))

cat("\nbp_outcomes frame:", nrow(bp_outcomes_out), "rows,",
    ncol(bp_outcomes_out), "cols\n")

write_rds(bp_outcomes_out, file.path(derived_dir, "02_bp_outcomes.rds"))
write_csv(bp_outcomes_out, file.path(derived_dir, "02_bp_outcomes.csv"))
cat("\nWrote:\n")
cat("  ", file.path(derived_dir, "02_bp_outcomes.rds"), "\n")
cat("  ", file.path(derived_dir, "02_bp_outcomes.csv"), "\n")
