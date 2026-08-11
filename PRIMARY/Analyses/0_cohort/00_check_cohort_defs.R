# 00_check_cohort_defs.R
# -----------------------------------------------------------------------------
# Slice-1 verification for the cohort-chapter refactor (branch: review-0_cohort).
# Sources cohort_defs.R, confirms each cohort's n reproduces the frozen Table 1
# numbers, and prints the hypertensive characteristic table with the new
# between-stratum P-value column. Writes NO output files — console only.
#
# Run with working directory = project root (the dID_Analysis.Rproj folder).
# -----------------------------------------------------------------------------

library(tidyverse)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")

frame <- read_rds("./PRIMARY/Build/derived/qc/user_analytic_frame_1X.rds")
cat("Loaded user frame:", nrow(frame), "rows x", ncol(frame), "cols\n\n")

cat("=== Cohort n's ===\n")
cat("(frozen reference: hypertensive = 20,021 ; full = 28,290 ; htn_naive = expect ~16,2xx)\n")
for (nm in names(COHORTS)) {
  d <- apply_cohort(frame, nm)
  cat(sprintf("\n  %-12s  n = %s   [%s]\n", nm, format(nrow(d), big.mark = ","),
              COHORTS[[nm]]$label))
  print(d %>% count(BASELINE_AHA, name = "n"))
}

cat("\n\n=== Hypertensive characteristic table WITH between-stratum P-values ===\n")
cat("(descriptive cells should match the frozen table_1_paper.md exactly;\n")
cat(" the P-value column is the only new content)\n\n")
ht <- apply_cohort(frame, "hypertensive")
print(make_char_table(ht, COHORTS$hypertensive$strata), n = Inf, width = Inf)

cat("\n=== check done ===\n")
