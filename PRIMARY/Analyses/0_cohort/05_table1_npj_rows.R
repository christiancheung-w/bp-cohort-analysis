# 05_table1_npj_rows.R
# Rows added to Table 1 for npj: antihypertensive medication status in four
# groups, and days from baseline to the 6-month reading. Same cohort as
# Table 1 (Elevated + Stage 1 + Stage 2, no sex restriction; N = 20,021).
#
# Medication groups (HTN_MED_BASELINE, Build/03_demographics_clinical.R):
#   stable_on_meds      first log >= 56 days before baseline
#   recent_initiator    first log 1-55 days before baseline
#   mid_study_initiator first log on/after baseline
#   never_logged        no log from 180 days before baseline to end of follow-up
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
# Output: Analyses/0_cohort/derived/results/table_1_npj_rows.csv

suppressPackageStartupMessages(library(tidyverse))

project_dir <- "./PRIMARY"
results_dir <- file.path(project_dir, "Analyses", "0_cohort", "derived", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

HYP <- c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")
u <- read_rds(file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")) %>%
  filter(BASELINE_AHA %in% HYP)
stopifnot(nrow(u) == 20021)

n_pct <- function(k, n) sprintf("%s (%.1f)", format(k, big.mark = ",", trim = TRUE), 100 * k / n)
med_iqr <- function(x) {
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  sprintf("%.0f [%.0f–%.0f]", q[2], q[1], q[3])
}
rows_for <- function(d) {
  n <- nrow(d)
  tibble(
    row = c("Antihypertensive medication: established use (first logged ≥8 weeks before baseline), n (%)",
            "Antihypertensive medication: recent initiation (<8 weeks before baseline), n (%)",
            "Antihypertensive medication: initiation during follow-up, n (%)",
            "Antihypertensive medication: none recorded, n (%)",
            "Days from baseline to 6-month reading, median [IQR]"),
    value = c(n_pct(sum(d$HTN_MED_BASELINE == "stable_on_meds"), n),
              n_pct(sum(d$HTN_MED_BASELINE == "recent_initiator"), n),
              n_pct(sum(d$HTN_MED_BASELINE == "mid_study_initiator"), n),
              n_pct(sum(d$HTN_MED_BASELINE == "never_logged"), n),
              med_iqr(d$MO6_DAY)))
}

out <- bind_rows(
  rows_for(u) %>% mutate(column = sprintf("Overall (%s)", format(nrow(u), big.mark = ","))),
  map_dfr(HYP, function(s) {
    d <- u %>% filter(BASELINE_AHA == s)
    rows_for(d) %>% mutate(column = sprintf("%s (%s)", sub(" Hypertension", "", s), format(nrow(d), big.mark = ",")))
  })) %>%
  pivot_wider(names_from = column, values_from = value)

write_csv(out, file.path(results_dir, "table_1_npj_rows.csv"))
print(as.data.frame(out), row.names = FALSE)
cat("\nWrote:", file.path(results_dir, "table_1_npj_rows.csv"), "\n")
