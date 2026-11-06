# 07_table1_s1_npj_20003.R
# Table 1 and Table S1 (completeness) for the npj analytic cohort: Elevated,
# Stage 1 and Stage 2 participants with recorded sex (N = 20,003). The earlier
# Table 1 / Table S1 used the unrestricted 20,021; every model ran on 20,003.
#
# Inputs:  Build/derived/qc/user_analytic_frame_1X.rds
#          Build/derived/qc/bp_reading_analytic_frame_1X.rds
#          Build/derived/09_glp1_flag.rds
# Outputs (Analyses/0_cohort/derived/results/):
#   table_1_npj_20003.csv          — Table 1 cells in the manuscript's format
#   table_s1_npj_20003.csv         — Table S1 rows (user level + reading level)

suppressPackageStartupMessages(library(tidyverse))

project_dir <- "./PRIMARY"
results_dir <- file.path(project_dir, "Analyses", "0_cohort", "derived", "results")
HYP <- c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")

u <- read_rds(file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")) %>%
  left_join(read_rds(file.path(project_dir, "Build", "derived", "09_glp1_flag.rds")), by = "STUDY_ID") %>%
  filter(BASELINE_AHA %in% HYP, SEX %in% c("male", "female"))
stopifnot(nrow(u) == 20003)

# --- Table 1 ---------------------------------------------------------------
msd     <- function(x) sprintf("%.1f±%.1f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
npct    <- function(k, n) sprintf("%s [%.1f]", format(k, big.mark = ",", trim = TRUE), 100 * k / n)
med_iqr <- function(x) { q <- quantile(x, c(.25, .5, .75), na.rm = TRUE, names = FALSE)
                         sprintf("%.0f [%.0f–%.0f]", q[2], q[1], q[3]) }
rows_for <- function(d) {
  n <- nrow(d)
  c(N          = format(n, big.mark = ","),
    Age        = msd(d$AGE_AT_START),
    Sex        = sprintf("%s/%s", format(sum(d$SEX == "male"), big.mark = ","),
                         format(sum(d$SEX == "female"), big.mark = ",")),
    BMI        = msd(d$BASELINE_BMI),
    SBP        = msd(d$BASELINE_SBP),
    DBP        = msd(d$BASELINE_DBP),
    Period     = med_iqr(d$MO6_DAY),
    Tobacco    = npct(sum(d$TOBACCO_NICOTINE_EVER_USED), n),
    HTN_est    = npct(sum(d$HTN_MED_BASELINE == "stable_on_meds"), n),
    HTN_rec    = npct(sum(d$HTN_MED_BASELINE == "recent_initiator"), n),
    HTN_mid    = npct(sum(d$HTN_MED_BASELINE == "mid_study_initiator"), n),
    GLP_prior  = npct(sum(d$GLP1_LOG_COUNT_BASELINE > 0), n),
    GLP_during = npct(sum(d$GLP1_LOG_COUNT_BASELINE == 0 & d$GLP1_LOG_COUNT_STUDY > 0), n))
}
t1 <- bind_cols(row = names(rows_for(u)), Overall = rows_for(u),
                map_dfc(set_names(HYP, c("Elevated", "Stage 1", "Stage 2")),
                        ~ rows_for(filter(u, BASELINE_AHA == .x))))
write_csv(t1, file.path(results_dir, "table_1_npj_20003.csv"))
cat("--- Table 1 (N = 20,003) ---\n"); print(as.data.frame(t1), row.names = FALSE)

# --- Table S1: completeness ------------------------------------------------
comp <- function(d, cols, denom) {
  k <- sum(stats::complete.cases(d[, cols, drop = FALSE]))
  c(n = format(k, big.mark = ","), pct = sprintf("%.1f", 100 * k / denom))
}
wear <- c(RHR = "RHR", HRV = "HRV", `Physical Activity` = "STRAIN",
          `Sleep Duration` = "SLEEP_DURATION", `Sleep Consistency` = "SLEEP_CONSISTENCY")
win  <- c(Baseline = "_BASELINE_MEAN", `3-Month` = "_MO3_MEAN", `6-Month` = "_MO6_MEAN")
user_spec <- c(
  list(`3-Month Blood Pressure` = c("MO3_SBP", "MO3_DBP"),
       `6-Month Blood Pressure` = c("MO6_SBP", "MO6_DBP"),
       `12-Month Blood Pressure` = c("MO12_SBP", "MO12_DBP"),
       Age = "AGE_AT_START", Sex = "SEX", `Baseline BMI` = "BASELINE_BMI"),
  set_names(as.list(c(outer(wear, win, paste0))),
            c(outer(names(wear), names(win), function(m, w) paste(w, m)))),
  list(`Days Journaled` = "PCT_DAYS_JOURNALED",
       `Average Daily App Time` = "LOG_AVG_DAILY_TIME_IN_APP_MIN",
       `Days Engaged with Notifications` = "PCT_DAYS_ENGAGED_NOTIFICATIONS",
       `Days Used AI Coach` = "PCT_DAYS_USED_COACH",
       `Days Fulfilled Sleep Need` = "AVG_SLEEP_NEED_PCT",
       `Days Met Physical Activity Recommendation` = "PCT_DAYS_MET_OPTIMAL_STRAIN",
       `Number of Blood Pressure Readings` = "LOG_N_BP_READINGS_STUDY"))
s1_user <- imap_dfr(user_spec, ~ tibble(level = "user", variable = .y, !!!as.list(comp(u, .x, nrow(u)))))

bp <- read_rds(file.path(project_dir, "Build", "derived", "qc", "bp_reading_analytic_frame_1X.rds")) %>%
  filter(BASELINE_AHA %in% HYP, STUDY_ID %in% u$STUDY_ID)
read_spec <- set_names(paste0(wear, "_LOOKBACK_MEAN"), names(wear))
s1_read <- imap_dfr(read_spec, ~ tibble(level = "reading", variable = .y, !!!as.list(comp(bp, .x, nrow(bp)))))
s1 <- bind_rows(tibble(level = "user", variable = "N", n = format(nrow(u), big.mark = ","), pct = ""),
                s1_user,
                tibble(level = "reading", variable = "N", n = format(nrow(bp), big.mark = ","), pct = ""),
                s1_read)
write_csv(s1, file.path(results_dir, "table_s1_npj_20003.csv"))
cat("\n--- Table S1 (N = 20,003 users;", format(nrow(bp), big.mark = ","), "readings) ---\n")
print(as.data.frame(s1), row.names = FALSE)
