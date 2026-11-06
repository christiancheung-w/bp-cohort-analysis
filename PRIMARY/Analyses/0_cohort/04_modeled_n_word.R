# 04_modeled_n_word.R
# -----------------------------------------------------------------------------
# Word-pasteable rendering of Table A (modeled n) + Table A.2 (sensitivity).
# Reads the folded CSVs from 03_modeled_n_supp.R and emits a tab-separated file:
# open in Excel (or paste + Convert Text to Table in Word) to get a clean grid.
#
# Transformations (presentation only — no numbers change):
#   - drop `section`
#   - consistent cohort labels (Antihypertensive Medication Naive; Full cohort (including normal BP))
#   - presentable outcome / predictor-set titles (engagement labels reuse the
#     canonical pred_label_map from 8_engagement_logistic.R)
#   - merge n_events_SBP / n_events_DBP into one "Events (SBP / DBP)" column
#   - every blank / NA cell rendered as "-"
#
# Run with working directory = project root.
# Output: Analyses/0_cohort/derived/results/table_A_A2_word_pasteable.tsv
# -----------------------------------------------------------------------------

library(tidyverse)

results_dir <- "./PRIMARY/Analyses/0_cohort/derived/results"

# --- Label maps ---------------------------------------------------------------
cohort_map <- c(
  "hypertensive"                  = "Hypertensive",
  "hypertensive (MO12 completers)" = "Hypertensive (12-mo completers)",
  "hypertensive MO12 completers"   = "Hypertensive (12-mo completers)",
  "full"                          = "Full cohort (including normal BP)",
  "FullCohort"                    = "Full cohort (including normal BP)",
  "htn_naive"                     = "Antihypertensive Medication Naive",
  "HTN_naive"                     = "Antihypertensive Medication Naive",
  "2xWeekly"                      = "2×Weekly"
)  # age/sex strata (Age 18-39, Female, ...) pass through unchanged

outcome_map <- c(
  "responder"        = "BP responder",
  "HRV"              = "Heart rate variability (HRV)",
  "RHR"              = "Resting heart rate (RHR)",
  "SleepConsistency" = "Sleep consistency",
  "SleepDuration"    = "Sleep duration",
  "Strain"           = "Strain"
)

pred_map <- c(
  "adjusted"              = "Adjusted",
  "unadjusted"            = "Unadjusted",
  "Primary"               = "Primary",
  "anchored"              = "Anchored",
  "PlusDHDL"              = "Primary + Δ HDL",
  "PlusDHbA1c"            = "Primary + Δ HbA1c",
  "PlusDNonHDL"           = "Primary + Δ Non-HDL",
  "PlusDTG"               = "Primary + Δ Triglycerides",
  "PlusDweight"           = "Primary + Δ Weight",
  "RTM-adjusted Δ"        = "RTM-adjusted change",
  "RTM-adjusted Δ (12mo)" = "RTM-adjusted change (12-mo)",
  "median paired Δ"       = "Median paired change",
  "joint"                 = "Joint (all 7 predictors)",
  # engagement predictors — clean identity (scaling/transform annotations dropped;
  # this is a sample-size table, not an effect-size table)
  "PCT_DAYS_JOURNALED_per10"             = "Days journaled",
  "LOG_AVG_DAILY_TIME_IN_APP_MIN"        = "Daily app time",
  "PCT_DAYS_ENGAGED_NOTIFICATIONS_per10" = "Days engaged with notifications",
  "PCT_DAYS_USED_COACH_per10"            = "Days used AI Coach",
  "AVG_SLEEP_NEED_PCT_per10"             = "Sleep need achieved",
  "PCT_DAYS_MET_OPTIMAL_STRAIN_per10"    = "Days met optimal strain",
  "LOG_N_BP_READINGS_STUDY"              = "BP readings logged"
)

# analysis-label cleanup so cohort terminology matches the Cohort column
clean_analysis <- function(x) {
  x <- gsub("full-cohort, Normal ref", "Full cohort (including normal BP)", x, fixed = TRUE)
  x <- gsub("full-cohort",             "Full cohort (including normal BP)", x, fixed = TRUE)
  x <- gsub("FullCohort",              "Full cohort (including normal BP)", x, fixed = TRUE)
  x <- gsub("HTN-naive subset", "Antihypertensive Medication Naive", x, fixed = TRUE)
  x <- gsub("HTN-naive",        "Antihypertensive Medication Naive", x, fixed = TRUE)
  x <- gsub("HTN_naive",        "Antihypertensive Medication Naive", x, fixed = TRUE)
  x
}

remap <- function(x, map) ifelse(x %in% names(map), map[x], x)
dash  <- function(x) { x <- as.character(x); ifelse(is.na(x) | x == "", "-", x) }
fmtn  <- function(x) ifelse(is.na(x), "-", format(x, big.mark = ",", trim = TRUE))
fmt_events <- function(sbp, dbp) ifelse(is.na(sbp) & is.na(dbp), "-",
                                        paste0(fmtn(sbp), " / ", fmtn(dbp)))

HDR <- c("Analysis", "Cohort", "Outcome", "Predictor set", "Participants (n)",
         "Observations (n)", "Events (SBP / DBP)", "Cohort ceiling (n)", "% of ceiling")

presentable <- function(f) {
  df <- read_csv(file.path(results_dir, f), show_col_types = FALSE)
  tibble(
    Analysis            = clean_analysis(df$analysis),
    Cohort              = dash(remap(df$cohort, cohort_map)),
    Outcome             = dash(remap(df$outcome, outcome_map)),
    `Predictor set`     = dash(remap(df$predictor_set, pred_map)),
    `Participants (n)`  = fmtn(df$n_users),
    `Observations (n)`  = fmtn(df$n_obs),
    `Events (SBP / DBP)` = fmt_events(df$n_events_SBP, df$n_events_DBP),
    `Cohort ceiling (n)` = fmtn(df$ceiling_n),
    `% of ceiling`      = ifelse(is.na(df$`% of ceiling`), "-",
                                  sprintf("%.1f", as.numeric(df$`% of ceiling`)))
  )
}

primary <- presentable("table_A_modeled_n.csv")
sens    <- presentable("table_A2_sensitivity_n.csv")

emit <- function(df) c(paste(HDR, collapse = "\t"),
                       apply(df, 1, function(r) paste(r, collapse = "\t")))

out <- c(
  "Table A. Modeled N — Primary analyses",
  emit(primary),
  "",
  "Table A.2. Modeled N — Sensitivity analyses",
  emit(sens)
)

out_path <- file.path(results_dir, "table_A_A2_word_pasteable.tsv")
writeLines(out, out_path)
cat("Wrote", out_path, "\n")
cat("Primary rows:", nrow(primary), " Sensitivity rows:", nrow(sens), "\n")
