# 02_completeness_supp.R
# -----------------------------------------------------------------------------
# Supplementary "Table B": analytic-variable completeness (missingness) for the
# §2.2 primary (hypertensive) cohort. Reports, for every modeled outcome /
# covariate / predictor, the n and % of users non-missing — i.e. the CEILING n
# any model using that variable could run on. The per-analysis n a specific
# model actually ran on (always ≤ this ceiling, due to joint complete-case
# deletion) is the separate "Table A" supplement.
#
# Two tables, one markdown:
#   B    — user-level completeness (primary cohort, N from user frame)
#   B.2  — reading-level completeness for the §7.2 predictors (different grain:
#          one row per BP reading), since those live in the reading-level frame.
#
# Sources cohort_defs.R. Run with working directory = project root.
#
# Outputs (Analyses/0_cohort/derived/results/):
#   02_completeness_user_level.csv
#   02_completeness_reading_level.csv
#   02_completeness_paper.md
#   02_completeness_supp.log
# -----------------------------------------------------------------------------

library(tidyverse)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")

project_dir <- "./PRIMARY"
user_path   <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
read_path   <- file.path(project_dir, "Build", "derived", "qc", "bp_reading_analytic_frame_1X.rds")
results_dir <- file.path(project_dir, "Analyses", "0_cohort", "derived", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "02_completeness_supp.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)
cat("02_completeness_supp —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Completeness helper -----------------------------------------------------
# A row is "complete" for a user/reading iff non-missing on ALL listed columns
# (so BP timepoints require both SBP and DBP). Errors loudly if a column is
# absent — we never want to silently report a missing variable as 0% complete.
completeness <- function(d, spec, denom = nrow(d)) {
  map_dfr(spec, function(r) {
    miss_cols <- setdiff(r$cols, names(d))
    if (length(miss_cols) > 0) stop("Column(s) absent from frame: ", paste(miss_cols, collapse = ", "))
    n_ok <- sum(stats::complete.cases(d[, r$cols, drop = FALSE]))
    tibble(Section = r$section, Variable = r$label,
           `n complete` = format(n_ok, big.mark = ","),
           `% complete` = sprintf("%.1f", 100 * n_ok / denom))
  })
}

# --- WHOOP trajectory spec (§7.1): 5 co-measured metrics x 3 windows ---------
whoop_metrics <- c(RHR = "RHR", HRV = "HRV", Strain = "STRAIN",
                   `Sleep duration` = "SLEEP_DURATION",
                   `Sleep consistency` = "SLEEP_CONSISTENCY")
whoop_windows <- c(baseline = "_BASELINE_MEAN", `3-mo` = "_MO3_MEAN", `6-mo` = "_MO6_MEAN")
whoop_spec <- unlist(lapply(names(whoop_metrics), function(m) {
  lapply(names(whoop_windows), function(w) {
    list(section = "Predictor — WHOOP daily-metric trajectory (§7.1)",
         label   = paste0(m, " — ", w, " mean"),
         cols    = paste0(whoop_metrics[[m]], whoop_windows[[w]]))
  })
}), recursive = FALSE)

# --- Full user-level variable spec -------------------------------------------
user_spec <- c(
  list(
    list(section = "Outcome", label = "BP @ 3 months (SBP & DBP)",            cols = c("MO3_SBP", "MO3_DBP")),
    list(section = "Outcome", label = "BP @ 6 months (SBP & DBP) ‖",          cols = c("MO6_SBP", "MO6_DBP")),
    list(section = "Outcome", label = "BP @ 12 months (SBP & DBP)",           cols = c("MO12_SBP", "MO12_DBP")),
    list(section = "Outcome", label = "Δ weight (kg)",                        cols = "DELTA_WEIGHT_KG"),
    list(section = "Outcome", label = "Δ HDL cholesterol",                    cols = "DELTA_HDL_CHOLESTEROL"),
    list(section = "Outcome", label = "Δ non-HDL cholesterol",                cols = "DELTA_NON_HDL_CHOLESTEROL"),
    list(section = "Outcome", label = "Δ triglycerides",                      cols = "DELTA_TRIGLYCERIDES"),
    list(section = "Outcome", label = "Δ HbA1c",                              cols = "DELTA_HBA1C"),
    list(section = "Covariate (shared)", label = "Age at baseline",           cols = "AGE_AT_START"),
    list(section = "Covariate (shared)", label = "Sex",                       cols = "SEX"),
    list(section = "Covariate (shared)", label = "Baseline BMI",              cols = "BASELINE_BMI")
  ),
  whoop_spec,
  list(
    list(section = "Predictor — behavioral study-window mean (§6.1)", label = "Sleep duration (study mean)",       cols = "SLEEP_DURATION_STUDY_MEAN"),
    list(section = "Predictor — behavioral study-window mean (§6.1)", label = "Strain (study mean)",               cols = "STRAIN_STUDY_MEAN"),
    list(section = "Predictor — behavioral study-window mean (§6.1)", label = "Tobacco/nicotine use days (log)",   cols = "LOG_TOBACCO_NICOTINE_USE_DAYS"),
    list(section = "Predictor — behavioral study-window mean (§6.1)", label = "Tobacco/nicotine ever logged",      cols = "TOBACCO_NICOTINE_EVER_LOGGED"),
    list(section = "Predictor — engagement (§8)", label = "% days journaled",                cols = "PCT_DAYS_JOURNALED"),
    list(section = "Predictor — engagement (§8)", label = "Avg daily time in app (log)",      cols = "LOG_AVG_DAILY_TIME_IN_APP_MIN"),
    list(section = "Predictor — engagement (§8)", label = "% days engaged w/ notifications",  cols = "PCT_DAYS_ENGAGED_NOTIFICATIONS"),
    list(section = "Predictor — engagement (§8)", label = "% days used AI Coach",             cols = "PCT_DAYS_USED_COACH"),
    list(section = "Predictor — engagement (§8)", label = "Sleep % of target",                cols = "SLEEP_PCT_OF_TARGET"),
    list(section = "Predictor — engagement (§8)", label = "% days met optimal strain",        cols = "PCT_DAYS_MET_OPTIMAL_STRAIN"),
    list(section = "Predictor — engagement (§8)", label = "N BP readings, study (log)",        cols = "LOG_N_BP_READINGS_STUDY")
  )
)

# --- §7.2 reading-level spec -------------------------------------------------
reading_spec <- list(
  list(section = "Predictor — WHOOP lookback (§7.2)", label = "RHR (lookback mean)",                cols = "RHR_LOOKBACK_MEAN"),
  list(section = "Predictor — WHOOP lookback (§7.2)", label = "HRV (lookback mean)",                cols = "HRV_LOOKBACK_MEAN"),
  list(section = "Predictor — WHOOP lookback (§7.2)", label = "Strain (lookback mean)",             cols = "STRAIN_LOOKBACK_MEAN"),
  list(section = "Predictor — WHOOP lookback (§7.2)", label = "Sleep duration (lookback mean)",     cols = "SLEEP_DURATION_LOOKBACK_MEAN"),
  list(section = "Predictor — WHOOP lookback (§7.2)", label = "Sleep consistency (lookback mean)",  cols = "SLEEP_CONSISTENCY_LOOKBACK_MEAN")
)

# --- Compute -----------------------------------------------------------------
uf  <- read_rds(user_path)
hyp <- apply_cohort(uf, "hypertensive")
N_user <- nrow(hyp)
tab_user <- completeness(hyp, user_spec, denom = N_user)
cat(sprintf("User-level: hypertensive N = %s\n", format(N_user, big.mark = ",")))
print(tab_user, n = Inf, width = Inf)

# Reading-level: subset the reading frame to the hypertensive strata to match
# the primary-cohort framing. (§7.2's modeled n applies further filters —
# PERIOD == DURING, >=2 readings/user — which appear in Table A, not here.)
bp  <- read_rds(read_path) %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"))
N_read  <- nrow(bp)
N_rusers <- n_distinct(bp$STUDY_ID)
tab_read <- completeness(bp, reading_spec, denom = N_read)
cat(sprintf("\nReading-level: %s hypertensive readings from %s users\n",
            format(N_read, big.mark = ","), format(N_rusers, big.mark = ",")))
print(tab_read, n = Inf, width = Inf)

# --- Write CSVs --------------------------------------------------------------
write_csv(tab_user, file.path(results_dir, "02_completeness_user_level.csv"))
write_csv(tab_read, file.path(results_dir, "02_completeness_reading_level.csv"))

# --- Markdown ----------------------------------------------------------------
# Render a completeness table with section sub-headers (blank value cells).
render_completeness_md <- function(tab, denom_label) {
  lines <- c(
    paste0("| Variable | n complete | % of ", denom_label, " |"),
    "|:---|--:|--:|"
  )
  for (sec in unique(tab$Section)) {
    lines <- c(lines, paste0("| **", sec, "** |  |  |"))
    sub <- tab %>% filter(Section == sec)
    for (i in seq_len(nrow(sub))) {
      lines <- c(lines, paste0("| ", sub$Variable[i], " | ",
                               sub$`n complete`[i], " | ", sub$`% complete`[i], " |"))
    }
  }
  lines
}

md <- c(
  "# Table B. Analytic-variable completeness (§2.2 primary hypertensive cohort)",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  paste0("Completeness = users non-missing on the listed variable(s); it is the **maximum n** any model using that variable could run on. The per-analysis n actually used (≤ this, after joint complete-case deletion across a model's full variable set) is reported in Table A."),
  "",
  paste0("## B.1 User-level (N = ", format(N_user, big.mark = ","), " users)"),
  "",
  render_completeness_md(tab_user, "users"),
  "",
  paste0("## B.2 Reading-level — §7.2 predictors (", format(N_read, big.mark = ","),
         " readings from ", format(N_rusers, big.mark = ","), " users)"),
  "",
  render_completeness_md(tab_read, "readings"),
  "",
  "<sub>‖ BP responder outcomes (SBP / DBP at 6 mo; §6.1 / §6.2 / §7.1) are derived from 6-month BP and share its completeness.</sub>",
  "",
  "<sub>WHOOP daily metrics within a window are co-measured, so their completeness is identical up to a handful of users. §7.2 lookback predictors are summarized at the reading grain (B.2) because they vary per reading; all other variables are user-level (B.1).</sub>"
)
writeLines(md, file.path(results_dir, "02_completeness_paper.md"))
cat("\nWrote: 02_completeness_user_level.csv, 02_completeness_reading_level.csv, 02_completeness_paper.md\n")
cat("\n=== 02_completeness_supp done ===\n")
