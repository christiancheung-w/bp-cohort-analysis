# run_all.R
# Runs the pipeline in dependency order, each script in a fresh Rscript
# process with the project root as working directory.
#
# Usage (from the project root):
#   Rscript run_all.R                 # full run
#   Rscript run_all.R <script path>   # resume from that script
#
# Logs:    derived/run_logs/<script>.log
# Summary: derived/run_logs/run_all_summary.tsv
# Stops at the first script that exits non-zero.

scripts <- c(
  # Build
  "PRIMARY/Build/01_bp_timepoints.R",
  "PRIMARY/Build/02_bp_outcomes.R",
  "PRIMARY/Build/03_demographics_clinical.R",
  "PRIMARY/Build/04_lab_weight_outcomes.R",
  "PRIMARY/Build/05_whoop_window_means.R",
  "PRIMARY/Build/06_le8_tobacco_nicotine.R",
  "PRIMARY/Build/07_engagement.R",
  "PRIMARY/Build/08_lmm_covariates.R",
  "PRIMARY/Build/09_glp1_flag.R",
  "PRIMARY/Build/98_qc_user_frame.R",
  "PRIMARY/Build/99_assemble_user_frame.R",
  "PRIMARY/Build/bp_reading_long_frame.R",
  "PRIMARY/Build/Sensitivity/2xWeekly/build_user_frame.R",

  # Cohort
  "PRIMARY/Analyses/0_cohort/00_check_cohort_defs.R",
  "PRIMARY/Analyses/0_cohort/01_table1.R",
  "PRIMARY/Analyses/0_cohort/02_completeness_supp.R",
  "PRIMARY/Analyses/0_cohort/05_table1_npj_rows.R",

  # 5.1 six-month trajectories
  "PRIMARY/Analyses/5_1/5_1_primary_lmm.R",
  "PRIMARY/Analyses/5_1/5_1_diagnostics.R",
  "PRIMARY/Analyses/5_1/5_1_figure_1a.R",
  "PRIMARY/Analyses/5_1/5_1_table_2.R",
  "PRIMARY/Analyses/5_1/Sensitivity/2xWeekly/5_1_sensitivity_2xWeekly.R",
  "PRIMARY/Analyses/5_1/Sensitivity/age_sex/5_1_sensitivity_age_sex.R",
  "PRIMARY/Analyses/5_1/Sensitivity/age_sex/5_1_recover_age_sex_post.R",
  "PRIMARY/Analyses/5_1/Sensitivity/full_cohort/5_1_sensitivity_full_cohort.R",
  "PRIMARY/Analyses/5_1/Sensitivity/htn_naive/5_1_sensitivity_htn_naive.R",
  "PRIMARY/Analyses/5_1/Sensitivity/rtm_barnett/5_1_sensitivity_rtm_barnett.R",
  "PRIMARY/Analyses/5_1/5_1_comparison_1X_vs_2xWeekly.R",

  # 5.2-5.4 six-month descriptives
  "PRIMARY/Analyses/5_2/5_2_clinical_bins.R",
  "PRIMARY/Analyses/5_3/5_3_distributional_comparison.R",
  "PRIMARY/Analyses/5_4/5_4_aha_category_shift.R",

  # 5.5 twelve-month trajectories
  "PRIMARY/Analyses/5_5/5_5_primary_lmm.R",
  "PRIMARY/Analyses/5_5/5_5_diagnostics.R",
  "PRIMARY/Analyses/5_5/5_5_figure_2.R",
  "PRIMARY/Analyses/5_5/5_5_table_s3.R",
  "PRIMARY/Analyses/5_5/5_5_aha_distribution_12mo.R",
  "PRIMARY/Analyses/5_5/Sensitivity/rtm_barnett/5_5_sensitivity_rtm_barnett.R",
  "PRIMARY/Analyses/5_1/Sensitivity/rtm_barnett_bootstrap/rtm_barnett_bootstrap.R",
  "PRIMARY/Analyses/5_1/5_1_table_2_npj.R",
  "PRIMARY/Analyses/5_5/5_5_table_s3_npj.R",

  # 6.2 wearable-signal logistic
  "PRIMARY/Analyses/6_2/6_2_logistic.R",

  # 7.2 reading-level within/between model
  "PRIMARY/Analyses/7_2/7_2_wb_lmm.R",
  "PRIMARY/Analyses/7_2/7_2_diagnostics.R",
  "PRIMARY/Analyses/7_2/7_2_sensitivities.R",
  "PRIMARY/Analyses/7_2/7_2_figure_4.R",

  # 8 engagement
  "PRIMARY/Analyses/8_engagement/8_engagement_logistic.R",
  "PRIMARY/Analyses/8_engagement/8_figure_forest.R",

  # Cross-chapter (read results from the chapters above)
  "PRIMARY/Analyses/glp1_naive_sensitivity/glp1_naive_sensitivity.R",
  "PRIMARY/Analyses/mid_excluded_sensitivity/mid_excluded_sensitivity.R",
  "PRIMARY/Analyses/5_1/5_1_table_sensitivity.R",
  "PRIMARY/Analyses/8_engagement/8_table_sensitivity.R",
  "PRIMARY/Analyses/6_2/6_2_table_sensitivity.R",
  "PRIMARY/Analyses/7_2/7_2_table_sensitivity.R",
  "PRIMARY/Analyses/0_cohort/03_modeled_n_supp.R",
  "PRIMARY/Analyses/0_cohort/04_modeled_n_word.R",
  "PRIMARY/Analyses/0_cohort/06_table_s2_npj.R",
  "PRIMARY/Analyses/0_cohort/07_table1_s1_npj_20003.R",
  "PRIMARY/Analyses/0_cohort/08_table_s2_npj_add_cohorts.R"
)

missing <- scripts[!file.exists(scripts)]
if (length(missing)) stop("Scripts not found:\n  ", paste(missing, collapse = "\n  "))

args  <- commandArgs(trailingOnly = TRUE)
start <- if (length(args)) match(args[1], scripts) else 1L
if (is.na(start)) stop("Unknown start script: ", args[1])

log_dir <- file.path("derived", "run_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
summary_path <- file.path(log_dir, "run_all_summary.tsv")
if (start == 1L || !file.exists(summary_path)) {
  writeLines(paste("script", "exit_status", "seconds", "finished", sep = "\t"), summary_path)
}

cat("R:", R.version.string, "\n")
cat("Running", length(scripts) - start + 1L, "of", length(scripts), "scripts\n\n")

for (s in scripts[start:length(scripts)]) {
  log_path <- file.path(log_dir, paste0(gsub("/", "__", sub("\\.R$", "", s)), ".log"))
  t0 <- Sys.time()
  cat(format(t0, "%H:%M:%S"), " ", s, " ... ", sep = "")
  status <- system2("Rscript", shQuote(s), stdout = log_path, stderr = log_path)
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  cat(if (status == 0) "ok" else paste("FAILED (exit", status, ")"), " [", secs, "s]\n", sep = "")
  cat(paste(s, status, secs, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sep = "\t"),
      "\n", file = summary_path, append = TRUE, sep = "")
  if (status != 0) stop("Stopped at ", s, ". See ", log_path)
}

cat("\nAll scripts finished.\n")
