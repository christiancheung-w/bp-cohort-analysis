# 08_table_s2_npj_add_cohorts.R
# Table S2 rows for the two sensitivity cohorts reported in Tables S3, S5, S6
# and S7 but absent from the modeled-n table: no recorded GLP-1 receptor
# agonist use, and excluding antihypertensive initiation during follow-up.
# Those scripts did not emit *_sizes_std.csv, so n's are taken from their
# result tables where available and otherwise recounted with the same
# complete-case rules, each recount gated against a model n already reported.
#
# Inputs:  Analyses/0_cohort/derived/results/table_s2_npj_paper.csv
#          Analyses/6_2/derived/results/6_2_table_sensitivity.csv
#          Analyses/7_2/derived/results/7_2_table_sensitivity.csv
#          Analyses/8_engagement/derived/results/8_table_sensitivity.csv
#          Analyses/mid_excluded_sensitivity/derived/results/5_1_mid_excluded_n_by_stratum.csv
#          Build/derived/qc/user_analytic_frame_1X.rds, Build/derived/09_glp1_flag.rds
#          Fit / diagnostics column: 5_1, 5_5 and sensitivity r2.csv / moments; 7_2_model_summaries.csv;
#          6_2 and 8_engagement calibration and VIF CSVs (primary-cohort fits; others shown as —)
# Output:  Analyses/0_cohort/derived/results/table_s2_npj_paper_full.csv

suppressPackageStartupMessages(library(tidyverse))

project_dir <- "./PRIMARY"
A <- file.path(project_dir, "Analyses")
results_dir <- file.path(A, "0_cohort", "derived", "results")
HYP <- c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")
fmt <- function(x) format(x, big.mark = ",", trim = TRUE)

u <- read_rds(file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")) %>%
  left_join(read_rds(file.path(project_dir, "Build", "derived", "09_glp1_flag.rds")) %>%
              select(STUDY_ID, GLP1_EVER), by = "STUDY_ID") %>%
  filter(BASELINE_AHA %in% HYP, SEX %in% c("male", "female"))
stopifnot(nrow(u) == 20003)
cohorts <- list(primary      = u,
                glp1_naive   = u %>% filter(!GLP1_EVER),
                mid_excluded = u %>% filter(HTN_MED_BASELINE != "mid_study_initiator"))
LABEL <- c(glp1_naive = "No recorded GLP-1 receptor agonist use",
           mid_excluded = "Excluding antihypertensive initiation during follow-up")
ceiling_n <- map_int(cohorts, nrow)
cat("Cohort ceilings:\n"); print(ceiling_n)

# --- 6-month trajectory LMM: recount users and rows with the models' complete-case rule ---
traj_n <- function(d, adjusted) {
  times <- c("BASELINE", "MO3", "MO6")
  cols  <- c(outer(times, c("SBP", "TIME_OF_DAY", "SEASON"), paste, sep = "_"))
  long <- d %>%
    select(STUDY_ID, AGE_AT_START, SEX, BASELINE_BMI, all_of(cols)) %>%
    pivot_longer(all_of(cols), names_to = c("Time", ".value"),
                 names_pattern = "(BASELINE|MO3|MO6)_(.+)") %>%
    filter(!is.na(SBP))
  if (adjusted) long <- long %>% filter(!is.na(AGE_AT_START), !is.na(SEX), !is.na(BASELINE_BMI),
                                        !is.na(TIME_OF_DAY), !is.na(SEASON))
  c(n_users = n_distinct(long$STUDY_ID), n_obs = nrow(long))
}
chk_primary <- c(traj_n(u, TRUE), traj_n(u, FALSE))
stopifnot(all(chk_primary == c(19986, 50209, 20003, 50252)))   # Table S2 primary rows
ref_mid <- read_csv(file.path(A, "mid_excluded_sensitivity", "derived", "results",
                              "5_1_mid_excluded_n_by_stratum.csv"), show_col_types = FALSE) %>%
  filter(cohort == "mid_excluded", model %in% c("sbp_adj", "sbp_unadj")) %>%
  group_by(model) %>% summarise(n = sum(n_users), .groups = "drop")
stopifnot(traj_n(cohorts$mid_excluded, TRUE)[["n_users"]]  == ref_mid$n[ref_mid$model == "sbp_adj"],
          traj_n(cohorts$mid_excluded, FALSE)[["n_users"]] == ref_mid$n[ref_mid$model == "sbp_unadj"])
cat("Trajectory recount gates (primary, mid-excluded): OK\n")
traj_rows <- map_dfr(names(LABEL), function(k) {
  map_dfr(c(Adjusted = TRUE, Unadjusted = FALSE), function(adj) {
    v <- traj_n(cohorts[[k]], adj)
    tibble(Analysis = "6-month BP trajectory (LMM)", Cohort = LABEL[[k]], Outcome = "-",
           `Predictor set` = names(which(c(Adjusted = TRUE, Unadjusted = FALSE) == adj)),
           n_users = v[["n_users"]], n_obs = v[["n_obs"]], events = "-", ceiling = ceiling_n[[k]])
  })
})

# --- Wearable logistic: n and events from Table S6's source; ceiling = cohort with >= 8 baseline wear days ---
w <- read_csv(file.path(A, "6_2", "derived", "results", "6_2_table_sensitivity.csv"), show_col_types = FALSE) %>%
  distinct(cohort, outcome, n, n_events)
wear_ceiling <- map_int(cohorts, ~ sum(.x$BASELINE_N_DAYS >= 8, na.rm = TRUE))
stopifnot(wear_ceiling[["primary"]] == 14230)                     # Table S2 primary wearable ceiling
wear_rows <- map_dfr(names(LABEL), function(k) {
  ww <- w %>% filter(cohort == k); stopifnot(n_distinct(ww$n) == 1)
  tibble(Analysis = "Wearable metrics (logistic)", Cohort = LABEL[[k]], Outcome = "≥5 mmHg reduction",
         `Predictor set` = "Primary", n_users = ww$n[1], n_obs = ww$n[1],
         events = sprintf("%s / %s", fmt(ww$n_events[ww$outcome == "SBP"]), fmt(ww$n_events[ww$outcome == "DBP"])),
         ceiling = wear_ceiling[[k]])
})

# --- Reading-level LMM: from Table S7's source; ceiling = cohort n ---
r <- read_csv(file.path(A, "7_2", "derived", "results", "7_2_table_sensitivity.csv"), show_col_types = FALSE) %>%
  distinct(cohort, n_users, n_obs)
stopifnot(r$n_users[r$cohort == "primary"] == 16161)
read_rows <- map_dfr(names(LABEL), function(k) {
  rr <- r %>% filter(cohort == k); stopifnot(nrow(rr) == 1)
  tibble(Analysis = "Reading-level within/between-person model (LMM)", Cohort = LABEL[[k]], Outcome = "-",
         `Predictor set` = "Anchored", n_users = rr$n_users, n_obs = rr$n_obs, events = "-", ceiling = ceiling_n[[k]])
})

# --- Engagement logistic: n from Table S5's source; events recounted with the model's complete-case rule ---
PRED <- c(AVG_SLEEP_NEED_PCT = "Sleep need achieved", LOG_AVG_DAILY_TIME_IN_APP_MIN = "Daily app time",
          LOG_N_BP_READINGS_STUDY = "Blood pressure readings", PCT_DAYS_ENGAGED_NOTIFICATIONS = "Notification engagement",
          PCT_DAYS_JOURNALED = "Journal completion", PCT_DAYS_MET_OPTIMAL_STRAIN = "Physical activity recommendation",
          PCT_DAYS_USED_COACH = "AI Coach engagement")
e <- read_csv(file.path(A, "8_engagement", "derived", "results", "8_table_sensitivity.csv"), show_col_types = FALSE) %>%
  mutate(predictor = sub("_per10$", "", predictor)) %>% distinct(cohort, predictor, n)
eng_count <- function(d, pred) {
  dd <- d %>% filter(!is.na(.data[[pred]]), !is.na(AGE_AT_START), !is.na(SEX), !is.na(BASELINE_BMI),
                     !is.na(BASELINE_SBP), !is.na(BASELINE_DBP), !is.na(MO6_SBP), !is.na(MO6_DBP))
  c(n = nrow(dd), sbp = sum(dd$MO6_SBP - dd$BASELINE_SBP <= -5), dbp = sum(dd$MO6_DBP - dd$BASELINE_DBP <= -5))
}
chk <- eng_count(u, "PCT_DAYS_JOURNALED"); stopifnot(all(chk == c(19986, 10029, 8284)))   # Table S2 primary row
eng_rows <- map_dfr(names(LABEL), function(k) map_dfr(names(PRED), function(p) {
  v <- eng_count(cohorts[[k]], p)
  n_ref <- e$n[e$cohort == k & e$predictor == p]; stopifnot(length(n_ref) == 1, v[["n"]] == n_ref)
  tibble(Analysis = "Engagement (logistic, separate models)", Cohort = LABEL[[k]], Outcome = "≥5 mmHg reduction",
         `Predictor set` = PRED[[p]], n_users = v[["n"]], n_obs = v[["n"]],
         events = sprintf("%s / %s", fmt(v[["sbp"]]), fmt(v[["dbp"]])), ceiling = ceiling_n[[k]])
}))
cat("Engagement n gates (all predictors, both cohorts): OK\n")

new <- bind_rows(traj_rows, wear_rows, read_rows, eng_rows) %>%
  transmute(section = "Sensitivity analyses", Analysis, Cohort, Outcome, `Predictor set`,
            `Participants (n)` = fmt(n_users), `Observations (n)` = fmt(n_obs), `Events (SBP / DBP)` = events,
            `Cohort ceiling (n)` = fmt(ceiling), `% of ceiling` = sprintf("%.1f", 100 * n_users / ceiling))

# Insert after the matching analysis block of the existing sensitivity rows (same column set)
base <- read_csv(file.path(results_dir, "table_s2_npj_paper.csv"), show_col_types = FALSE,
                 col_types = cols(.default = col_character()))
stopifnot(identical(names(base), names(new)))
out <- base[0, ]
for (i in seq_len(nrow(base))) {
  out <- bind_rows(out, base[i, ])
  nxt <- if (i < nrow(base)) base$Analysis[i + 1] else ""
  if (base$section[i] == "Sensitivity analyses" && nxt != base$Analysis[i]) {
    out <- bind_rows(out, new %>% filter(Analysis == base$Analysis[i]))
    new <- new %>% filter(Analysis != base$Analysis[i])
  }
}
stopifnot(nrow(new) == 0)
# --- Fit / diagnostics column -------------------------------------------------
r2_cell <- function(path, mod) {
  r <- read_csv(path, show_col_types = FALSE)
  v <- function(o) r[r$model == paste0(o, "_", mod), ]
  sprintf("R² %.2f / %.2f (SBP), %.2f / %.2f (DBP)",
          v("sbp")$marginal_r2, v("sbp")$conditional_r2, v("dbp")$marginal_r2, v("dbp")$conditional_r2)
}
lmm_r2 <- list(
  "6-month BP trajectory (LMM) Primary cohort"                        = file.path(A, "5_1", "derived", "results", "5_1_r2.csv"),
  "12-month BP trajectory (LMM, 4 timepoints) Primary cohort"         = file.path(A, "5_5", "derived", "results", "5_5_r2.csv"),
  "6-month BP trajectory (LMM) Full cohort (including normal BP)"      = file.path(A, "5_1", "Sensitivity", "full_cohort", "derived", "results", "r2.csv"),
  "6-month BP trajectory (LMM) No recorded antihypertensive medication" = file.path(A, "5_1", "Sensitivity", "htn_naive", "derived", "results", "r2.csv"),
  "6-month BP trajectory (LMM) Twice-weekly readings"                  = file.path(A, "5_1", "Sensitivity", "2xWeekly", "derived", "results", "r2.csv"))
m72 <- read_csv(file.path(A, "7_2", "derived", "results", "7_2_model_summaries.csv"), show_col_types = FALSE)
rl  <- sprintf("R² %.3f / %.2f (SBP), %.3f / %.2f (DBP)",
               m72$marginal_r2[m72$model == "SBP-baseline"], m72$conditional_r2[m72$model == "SBP-baseline"],
               m72$marginal_r2[m72$model == "DBP-baseline"], m72$conditional_r2[m72$model == "DBP-baseline"])
cal_max <- function(d) 100 * max(abs(d$mean_pred - d$mean_obs))
c62 <- read_csv(file.path(A, "6_2", "derived", "results", "6_2_calibration_data.csv"), show_col_types = FALSE)
v62 <- read_csv(file.path(A, "6_2", "derived", "results", "6_2_vif_primary.csv"), show_col_types = FALSE)
wear_cell <- sprintf("Calibration ≤%.1f / %.1f pp; VIF ≤%.1f",
                     cal_max(filter(c62, outcome == "SBP_responder")), cal_max(filter(c62, outcome == "DBP_responder")), max(v62$VIF))
c8  <- read_csv(file.path(A, "8_engagement", "derived", "results", "8_calibration_data.csv"), show_col_types = FALSE) %>% filter(frame == "sep")
v8  <- read_csv(file.path(A, "8_engagement", "derived", "results", "8_vif_joint.csv"), show_col_types = FALSE)
PREDKEY <- c("Sleep need achieved" = "AVG_SLEEP_NEED_PCT_per10", "Daily app time" = "LOG_AVG_DAILY_TIME_IN_APP_MIN",
             "Blood pressure readings" = "LOG_N_BP_READINGS_STUDY", "Notification engagement" = "PCT_DAYS_ENGAGED_NOTIFICATIONS_per10",
             "Journal completion" = "PCT_DAYS_JOURNALED_per10", "Physical activity recommendation" = "PCT_DAYS_MET_OPTIMAL_STRAIN_per10",
             "AI Coach engagement" = "PCT_DAYS_USED_COACH_per10")
eng_cell <- function(pred_label, cohort_label) {
  flavor <- c("Primary cohort" = "Primary", "No recorded antihypertensive medication" = "HTNNaive")[cohort_label]
  if (is.na(flavor)) return("—")
  k <- PREDKEY[[pred_label]]
  d <- c8 %>% filter(flavor == !!flavor, predictor == k)
  stopifnot(nrow(d) == 20)
  cell <- sprintf("Calibration ≤%.1f / %.1f pp", cal_max(filter(d, outcome == "SBP")), cal_max(filter(d, outcome == "DBP")))
  if (flavor == "Primary") cell <- sprintf("%s; VIF %.2f", cell, max(v8$VIF[v8$predictor == k]))
  cell
}
out <- out %>% rowwise() %>% mutate(`Fit / diagnostics` = {
  if (!is.null(lmm_r2[[paste(Analysis, Cohort)]])) r2_cell(lmm_r2[[paste(Analysis, Cohort)]], if (`Predictor set` == "Adjusted") "adj" else "unadj")
  else if (Analysis == "Reading-level within/between-person model (LMM)" && Cohort == "Primary cohort") rl
  else if (Analysis == "Wearable metrics (logistic)" && Cohort == "Primary cohort") wear_cell
  else if (Analysis == "Engagement (logistic, separate models)") eng_cell(`Predictor set`, Cohort)
  else "—" }) %>% ungroup()
write_csv(out, file.path(results_dir, "table_s2_npj_paper_full.csv"))
cat("\nRows:", nrow(base), "->", nrow(out), "\n")
print(as.data.frame(out %>% filter(Cohort %in% LABEL)), row.names = FALSE)
cat("\nWrote:", file.path(results_dir, "table_s2_npj_paper_full.csv"), "\n")
