# 5_1_table_sensitivity.R
# Supplementary table: adjusted 6-month change from baseline in SBP and DBP by
# baseline stratum, primary cohort vs each sensitivity cohort / subgroup.
# Reads existing result CSVs (no model refits); n = participants contributing to
# the adjusted SBP model, recounted with that model's complete-case rule.
#
# Inputs:
#   Analyses/5_1/derived/results/5_1_delta_contrasts.csv
#   Analyses/5_1/Sensitivity/full_cohort/derived/results/delta_contrasts.csv
#   Analyses/5_1/Sensitivity/htn_naive/derived/results/delta_contrasts.csv
#   Analyses/glp1_naive_sensitivity/derived/results/5_1_glp1_naive_delta_contrasts.csv
#   Analyses/mid_excluded_sensitivity/derived/results/5_1_mid_excluded_delta_contrasts.csv
#   Analyses/5_1/Sensitivity/age_sex/derived/results/delta_contrasts_by_subgroup.csv
#   Build/derived/qc/user_analytic_frame_1X.rds, Build/derived/09_glp1_flag.rds
# Outputs (Analyses/5_1/derived/results/):
#   5_1_table_sensitivity.csv        - numeric, one row per analysis x stratum
#   5_1_table_sensitivity_paper.csv  - formatted, 1 decimal, rows grouped by stratum

suppressPackageStartupMessages(library(tidyverse))
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")   # htn_baseline_on()

project_dir <- "./PRIMARY"
A <- file.path(project_dir, "Analyses")
results_dir <- file.path(A, "5_1", "derived", "results")

HYP   <- c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")
SHORT <- c("Elevated", "Stage 1", "Stage 2")

user <- read_rds(file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")) %>%
  left_join(read_rds(file.path(project_dir, "Build", "derived", "09_glp1_flag.rds")) %>%
              select(STUDY_ID, GLP1_EVER), by = "STUDY_ID") %>%
  filter(BASELINE_AHA %in% HYP, SEX %in% c("male", "female")) %>%
  mutate(age_band = case_when(AGE_AT_START >= 18 & AGE_AT_START <= 39 ~ "18–39",
                              AGE_AT_START >= 40 & AGE_AT_START <= 59 ~ "40–59",
                              AGE_AT_START >= 60 & AGE_AT_START <= 79 ~ "60–79"))

# --- n: participants with >= 1 complete row in the adjusted SBP model ---
adjusted_n <- function(u) {
  times <- c("BASELINE", "MO3", "MO6")
  cols  <- c(outer(times, c("SBP", "TIME_OF_DAY", "SEASON"), paste, sep = "_"))
  u %>%
    select(STUDY_ID, BASELINE_AHA, AGE_AT_START, SEX, BASELINE_BMI, HTN_MED_BASELINE, all_of(cols)) %>%
    pivot_longer(all_of(cols), names_to = c("Time", ".value"),
                 names_pattern = "(BASELINE|MO3|MO6)_(.+)") %>%
    filter(!is.na(SBP), !is.na(AGE_AT_START), !is.na(SEX), !is.na(BASELINE_BMI),
           !is.na(TIME_OF_DAY), !is.na(SEASON)) %>%
    distinct(STUDY_ID, BASELINE_AHA) %>%
    count(BASELINE_AHA, name = "n")
}

cohorts <- list(
  primary      = user,
  full_cohort  = user,   # the three strata's participants are the same in the full-cohort model
  htn_naive    = user %>% filter(HTN_MED_BASELINE == "never_logged"),
  glp1_naive   = user %>% filter(!GLP1_EVER),
  mid_excluded = user %>% filter(!(HTN_MED_BASELINE %in% "mid_study_initiator")),
  female       = user %>% filter(SEX == "female"),
  male         = user %>% filter(SEX == "male"),
  age_18_39    = user %>% filter(age_band == "18–39"),
  age_40_59    = user %>% filter(age_band == "40–59"),
  age_60_79    = user %>% filter(age_band == "60–79"))
n_tbl <- imap_dfr(cohorts, ~ adjusted_n(.x) %>% mutate(analysis = .y))

# Check the recount against n written by the model scripts
ref_n <- read_csv(file.path(A, "mid_excluded_sensitivity", "derived", "results", "5_1_mid_excluded_n_by_stratum.csv"),
                  show_col_types = FALSE) %>%
  filter(model == "sbp_adj") %>% transmute(analysis = cohort, BASELINE_AHA, n_ref = n_users)
chk <- inner_join(n_tbl, ref_n, by = c("analysis", "BASELINE_AHA"))
stopifnot(nrow(chk) == 6, all(chk$n == chk$n_ref))
cat("n recount matches model n (primary, mid_excluded): OK\n")

# --- Estimates: adjusted MO6 - BASELINE ---
from_lmm <- function(path, label) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(model %in% c("sbp_adj", "dbp_adj"), contrast == "MO6 - BASELINE") %>%
    transmute(analysis = label, outcome = toupper(sub("_adj", "", model)),
              BASELINE_AHA, estimate, lcl = asymp.LCL, ucl = asymp.UCL, p = p.value)
}
est <- bind_rows(
  from_lmm(file.path(A, "5_1", "derived", "results", "5_1_delta_contrasts.csv"), "primary"),
  from_lmm(file.path(A, "5_1", "Sensitivity", "full_cohort", "derived", "results", "delta_contrasts.csv"), "full_cohort") %>%
    filter(BASELINE_AHA %in% HYP),
  from_lmm(file.path(A, "5_1", "Sensitivity", "htn_naive", "derived", "results", "delta_contrasts.csv"), "htn_naive"),
  from_lmm(file.path(A, "glp1_naive_sensitivity", "derived", "results", "5_1_glp1_naive_delta_contrasts.csv"), "glp1_naive"),
  from_lmm(file.path(A, "mid_excluded_sensitivity", "derived", "results", "5_1_mid_excluded_delta_contrasts.csv"), "mid_excluded"),
  read_csv(file.path(A, "5_1", "Sensitivity", "age_sex", "derived", "results", "delta_contrasts_by_subgroup.csv"),
           show_col_types = FALSE) %>%
    filter(timepoint == "MO6") %>%
    transmute(analysis = recode(subgroup, "Female" = "female", "Male" = "male",
                                "Age 18–39" = "age_18_39", "Age 40–59" = "age_40_59", "Age 60–79" = "age_60_79"),
              outcome, BASELINE_AHA = HYP[match(stratum, SHORT)],
              estimate, lcl, ucl, p = p.value))
stopifnot(nrow(est) == 10 * 2 * 3, !anyNA(est$BASELINE_AHA))

LABELS <- c(primary      = "Primary cohort",
            full_cohort  = "Full cohort (including normal BP)",
            htn_naive    = "No recorded antihypertensive medication",
            glp1_naive   = "No recorded GLP-1 receptor agonist use",
            mid_excluded = "Excluding antihypertensive initiation during follow-up",
            female       = "Female",
            male         = "Male",
            age_18_39    = "Age 18–39 years",
            age_40_59    = "Age 40–59 years",
            age_60_79    = "Age 60–79 years")

long <- est %>%
  left_join(n_tbl, by = c("analysis", "BASELINE_AHA")) %>%
  mutate(analysis_label = LABELS[analysis],
         stratum = SHORT[match(BASELINE_AHA, HYP)]) %>%
  arrange(match(BASELINE_AHA, HYP), match(analysis, names(LABELS)), outcome)
write_csv(long, file.path(results_dir, "5_1_table_sensitivity.csv"))

fmt1 <- function(x) sub("^-", "−", sprintf("%.1f", round(x, 1) + 0))   # + 0 turns -0 into 0
cell <- function(e, l, u) sprintf("%s (%s, %s)", fmt1(e), fmt1(l), fmt1(u))
paper <- long %>%
  transmute(stratum, analysis = analysis_label, n = format(n, big.mark = ",", trim = TRUE),
            outcome, value = cell(estimate, lcl, ucl)) %>%
  pivot_wider(names_from = outcome, values_from = value) %>%
  select(stratum, analysis, n, SBP, DBP)
write_csv(paper, file.path(results_dir, "5_1_table_sensitivity_paper.csv"))
print(as.data.frame(paper), row.names = FALSE)
cat("\nWrote:", file.path(results_dir, "5_1_table_sensitivity.csv"), "\n")
cat("       ", file.path(results_dir, "5_1_table_sensitivity_paper.csv"), "\n")
