# 7_2_table_sensitivity.R
# Supplementary table: reading-level within/between-person coefficients (beta, mmHg per
# 1 SD of the predictor) in the primary cohort and five sensitivity variants.
# Assembled from the frozen coefficient files; no refitting. All CIs are Wald intervals
# from the cluster-robust (CR1) z statistics, as in 7_2_wb_lmm.R:
#   primary / RawBP / FullCohort / HTN_naive: conf_low, conf_high stored in the files
#   glp1_naive / mid_excluded: SE = estimate / z, CI = estimate +/- 1.959964 * SE
# Cohort sizes for glp1_naive and mid_excluded are recounted here with the same
# cohort rule as the fitting scripts (>=2 non-baseline DURING readings, complete cases).
#
# Outputs (Analyses/7_2/derived/results/):
#   7_2_table_sensitivity.csv        - numeric, one row per cohort x outcome x term
#   7_2_table_sensitivity_paper.csv  - formatted beta (95% CI), * = p < 0.005

suppressPackageStartupMessages(library(tidyverse))
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")   # htn_baseline_on()

ALPHA <- 0.005
project_dir <- "./PRIMARY"
A <- file.path(project_dir, "Analyses")
results_dir <- file.path(A, "7_2", "derived", "results")
AHA3 <- c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")
Z <- qnorm(0.975)

METRICS <- c(RHR = "RHR", HRV = "HRV", SleepDuration = "Sleep duration",
             SleepConsistency = "Sleep consistency", Strain = "Physical activity")
TERMS <- c(paste0(names(METRICS), "_b"), paste0(names(METRICS), "_w"))
term_label <- function(t) sprintf("%s, %s-person", METRICS[sub("_[bw]$", "", t)], ifelse(grepl("_b$", t), "between", "within"))
COHORTS <- c(primary      = "Primary cohort",
             rawbp        = "Raw BP outcome",
             full         = "Full cohort (including normal BP)",
             htn_naive    = "No recorded antihypertensive medication",
             glp1_naive   = "No recorded GLP-1 receptor agonist use",
             mid_excluded = "Excluding antihypertensive initiation during follow-up")

rd <- function(p) read_csv(p, show_col_types = FALSE)
otc <- function(x) substr(x, 1, 3)   # "SBP-baseline" / "SBP" -> "SBP"

prim <- rd(file.path(results_dir, "7_2_coefficients.csv")) %>%
  transmute(cohort = "primary", outcome = otc(outcome), term, estimate, lcl = conf_low, ucl = conf_high, p = p_value)
sens <- rd(file.path(results_dir, "7_2_sensitivity_coefficients.csv")) %>%
  mutate(cohort = c(RawBP = "rawbp", FullCohort = "full", HTN_naive = "htn_naive")[sensitivity]) %>%
  transmute(cohort, outcome = otc(outcome), term, estimate, lcl = conf_low, ucl = conf_high, p = p_value, n_users, n_obs)
wald <- function(p, ck) rd(p) %>%
  transmute(cohort = ck, outcome = otc(outcome), term, estimate, se = estimate / z_stat,
            lcl = estimate - Z * se, ucl = estimate + Z * se, p = p_value) %>% select(-se)
glp <- wald(file.path(A, "glp1_naive_sensitivity", "derived", "results", "7_2_glp1_naive_coefficients.csv"), "glp1_naive")
mid <- wald(file.path(A, "mid_excluded_sensitivity", "derived", "results", "7_2_mid_excluded_coefficients.csv"), "mid_excluded")

# consistency check of the Wald construction against the stored CIs (primary file has both)
pchk <- rd(file.path(results_dir, "7_2_coefficients.csv")) %>%
  mutate(lcl_w = estimate - Z * SE_robust, ucl_w = estimate + Z * SE_robust)
cat(sprintf("Wald-CI construction vs stored CI (primary): max |diff| = %.2e\n",
            max(abs(c(pchk$lcl_w - pchk$conf_low, pchk$ucl_w - pchk$conf_high)))))
stopifnot(max(abs(c(pchk$lcl_w - pchk$conf_low, pchk$ucl_w - pchk$conf_high))) < 1e-9)

long <- bind_rows(prim, sens %>% select(-n_users, -n_obs), glp, mid) %>% filter(term %in% TERMS)
stopifnot(nrow(long) == 6 * 2 * 10)

# --- cohort sizes ---
sizes_known <- bind_rows(
  rd(file.path(results_dir, "7_2_sizes_std.csv")) %>% transmute(cohort = "primary", outcome, n_users, n_obs),
  rd(file.path(results_dir, "7_2_sensitivities_sizes_std.csv")) %>%
    mutate(cohort = c("W/B daily LMM (raw-BP)" = "rawbp", "W/B daily LMM (full-cohort)" = "full",
                      "W/B daily LMM (HTN-naive)" = "htn_naive")[analysis]) %>% transmute(cohort, outcome, n_users, n_obs)) %>%
  distinct(cohort, n_users, n_obs)
daily72 <- c(RHR = "RHR_LOOKBACK_MEAN", HRV = "HRV_LOOKBACK_MEAN", Strain = "STRAIN_LOOKBACK_MEAN",
             SleepDuration = "SLEEP_DURATION_LOOKBACK_MEAN", SleepConsistency = "SLEEP_CONSISTENCY_LOOKBACK_MEAN")
bp <- read_rds(file.path(project_dir, "Build", "derived", "qc", "bp_reading_analytic_frame_1X.rds")) %>%
  left_join(read_rds(file.path(project_dir, "Build", "derived", "09_glp1_flag.rds")) %>% select(STUDY_ID, GLP1_EVER), by = "STUDY_ID")
stopifnot(!anyNA(bp$GLP1_EVER), !anyNA(bp$HTN_MED_BASELINE))
count72 <- function(df) {
  d <- df %>% filter(BASELINE_AHA %in% AHA3, SEX %in% c("male", "female"), !IS_BASELINE_READING, PERIOD == "DURING") %>%
    group_by(STUDY_ID) %>% filter(n() >= 2) %>% ungroup() %>%
    mutate(HTN_BASELINE_ON = htn_baseline_on(HTN_MED_BASELINE)) %>%
    filter(if_all(all_of(c(unname(daily72), "AGE_AT_START", "SEX", "BASELINE_BMI", "HTN_BASELINE_ON",
                           "SBP_MINUS_BASELINE", "DBP_MINUS_BASELINE")), ~ !is.na(.)))
  tibble(n_users = n_distinct(d$STUDY_ID), n_obs = nrow(d))
}
n_prim_chk <- count72(bp)
cat(sprintf("Recount check, primary: users %d (file %d), readings %d (file %d)\n",
            n_prim_chk$n_users, sizes_known$n_users[sizes_known$cohort == "primary"],
            n_prim_chk$n_obs,   sizes_known$n_obs[sizes_known$cohort == "primary"]))
stopifnot(n_prim_chk$n_users == sizes_known$n_users[sizes_known$cohort == "primary"],
          n_prim_chk$n_obs   == sizes_known$n_obs[sizes_known$cohort == "primary"])
sizes <- bind_rows(sizes_known,
                   count72(bp %>% filter(!GLP1_EVER)) %>% mutate(cohort = "glp1_naive"),
                   count72(bp %>% filter(!(HTN_MED_BASELINE %in% "mid_study_initiator"))) %>% mutate(cohort = "mid_excluded"))

long <- long %>% left_join(sizes, by = "cohort") %>%
  mutate(term_label = term_label(term), cohort_label = COHORTS[cohort])
write_csv(long, file.path(results_dir, "7_2_table_sensitivity.csv"))

f2 <- function(x) { s <- sprintf("%.2f", x); s[s == "-0.00"] <- "0.00"; sub("^-", "\u2212", s) }
paper_b <- long %>%
  mutate(value = sprintf("%s (%s, %s)%s", f2(estimate), f2(lcl), f2(ucl), ifelse(p < ALPHA, "*", "")), row = term_label) %>%
  select(cohort, outcome, row, value)
paper_n <- sizes %>% crossing(outcome = c("SBP", "DBP")) %>%
  transmute(cohort, outcome, row = "Participants (readings), n",
            value = sprintf("%s (%s)", format(n_users, big.mark = ",", trim = TRUE), format(n_obs, big.mark = ",", trim = TRUE)))
row_levels <- c("Participants (readings), n", term_label(TERMS))
paper <- bind_rows(paper_n, paper_b) %>%
  mutate(cohort = COHORTS[cohort]) %>%
  pivot_wider(names_from = cohort, values_from = value) %>%
  mutate(outcome = factor(outcome, levels = c("SBP", "DBP")), row = factor(row, levels = row_levels)) %>%
  arrange(outcome, row) %>% mutate(across(c(outcome, row), as.character)) %>%
  select(outcome, row, all_of(unname(COHORTS)))
write_csv(paper, file.path(results_dir, "7_2_table_sensitivity_paper.csv"))
print(as.data.frame(paper), row.names = FALSE)
cat("\nWrote:", file.path(results_dir, "7_2_table_sensitivity.csv"), "\n")
cat("       ", file.path(results_dir, "7_2_table_sensitivity_paper.csv"), "\n")
