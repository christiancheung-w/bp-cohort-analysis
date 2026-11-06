# 8_table_sensitivity.R
# Supplementary table: separate-model engagement ORs for >=5 mmHg SBP and DBP
# reduction at 6 months, primary cohort vs three cohort restrictions.
# Refits the 7 x 2 single-predictor logistic models per cohort (specs as
# 8_engagement_logistic.R) so every CI is profile likelihood; ORs are checked
# against the four existing result files, CIs against the two that store them.
#
# Cohorts:
#   primary       §2.2 cohort (Elevated / Stage 1 / Stage 2, sex male/female)
#   htn_naive     HTN_MED_BASELINE == never_logged (HTN_BASELINE_ON dropped: constant)
#   glp1_naive    GLP1_EVER == FALSE
#   mid_excluded  HTN_MED_BASELINE != mid_study_initiator
#
# Inputs:
#   Build/derived/qc/user_analytic_frame_1X.rds, Build/derived/09_glp1_flag.rds
#   Analyses/8_engagement/derived/results/8_coefficients_separate_{primary,htn_naive}.csv
#   Analyses/glp1_naive_sensitivity/derived/results/8_glp1_naive_coefficients_separate.csv
#   Analyses/mid_excluded_sensitivity/derived/results/8_mid_excluded_coefficients_separate.csv
# Outputs (Analyses/8_engagement/derived/results/):
#   8_table_sensitivity.csv        - numeric, one row per cohort x outcome x predictor
#   8_table_sensitivity_paper.csv  - formatted OR (95% CI), * = p < 0.005

suppressPackageStartupMessages(library(tidyverse))
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")   # htn_baseline_on()

ALPHA <- 0.005
project_dir <- "./PRIMARY"
A <- file.path(project_dir, "Analyses")
results_dir <- file.path(A, "8_engagement", "derived", "results")

# Labels match 8_figure_forest.R; units are stated in the table footnote.
PREDS <- c(LOG_AVG_DAILY_TIME_IN_APP_MIN        = "Daily app time",
           PCT_DAYS_JOURNALED_per10             = "Journal completion",
           PCT_DAYS_ENGAGED_NOTIFICATIONS_per10 = "Notification engagement",
           PCT_DAYS_USED_COACH_per10            = "AI Coach engagement",
           AVG_SLEEP_NEED_PCT_per10             = "Sleep need achieved",
           PCT_DAYS_MET_OPTIMAL_STRAIN_per10    = "Physical activity recommendation",
           LOG_N_BP_READINGS_STUDY              = "Blood pressure readings")
COV      <- c("AGE_AT_START", "SEX", "BASELINE_BMI", "HTN_BASELINE_ON")
COV_HTNN <- c("AGE_AT_START", "SEX", "BASELINE_BMI")
COHORTS  <- c(primary      = "Primary cohort",
              htn_naive    = "No recorded antihypertensive medication",
              glp1_naive   = "No recorded GLP-1 receptor agonist use",
              mid_excluded = "Excluding antihypertensive initiation during follow-up")

base <- read_rds(file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")) %>%
  left_join(read_rds(file.path(project_dir, "Build", "derived", "09_glp1_flag.rds")) %>%
              select(STUDY_ID, GLP1_EVER), by = "STUDY_ID") %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female")) %>%
  mutate(DBP_RESPONDER_6MO = as.integer((MO6_DBP - BASELINE_DBP) <= -5),
         SEX = factor(SEX),
         PCT_DAYS_JOURNALED_per10             = PCT_DAYS_JOURNALED / 10,
         PCT_DAYS_ENGAGED_NOTIFICATIONS_per10 = PCT_DAYS_ENGAGED_NOTIFICATIONS / 10,
         PCT_DAYS_USED_COACH_per10            = PCT_DAYS_USED_COACH / 10,
         AVG_SLEEP_NEED_PCT_per10             = AVG_SLEEP_NEED_PCT / 10,
         PCT_DAYS_MET_OPTIMAL_STRAIN_per10    = PCT_DAYS_MET_OPTIMAL_STRAIN / 10,
         HTN_BASELINE_ON = htn_baseline_on(HTN_MED_BASELINE))
stopifnot(!anyNA(base$GLP1_EVER))

cohort_data <- list(
  primary      = list(d = base,                                                          cov = COV),
  htn_naive    = list(d = base %>% filter(HTN_MED_BASELINE == "never_logged"),           cov = COV_HTNN),
  glp1_naive   = list(d = base %>% filter(!GLP1_EVER),                                   cov = COV),
  mid_excluded = list(d = base %>% filter(!(HTN_MED_BASELINE %in% "mid_study_initiator")), cov = COV))

fit_one <- function(df, cov, otc, pred) {
  bcol <- paste0("BASELINE_", otc); resp <- paste0(otc, "_RESPONDER_6MO")
  d <- df %>% filter(if_all(all_of(c(pred, COV, bcol, resp)), ~ !is.na(.)))
  fit <- glm(as.formula(paste0(resp, " ~ ", paste(c(bcol, pred, cov), collapse = " + "))),
             family = binomial(), data = d)
  suppressMessages(broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)) %>%
    filter(term == pred) %>%
    transmute(outcome = otc, predictor = pred, OR = estimate, OR_lcl = conf.low, OR_ucl = conf.high,
              p.value, n = nrow(d))
}

long <- imap_dfr(cohort_data, function(cd, ck) {
  map_dfr(c("SBP", "DBP"), function(otc) map_dfr(names(PREDS), ~ fit_one(cd$d, cd$cov, otc, .x))) %>%
    mutate(cohort = ck)
})

# --- Fidelity against existing outputs ---
ref <- bind_rows(
  read_csv(file.path(results_dir, "8_coefficients_separate_primary.csv"), show_col_types = FALSE) %>%
    filter(term == predictor) %>% mutate(cohort = "primary"),
  read_csv(file.path(results_dir, "8_coefficients_separate_htn_naive.csv"), show_col_types = FALSE) %>%
    filter(term == predictor) %>% mutate(cohort = "htn_naive"),
  read_csv(file.path(A, "glp1_naive_sensitivity", "derived", "results", "8_glp1_naive_coefficients_separate.csv"),
           show_col_types = FALSE) %>% mutate(cohort = "glp1_naive"),
  read_csv(file.path(A, "mid_excluded_sensitivity", "derived", "results", "8_mid_excluded_coefficients_separate.csv"),
           show_col_types = FALSE) %>% mutate(cohort = "mid_excluded")) %>%
  select(cohort, outcome, predictor, OR_ref = OR, lcl_ref = any_of("OR_lcl"), ucl_ref = any_of("OR_ucl"), n_ref = n)
chk <- inner_join(long, ref, by = c("cohort", "outcome", "predictor"))
d_or <- max(abs(chk$OR - chk$OR_ref)); d_ci <- max(abs(c(chk$OR_lcl - chk$lcl_ref, chk$OR_ucl - chk$ucl_ref)), na.rm = TRUE)
cat(sprintf("Fidelity: %d of %d cells matched | max |OR diff| = %.2e | max |CI diff| (primary, HTN-naive) = %.2e | n identical: %s\n",
            nrow(chk), nrow(long), d_or, d_ci, all(chk$n == chk$n_ref)))
stopifnot(nrow(chk) == nrow(long), d_or < 1e-9, d_ci < 1e-9, all(chk$n == chk$n_ref))

long <- long %>% mutate(predictor_label = PREDS[predictor], cohort_label = COHORTS[cohort])
write_csv(long, file.path(results_dir, "8_table_sensitivity.csv"))

f2 <- function(x) sprintf("%.2f", x)
paper_or <- long %>%
  mutate(value = sprintf("%s (%s, %s)%s", f2(OR), f2(OR_lcl), f2(OR_ucl), ifelse(p.value < ALPHA, "*", "")),
         row = predictor_label)
paper_n <- long %>%
  group_by(cohort, outcome) %>%
  summarise(value = if (min(n) == max(n)) format(min(n), big.mark = ",")
                    else paste0(format(min(n), big.mark = ","), "–", format(max(n), big.mark = ",")),
            .groups = "drop") %>%
  mutate(row = "Participants, n")
paper <- bind_rows(paper_n, paper_or %>% select(cohort, outcome, row, value)) %>%
  mutate(cohort = COHORTS[cohort]) %>%
  pivot_wider(names_from = cohort, values_from = value) %>%
  mutate(outcome = factor(outcome, levels = c("SBP", "DBP")),
         row = factor(row, levels = c("Participants, n", unname(PREDS)))) %>%
  arrange(outcome, row) %>%
  mutate(across(c(outcome, row), as.character)) %>%
  select(outcome, row, all_of(unname(COHORTS)))
write_csv(paper, file.path(results_dir, "8_table_sensitivity_paper.csv"))
print(as.data.frame(paper), row.names = FALSE)
cat("\nWrote:", file.path(results_dir, "8_table_sensitivity.csv"), "\n")
cat("       ", file.path(results_dir, "8_table_sensitivity_paper.csv"), "\n")
