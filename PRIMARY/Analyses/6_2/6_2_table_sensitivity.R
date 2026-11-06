# 6_2_table_sensitivity.R
# Supplementary table: wearable-metric logistic ORs for >=5 mmHg SBP and DBP
# reduction at 6 months, primary cohort vs four cohort variants.
# Refits the multivariable model per cohort (spec as 6_2_logistic.R) so every CI is
# profile likelihood; ORs and n are checked against the five existing result files.
#
# Cohorts:
#   primary       §2.2 cohort (Elevated / Stage 1 / Stage 2, sex male/female), BASELINE_N_DAYS >= 8
#   full          as primary plus Normal baseline BP
#   htn_naive     HTN_MED_BASELINE == never_logged (HTN_BASELINE_ON dropped: constant)
#   glp1_naive    GLP1_EVER == FALSE
#   mid_excluded  HTN_MED_BASELINE != mid_study_initiator
#
# Inputs:
#   Build/derived/qc/user_analytic_frame_1X.rds, Build/derived/09_glp1_flag.rds
#   Analyses/6_2/derived/results/6_2_coefficients_primary.csv, 6_2_sensitivity_{full_cohort,htn_naive}.csv
#   Analyses/glp1_naive_sensitivity/derived/results/6_2_glp1_naive_coefficients.csv
#   Analyses/mid_excluded_sensitivity/derived/results/6_2_mid_excluded_coefficients.csv
# Outputs (Analyses/6_2/derived/results/):
#   6_2_table_sensitivity.csv        - numeric, one row per cohort x outcome x predictor
#   6_2_table_sensitivity_paper.csv  - formatted OR (95% CI), * = p < 0.005

suppressPackageStartupMessages(library(tidyverse))
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")   # htn_baseline_on()

ALPHA <- 0.005
project_dir <- "./PRIMARY"
A <- file.path(project_dir, "Analyses")
results_dir <- file.path(A, "6_2", "derived", "results")
AHA3 <- c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")

PREDS <- c(RHR_BASELINE_MEAN = "RHR, baseline", RHR_STUDY_MEAN = "RHR, 6-month",
           HRV_BASELINE_MEAN = "HRV, baseline", HRV_STUDY_MEAN = "HRV, 6-month",
           SLEEP_DURATION_BASELINE_MEAN = "Sleep duration, baseline", SLEEP_DURATION_STUDY_MEAN = "Sleep duration, 6-month",
           SLEEP_CONSISTENCY_BASELINE_MEAN = "Sleep consistency, baseline", SLEEP_CONSISTENCY_STUDY_MEAN = "Sleep consistency, 6-month",
           STRAIN_BASELINE_MEAN = "Physical activity, baseline", STRAIN_STUDY_MEAN = "Physical activity, 6-month")
COV      <- c("AGE_AT_START", "SEX", "BASELINE_BMI", "HTN_BASELINE_ON")
COV_HTNN <- c("AGE_AT_START", "SEX", "BASELINE_BMI")
COHORTS  <- c(primary      = "Primary cohort",
              full         = "Full cohort (including normal BP)",
              htn_naive    = "No recorded antihypertensive medication",
              glp1_naive   = "No recorded GLP-1 receptor agonist use",
              mid_excluded = "Excluding antihypertensive initiation during follow-up")

base <- read_rds(file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")) %>%
  left_join(read_rds(file.path(project_dir, "Build", "derived", "09_glp1_flag.rds")) %>%
              select(STUDY_ID, GLP1_EVER), by = "STUDY_ID") %>%
  filter(SEX %in% c("male", "female"), BASELINE_N_DAYS >= 8) %>%
  mutate(DBP_RESPONDER_6MO = as.integer((MO6_DBP - BASELINE_DBP) <= -5),
         SEX = factor(SEX), HTN_BASELINE_ON = htn_baseline_on(HTN_MED_BASELINE))
stopifnot(!anyNA(base$GLP1_EVER))
hyp <- base %>% filter(BASELINE_AHA %in% AHA3)
cohort_data <- list(
  primary      = list(d = hyp,                                                          cov = COV),
  full         = list(d = base,                                                         cov = COV),
  htn_naive    = list(d = hyp %>% filter(HTN_MED_BASELINE == "never_logged"),           cov = COV_HTNN),
  glp1_naive   = list(d = hyp %>% filter(!GLP1_EVER),                                   cov = COV),
  mid_excluded = list(d = hyp %>% filter(!(HTN_MED_BASELINE %in% "mid_study_initiator")), cov = COV))

fit_one <- function(df, cov, otc) {
  bcol <- paste0("BASELINE_", otc); resp <- paste0(otc, "_RESPONDER_6MO")
  d <- df %>% filter(if_all(all_of(c(names(PREDS), COV, bcol, resp)), ~ !is.na(.)))
  fit <- glm(as.formula(paste0(resp, " ~ ", paste(c(bcol, names(PREDS), cov), collapse = " + "))),
             family = binomial(), data = d)
  suppressMessages(broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)) %>%
    filter(term %in% names(PREDS)) %>%
    transmute(outcome = otc, predictor = term, OR = estimate, OR_lcl = conf.low, OR_ucl = conf.high,
              p.value, n = nrow(d), n_events = sum(d[[resp]]))
}
long <- imap_dfr(cohort_data, function(cd, ck) map_dfr(c("SBP", "DBP"), ~ fit_one(cd$d, cd$cov, .x)) %>% mutate(cohort = ck))

# --- Fidelity against existing outputs (OR and n) ---
rd <- function(p) read_csv(p, show_col_types = FALSE)
ref <- bind_rows(
  rd(file.path(results_dir, "6_2_coefficients_primary.csv")) %>% mutate(cohort = "primary"),
  rd(file.path(results_dir, "6_2_sensitivity_full_cohort.csv")) %>% mutate(cohort = "full"),
  rd(file.path(results_dir, "6_2_sensitivity_htn_naive.csv")) %>% mutate(cohort = "htn_naive"),
  rd(file.path(A, "glp1_naive_sensitivity", "derived", "results", "6_2_glp1_naive_coefficients.csv")) %>% mutate(cohort = "glp1_naive"),
  rd(file.path(A, "mid_excluded_sensitivity", "derived", "results", "6_2_mid_excluded_coefficients.csv")) %>% mutate(cohort = "mid_excluded")) %>%
  filter(term %in% names(PREDS)) %>%
  transmute(cohort, outcome = sub("_responder", "", outcome), predictor = term, OR_ref = OR, n_ref = n)
chk <- inner_join(long, ref, by = c("cohort", "outcome", "predictor"))
d_or <- max(abs(chk$OR - chk$OR_ref))
cat(sprintf("Fidelity: %d of %d cells matched | max |OR diff| = %.2e | n identical: %s\n",
            nrow(chk), nrow(long), d_or, all(chk$n == chk$n_ref)))
stopifnot(nrow(chk) == nrow(long), d_or < 1e-9, all(chk$n == chk$n_ref))

# units check (printed for the footnote)
cat("Predictor means (primary cohort):\n")
print(hyp %>% summarise(across(all_of(names(PREDS)), ~ round(mean(.x, na.rm = TRUE), 2))) %>% pivot_longer(everything()), n = 20)

long <- long %>% mutate(predictor_label = PREDS[predictor], cohort_label = COHORTS[cohort])
write_csv(long, file.path(results_dir, "6_2_table_sensitivity.csv"))

f3 <- function(x) sprintf("%.3f", x)
paper_or <- long %>%
  mutate(value = sprintf("%s (%s, %s)%s", f3(OR), f3(OR_lcl), f3(OR_ucl), ifelse(p.value < ALPHA, "*", "")), row = predictor_label)
paper_n <- long %>% distinct(cohort, outcome, n, n_events) %>%
  transmute(cohort, outcome, row = "Participants (responders), n",
            value = sprintf("%s (%s)", format(n, big.mark = ",", trim = TRUE), format(n_events, big.mark = ",", trim = TRUE)))
paper <- bind_rows(paper_n, paper_or %>% select(cohort, outcome, row, value)) %>%
  mutate(cohort = COHORTS[cohort]) %>%
  pivot_wider(names_from = cohort, values_from = value) %>%
  mutate(outcome = factor(outcome, levels = c("SBP", "DBP")),
         row = factor(row, levels = c("Participants (responders), n", unname(PREDS)))) %>%
  arrange(outcome, row) %>% mutate(across(c(outcome, row), as.character)) %>%
  select(outcome, row, all_of(unname(COHORTS)))
write_csv(paper, file.path(results_dir, "6_2_table_sensitivity_paper.csv"))
print(as.data.frame(paper), row.names = FALSE)
cat("\nWrote:", file.path(results_dir, "6_2_table_sensitivity.csv"), "\n")
cat("       ", file.path(results_dir, "6_2_table_sensitivity_paper.csv"), "\n")
