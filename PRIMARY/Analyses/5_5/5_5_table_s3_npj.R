# 5_5_table_s3_npj.R
# Table S3 for npj: 12-month trajectory table (5_5_table_s3.R, four-timepoint
# LMM) plus the RTM-corrected 12-month change with bootstrap 95% CI, an n column
# and BP to 1 decimal. Reads CSVs and the user frame only.
#
# Inputs:
#   Analyses/5_5/derived/results/5_5_table_s3.csv
#   Analyses/5_1/Sensitivity/rtm_barnett_bootstrap/derived/results/rtm_bootstrap_summary.csv
#   Build/derived/qc/user_analytic_frame_1X.rds
#   Analyses/5_5/derived/results/5_5_r2.csv, 5_5_moments.csv   (model-fit rows)
# Output:
#   Analyses/5_5/derived/results/5_5_table_s3_npj_paper.csv
#
# n = participants with a baseline and a reading at that timepoint (SBP).
# RTM-corrected = adjusted 12-month change + expected RTM (reference: 12-month
# completers), bootstrap percentile 95% CI; 12-month rows only.

suppressPackageStartupMessages(library(tidyverse))

project_dir <- "./PRIMARY"
A <- file.path(project_dir, "Analyses")
results_dir <- file.path(A, "5_5", "derived", "results")
HYP   <- c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")
SHORT <- c("Elevated", "Stage 1", "Stage 2")
TP    <- c("3 months" = "MO3", "6 months" = "MO6", "12 months" = "MO12")

s3  <- read_csv(file.path(results_dir, "5_5_table_s3.csv"), show_col_types = FALSE)
rtm <- read_csv(file.path(A, "5_1", "Sensitivity", "rtm_barnett_bootstrap", "derived", "results",
                          "rtm_bootstrap_summary.csv"), show_col_types = FALSE) %>%
  filter(horizon == "MO12") %>%
  transmute(outcome, stratum = SHORT[match(stratum, HYP)], timepoint = "12 months",
            adj_check = delta, rtm_corrected, rtm_corrected_lcl, rtm_corrected_ucl)

u <- read_rds(file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")) %>%
  filter(BASELINE_AHA %in% HYP, SEX %in% c("male", "female"))
n_tbl <- map_dfr(names(TP), function(tp) {
  u %>% filter(!is.na(BASELINE_SBP), !is.na(.data[[paste0(TP[[tp]], "_SBP")]])) %>%
    count(BASELINE_AHA, name = "n") %>%
    transmute(stratum = SHORT[match(BASELINE_AHA, HYP)], timepoint = tp, n)
})

tab <- s3 %>%
  left_join(n_tbl, by = c("stratum", "timepoint")) %>%
  left_join(rtm, by = c("outcome", "stratum", "timepoint"))
chk <- tab %>% filter(timepoint == "12 months")
stopifnot(nrow(chk) == 6, max(abs(chk$adj_est - chk$adj_check)) < 1e-9, !anyNA(tab$n))

fmt1 <- function(x) sub("^-", "−", sprintf("%.1f", round(x, 1) + 0))
ci   <- function(e, l, u) ifelse(is.na(e), "—", sprintf("%s (%s, %s)", fmt1(e), fmt1(l), fmt1(u)))
fmtp <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

paper <- tab %>%
  mutate(outcome = factor(outcome, levels = c("SBP", "DBP")),
         stratum = factor(stratum, levels = SHORT),
         timepoint = factor(timepoint, levels = names(TP))) %>%
  arrange(outcome, stratum, timepoint) %>%
  transmute(outcome = as.character(outcome), stratum = as.character(stratum),
            months = sub(" months", "", as.character(timepoint)),
            n = format(n, big.mark = ",", trim = TRUE),
            unadjusted    = ci(unadj_est, unadj_lcl, unadj_ucl),
            adjusted      = ci(adj_est, adj_lcl, adj_ucl),
            rtm_corrected = ci(rtm_corrected, rtm_corrected_lcl, rtm_corrected_ucl),
            p = fmtp(p_holm),
            d = sprintf("%.2f", abs(d_baseline)))

# Model-fit rows under each outcome block (same layout as Table 2).
fit_rows <- function(r2_path, mom_path) {
  r2  <- read_csv(r2_path,  show_col_types = FALSE)
  mom <- read_csv(mom_path, show_col_types = FALSE) %>% filter(quantity == "residuals")
  map_dfr(c("SBP", "DBP"), function(o) {
    g <- function(tab, mod, col) tab[[col]][tab$model == paste0(tolower(o), "_", mod)]
    tibble(outcome = o,
           stratum = c("Marginal R²", "Conditional R²", "Residual skew / excess kurtosis"),
           months = "", n = "",
           unadjusted = c(sprintf("%.2f", g(r2, "unadj", "marginal_r2")), sprintf("%.2f", g(r2, "unadj", "conditional_r2")),
                          sprintf("%.2f / %.1f", g(mom, "unadj", "skewness"), g(mom, "unadj", "excess_kurtosis"))),
           adjusted   = c(sprintf("%.2f", g(r2, "adj", "marginal_r2")), sprintf("%.2f", g(r2, "adj", "conditional_r2")),
                          sprintf("%.2f / %.1f", g(mom, "adj", "skewness"), g(mom, "adj", "excess_kurtosis"))),
           rtm_corrected = "", p = "", d = "")
  })
}
fits  <- fit_rows(file.path(results_dir, "5_5_r2.csv"), file.path(results_dir, "5_5_moments.csv"))
paper <- bind_rows(paper %>% filter(outcome == "SBP"), fits %>% filter(outcome == "SBP"),
                   paper %>% filter(outcome == "DBP"), fits %>% filter(outcome == "DBP"))
write_csv(paper, file.path(results_dir, "5_5_table_s3_npj_paper.csv"))
print(as.data.frame(paper), row.names = FALSE)
cat("\nWrote:", file.path(results_dir, "5_5_table_s3_npj_paper.csv"), "\n")
