# 5_1_table_2_npj.R
# Table 2 for npj: Table 2 (5_1_table_2.R) plus the RTM-corrected 6-month change
# with bootstrap 95% CI, formatted to 1 decimal for BP. Reads CSVs only.
#
# Inputs:
#   Analyses/5_1/derived/results/5_1_table_2_within_stratum.csv
#   Analyses/5_1/Sensitivity/rtm_barnett_bootstrap/derived/results/rtm_bootstrap_summary.csv
#   Analyses/5_1/derived/results/5_1_r2.csv, 5_1_moments.csv   (model-fit rows)
# Output:
#   Analyses/5_1/derived/results/5_1_table_2_npj_paper.csv
#
# Columns: unadjusted and adjusted change (95% CI); RTM-corrected = adjusted
# change + expected RTM, bootstrap percentile 95% CI (6 months only);
# P = adjusted model, Holm-corrected; d = |adjusted change| / pooled baseline SD.

suppressPackageStartupMessages(library(tidyverse))

A <- "./PRIMARY/Analyses"
results_dir <- file.path(A, "5_1", "derived", "results")
HYP   <- c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")
SHORT <- c("Elevated", "Stage 1", "Stage 2")

t2  <- read_csv(file.path(results_dir, "5_1_table_2_within_stratum.csv"), show_col_types = FALSE)
rtm <- read_csv(file.path(A, "5_1", "Sensitivity", "rtm_barnett_bootstrap", "derived", "results",
                          "rtm_bootstrap_summary.csv"), show_col_types = FALSE) %>%
  filter(horizon == "MO6") %>%
  transmute(outcome, stratum = SHORT[match(stratum, HYP)], timepoint = "MO6",
            adj_check = delta, rtm_corrected, rtm_corrected_lcl, rtm_corrected_ucl)

tab <- t2 %>% left_join(rtm, by = c("outcome", "stratum", "timepoint"))
chk <- tab %>% filter(timepoint == "MO6")
stopifnot(nrow(chk) == 6, max(abs(chk$adj_delta - chk$adj_check)) < 1e-9)

fmt1 <- function(x) sub("^-", "−", sprintf("%.1f", round(x, 1) + 0))
ci   <- function(e, l, u) ifelse(is.na(e), "—", sprintf("%s (%s, %s)", fmt1(e), fmt1(l), fmt1(u)))
fmtp <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

paper <- tab %>%
  mutate(outcome = factor(outcome, levels = c("SBP", "DBP")),
         stratum = factor(stratum, levels = SHORT),
         months  = recode(timepoint, MO3 = "3", MO6 = "6")) %>%
  arrange(outcome, stratum, months) %>%
  transmute(outcome = as.character(outcome), stratum = as.character(stratum), months,
            n = format(n_paired, big.mark = ",", trim = TRUE),
            unadjusted    = ci(unadj_delta, unadj_lcl, unadj_ucl),
            adjusted      = ci(adj_delta, adj_lcl, adj_ucl),
            rtm_corrected = ci(rtm_corrected, rtm_corrected_lcl, rtm_corrected_ucl),
            p = fmtp(adj_p_holm),
            d = sprintf("%.2f", abs(d_baseline)))

# Model-fit rows under each outcome block: Nakagawa R² and residual moments,
# unadjusted and adjusted models in their own columns.
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
fits  <- fit_rows(file.path(results_dir, "5_1_r2.csv"), file.path(results_dir, "5_1_moments.csv"))
paper <- bind_rows(paper %>% filter(outcome == "SBP"), fits %>% filter(outcome == "SBP"),
                   paper %>% filter(outcome == "DBP"), fits %>% filter(outcome == "DBP"))
write_csv(paper, file.path(results_dir, "5_1_table_2_npj_paper.csv"))
print(as.data.frame(paper), row.names = FALSE)
cat("\nWrote:", file.path(results_dir, "5_1_table_2_npj_paper.csv"), "\n")
