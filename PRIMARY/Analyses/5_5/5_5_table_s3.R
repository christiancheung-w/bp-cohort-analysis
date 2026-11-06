# 5_5_table_s3.R
# SAP §5.5 Supplementary Table S3: per-stratum change from baseline at 3, 6, and 12 months.
#
# Reads from already-computed CSVs — does NOT refit models.
#
# Columns: Outcome | Stratum | Timepoint | Unadjusted Δ (95% CI) |
#          Adjusted Δ (95% CI) | p.holm (vs. baseline) | Cohen's d
#
# Cohen's d = adjusted LMM Δ / pooled baseline SD of the outcome (all strata
# combined, one value per outcome), following Kazis et al. (Med Care 1989).
# Pooled rather than per-stratum to avoid range restriction from baseline-defined
# strata; baseline SD uses fully-observed baseline BP, robust to MO12 attrition.
# Signed (negative = BP reduction). Harmonized with §5.1 Table 2.
#
# Input:  derived/results/5_5_delta_contrasts.csv
#         Build/derived/qc/user_analytic_frame_1X.rds  (baseline SD denominator)
#
# Output: derived/results/5_5_table_s3.csv        — tidy numeric
#         derived/results/5_5_table_s3_paper.csv   — formatted for manuscript

library(tidyverse)

project_dir <- "./PRIMARY"
results_dir <- file.path(project_dir, "Analyses", "5_5", "derived", "results")

log_path <- file.path(results_dir, "5_5_table_s3.log")
sink(log_path, split = TRUE)
cat("§5.5 Table S3 —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read inputs ---
delta   <- read_csv(file.path(results_dir, "5_5_delta_contrasts.csv"),
                    show_col_types = FALSE)
uf_path    <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
user_frame <- read_rds(uf_path)

cat("delta_contrasts rows:", nrow(delta), "\n")
cat("user_frame rows:", nrow(user_frame), "\n\n")

# --- Cohort baseline SD: denominator for Cohen's d (one value per outcome) ---
# Pooled across strata to avoid range restriction from baseline-defined strata;
# uses fully-observed baseline BP (robust to MO12 attrition). Kazis (1989).
baseline_sd <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female")) %>%
  summarise(SBP = sd(BASELINE_SBP, na.rm = TRUE),
            DBP = sd(BASELINE_DBP, na.rm = TRUE)) %>%
  unlist()

cat("Cohort baseline SD (Cohen's d denominator):\n")
cat("  SBP:", round(baseline_sd[["SBP"]], 4), " DBP:", round(baseline_sd[["DBP"]], 4), "\n\n")

# Label maps
tp_map <- c("MO3 - BASELINE"  = "3 months",
            "MO6 - BASELINE"  = "6 months",
            "MO12 - BASELINE" = "12 months")

stratum_map <- c("Elevated"             = "Elevated",
                 "Stage 1 Hypertension" = "Stage 1",
                 "Stage 2 Hypertension" = "Stage 2")

target_contrasts <- names(tp_map)

# --- Unadjusted estimates ---
unadj <- delta %>%
  filter(model %in% c("sbp_unadj", "dbp_unadj"),
         contrast %in% target_contrasts) %>%
  transmute(
    outcome   = if_else(startsWith(model, "sbp"), "SBP", "DBP"),
    stratum   = stratum_map[BASELINE_AHA],
    timepoint = tp_map[contrast],
    unadj_est = estimate,
    unadj_lcl = asymp.LCL,
    unadj_ucl = asymp.UCL
  )

# --- Adjusted estimates + p-values ---
adj <- delta %>%
  filter(model %in% c("sbp_adj", "dbp_adj"),
         contrast %in% target_contrasts) %>%
  transmute(
    outcome    = if_else(startsWith(model, "sbp"), "SBP", "DBP"),
    stratum    = stratum_map[BASELINE_AHA],
    timepoint  = tp_map[contrast],
    adj_est    = estimate,
    adj_lcl    = asymp.LCL,
    adj_ucl    = asymp.UCL,
    p_holm     = p.holm,
    d_baseline = estimate / if_else(startsWith(model, "sbp"),
                                    baseline_sd[["SBP"]], baseline_sd[["DBP"]])
  )

# --- Join and order ---
tbl <- left_join(unadj, adj, by = c("outcome", "stratum", "timepoint")) %>%
  mutate(
    outcome   = factor(outcome,   levels = c("SBP", "DBP")),
    stratum   = factor(stratum,   levels = c("Elevated", "Stage 1", "Stage 2")),
    timepoint = factor(timepoint, levels = c("3 months", "6 months", "12 months"))
  ) %>%
  arrange(outcome, stratum, timepoint)

cat("=== Tidy numeric table ===\n")
print(tbl, n = Inf, width = Inf)

write_csv(tbl, file.path(results_dir, "5_5_table_s3.csv"))
cat("\nWrote:", file.path(results_dir, "5_5_table_s3.csv"), "\n")

# --- Formatted paper version ---
fmt_change <- function(est, lcl, ucl) {
  sprintf("%+.2f (%+.2f, %+.2f)", est, lcl, ucl)
}
fmt_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}

paper <- tbl %>%
  transmute(
    Outcome              = as.character(outcome),
    Stratum              = as.character(stratum),
    Timepoint            = as.character(timepoint),
    `Unadjusted Δ (mmHg)` = fmt_change(unadj_est, unadj_lcl, unadj_ucl),
    `Adjusted Δ (mmHg)`   = fmt_change(adj_est,   adj_lcl,   adj_ucl),
    `p (vs. baseline)`    = fmt_p(p_holm),
    `Cohen's d`           = sprintf("%.3f", d_baseline)
  )

cat("\n=== Paper-formatted table ===\n")
print(paper, n = Inf, width = Inf)

write_csv(paper, file.path(results_dir, "5_5_table_s3_paper.csv"))
cat("\nWrote:", file.path(results_dir, "5_5_table_s3_paper.csv"), "\n")

cat("\n=== §5.5 Table S3 done ===\n")
sink()
