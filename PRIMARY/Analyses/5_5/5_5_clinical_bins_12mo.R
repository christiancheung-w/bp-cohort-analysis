# 5_5_clinical_bins_12mo.R
# SAP §5.5 companion (parallels §5.2 at MO12).
#
# Per-user 3-bin classification of ΔBP (MO12 − BASELINE), separately for SBP
# and DBP, tabulated by baseline AHA stratum:
#   improved:  Δ ≤ −5 mmHg
#   stable:    −5 < Δ < +5 mmHg
#   worsened:  Δ ≥ +5 mmHg
#
# Cohort: §2.2 hypertensive (Elevated + Stage 1 + Stage 2), restricted to the
# MO12 subset (users with both BASELINE and MO12 BP for the relevant outcome).
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/5_5/derived/):
#   results/5_5_clinical_bins_12mo_counts.csv
#   results/5_5_clinical_bins_12mo_percentages.csv
#   results/5_5_clinical_bins_12mo_paper.md
#   results/5_5_clinical_bins_12mo.log

library(tidyverse)

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "5_5", "derived")
results_dir <- file.path(out_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "5_5_clinical_bins_12mo.log")
sink(log_path, split = TRUE)
cat("§5.5 clinical bins (MO12) —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

user_frame <- read_rds(input_path)
cohort <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"))
cat("§2.2 cohort:", nrow(cohort), "users\n")

binned <- cohort %>%
  mutate(
    DELTA_SBP_MO12 = MO12_SBP - BASELINE_SBP,
    DELTA_DBP_MO12 = MO12_DBP - BASELINE_DBP,
    bin_SBP = case_when(
      is.na(DELTA_SBP_MO12) ~ NA_character_,
      DELTA_SBP_MO12 <= -5  ~ "improved",
      DELTA_SBP_MO12 >=  5  ~ "worsened",
      TRUE                  ~ "stable"
    ),
    bin_DBP = case_when(
      is.na(DELTA_DBP_MO12) ~ NA_character_,
      DELTA_DBP_MO12 <= -5  ~ "improved",
      DELTA_DBP_MO12 >=  5  ~ "worsened",
      TRUE                  ~ "stable"
    ),
    stratum = factor(BASELINE_AHA,
                     levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
                     labels = c("Elevated", "Stage 1", "Stage 2"))
  )

bin_long <- bind_rows(
  binned %>% filter(!is.na(bin_SBP)) %>%
    count(stratum, bin = bin_SBP) %>% mutate(outcome = "SBP"),
  binned %>% filter(!is.na(bin_DBP)) %>%
    count(stratum, bin = bin_DBP) %>% mutate(outcome = "DBP")
) %>%
  mutate(bin = factor(bin, levels = c("improved", "stable", "worsened"))) %>%
  group_by(outcome, stratum) %>%
  mutate(n_total = sum(n),
         prop    = n / n_total) %>%
  ungroup() %>%
  arrange(outcome, stratum, bin)

wilson_ci <- function(x, n, z = 1.959964) {
  p <- x / n
  denom  <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  list(lcl = center - half, ucl = center + half)
}
bin_long <- bin_long %>%
  rowwise() %>%
  mutate(ci = list(wilson_ci(n, n_total)),
         lcl = ci$lcl, ucl = ci$ucl) %>%
  ungroup() %>%
  select(-ci) %>%
  select(outcome, stratum, bin, n, n_total, prop, lcl, ucl)

write_csv(bin_long, file.path(results_dir, "5_5_clinical_bins_12mo_counts.csv"))

bin_wide <- bin_long %>%
  mutate(cell = sprintf("%d / %d (%.1f%% [%.1f, %.1f])",
                        n, n_total, 100*prop, 100*lcl, 100*ucl)) %>%
  select(outcome, stratum, bin, cell) %>%
  pivot_wider(names_from = bin, values_from = cell) %>%
  arrange(outcome, stratum)

write_csv(bin_wide, file.path(results_dir, "5_5_clinical_bins_12mo_percentages.csv"))

cat("\n--- Bin tabulation at MO12 ---\n")
print(bin_long, n = Inf, width = Inf)

# --- Paper-ready markdown narrative ---
narrative_row <- function(outcome_label, stratum_label) {
  r <- bin_long %>% filter(outcome == outcome_label, stratum == stratum_label)
  if (nrow(r) == 0) return("")
  improved <- r %>% filter(bin == "improved")
  stable   <- r %>% filter(bin == "stable")
  worsened <- r %>% filter(bin == "worsened")
  sprintf(
    "%s × %s: improved %.1f%% (n = %d), stable %.1f%% (n = %d), worsened %.1f%% (n = %d). Total n = %d.",
    outcome_label, stratum_label,
    100*improved$prop, improved$n,
    100*stable$prop,   stable$n,
    100*worsened$prop, worsened$n,
    improved$n_total
  )
}

md_lines <- c(
  "# §5.5 — Clinical-threshold bin proportions at 12 months",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "## Per-stratum bin percentages (MO12 subset)",
  ""
)
for (otc in c("SBP", "DBP")) {
  md_lines <- c(md_lines, paste0("### ", otc))
  for (st in c("Elevated", "Stage 1", "Stage 2")) {
    md_lines <- c(md_lines, paste0("- ", narrative_row(otc, st)))
  }
  md_lines <- c(md_lines, "")
}
md_lines <- c(md_lines,
  "## Notes",
  "",
  "- Bin thresholds at ±5 mmHg (per SAP §5.2; same as MO6 companion).",
  "- Cohort = §2.2 hypertensive, restricted to MO12 subset (~7% of full cohort). Per-stratum n reflects the analyzable MO12 paired subset.",
  "- MO12 subset is selection-enriched for sustained logging engagement. Results are conditional on this subset and not directly comparable to the full-cohort MO6 percentages without acknowledging the selection.",
  "- Per SAP, no formal hypothesis test. Wilson 95% CIs reported.",
  "- Pairs with §5.2 (MO6 clinical bins) for the 6-vs-12-month comparison."
)
writeLines(md_lines, file.path(results_dir, "5_5_clinical_bins_12mo_paper.md"))

cat("\nWrote markdown:\n  ", file.path(results_dir, "5_5_clinical_bins_12mo_paper.md"), "\n")
cat("\n=== §5.5 12mo clinical bins done ===\n")
sink()
