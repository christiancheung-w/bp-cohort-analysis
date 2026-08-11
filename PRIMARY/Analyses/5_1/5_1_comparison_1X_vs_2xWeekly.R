# 5_1_comparison_1X_vs_2xWeekly.R
# Side-by-side comparison of §5.1 primary (1X cohort) vs §5.1 sensitivity
# (2xWeekly cohort). Produces a single table with per-cell Δ from each cohort,
# attenuation ratio (2xWeekly / 1X), and a concordance flag.
#
# The point of this output: a single artifact the manuscript can reference to
# document the RTM-related attenuation pattern without re-running anything.
#
# Inputs:
#   Analyses/5_1/derived/emmeans_results.rds                          (1X primary)
#   Analyses/5_1/Sensitivity/2xWeekly/derived/emmeans_results.rds     (2xWeekly)
#
# Outputs (Analyses/5_1/derived/results/):
#   5_1_comparison_1X_vs_2xWeekly.csv
#   5_1_comparison_paper.md

library(tidyverse)

project_dir   <- "./PRIMARY"
primary_path  <- file.path(project_dir, "Analyses", "5_1", "derived",
                           "emmeans_results.rds")
sens_path     <- file.path(project_dir, "Analyses", "5_1", "Sensitivity",
                           "2xWeekly", "derived", "emmeans_results.rds")
results_dir   <- file.path(project_dir, "Analyses", "5_1", "derived", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

primary_emm <- read_rds(primary_path)
sens_emm    <- read_rds(sens_path)

# --- Helper: derive Wald CIs from estimate ± qt(0.975, df)*SE (z when df=Inf) ---
ensure_ci <- function(dat) {
  nm <- names(dat)
  if (all(c("asymp.LCL", "asymp.UCL") %in% nm)) {
    dat$.lcl <- dat$asymp.LCL; dat$.ucl <- dat$asymp.UCL
  } else if (all(c("lower.CL", "upper.CL") %in% nm)) {
    dat$.lcl <- dat$lower.CL;  dat$.ucl <- dat$upper.CL
  } else if (all(c("estimate", "SE") %in% nm)) {
    crit <- if ("df" %in% nm) ifelse(is.finite(dat$df), qt(0.975, dat$df), 1.959964) else 1.959964
    dat$.lcl <- dat$estimate - crit * dat$SE
    dat$.ucl <- dat$estimate + crit * dat$SE
  } else {
    stop("Could not derive CI from columns: ", paste(nm, collapse = ", "))
  }
  dat
}

# --- Long-format extraction from one emmeans_results list ---
extract_long <- function(emm_list, cohort_label) {
  bind_rows(
    extract_within(emm_list$sbp_unadj$within_stratum, "SBP", "unadj"),
    extract_within(emm_list$sbp_adj$within_stratum,   "SBP", "adj"),
    extract_within(emm_list$dbp_unadj$within_stratum, "DBP", "unadj"),
    extract_within(emm_list$dbp_adj$within_stratum,   "DBP", "adj")
  ) %>%
    mutate(cohort = cohort_label)
}

extract_within <- function(tib, outcome_label, model_label) {
  tib <- ensure_ci(tib)
  tibble(
    outcome   = outcome_label,
    stratum   = factor(tib$BASELINE_AHA,
                       levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
                       labels = c("Elevated", "Stage 1", "Stage 2")),
    timepoint = sub(" - BASELINE", "", as.character(tib$contrast)),
    model     = model_label,
    delta     = tib$estimate,
    lcl       = tib$.lcl,
    ucl       = tib$.ucl,
    p         = tib$p.holm   # Holm-adjusted within family (v0.79); drives the concordance logic
  )
}

primary_long <- extract_long(primary_emm, "1X")
sens_long    <- extract_long(sens_emm,    "2xWeekly")

# --- Join into side-by-side wide ---
comp_wide <- bind_rows(primary_long, sens_long) %>%
  pivot_wider(names_from = cohort,
              values_from = c(delta, lcl, ucl, p),
              names_glue  = "{cohort}_{.value}") %>%
  mutate(
    # Retention = how much of 1X magnitude is preserved in 2xWeekly.
    # Only meaningful when 1X magnitude is non-trivial (|1X_delta| > 0.5 mmHg).
    abs_1X        = abs(`1X_delta`),
    retention_pct = ifelse(abs_1X > 0.5,
                           round(100 * abs(`2xWeekly_delta`) / abs_1X, 1),
                           NA_real_),
    # Concordance: same sign AND both significant at α = 0.005, OR both NS.
    same_sign     = sign(`1X_delta`) == sign(`2xWeekly_delta`),
    both_sig      = `1X_p` < 0.005 & `2xWeekly_p` < 0.005,
    both_ns       = `1X_p` >= 0.005 & `2xWeekly_p` >= 0.005,
    concordant    = (same_sign & both_sig) | both_ns
  ) %>%
  select(outcome, stratum, timepoint, model,
         `1X_delta`, `1X_lcl`, `1X_ucl`, `1X_p`,
         `2xWeekly_delta`, `2xWeekly_lcl`, `2xWeekly_ucl`, `2xWeekly_p`,
         retention_pct, concordant) %>%
  arrange(model, outcome, stratum, timepoint)

write_csv(comp_wide, file.path(results_dir, "5_1_comparison_1X_vs_2xWeekly.csv"))

# --- Paper-ready markdown rendering (adjusted only — primary reporting frame) ---
fmt_ci <- function(est, lcl, ucl, digits = 2) {
  sprintf("%.*f (%.*f, %.*f)", digits, est, digits, lcl, digits, ucl)
}
fmt_p <- function(p) {
  ifelse(p < 1e-300, "<1e-300",
  ifelse(p < 1e-10,  formatC(p, format = "e", digits = 1),
  ifelse(p < 0.001,  formatC(p, format = "e", digits = 1),
                     sprintf("%.3f", p))))
}

adj_md <- comp_wide %>%
  filter(model == "adj") %>%
  mutate(
    `1X Δ (95% CI)`       = fmt_ci(`1X_delta`,       `1X_lcl`,       `1X_ucl`),
    `2xWeekly Δ (95% CI)` = fmt_ci(`2xWeekly_delta`, `2xWeekly_lcl`, `2xWeekly_ucl`),
    `Retention %`         = ifelse(is.na(retention_pct), "—",
                                   sprintf("%.0f%%", retention_pct)),
    Concordant            = ifelse(concordant, "✓", "✗"),
    `1X p`                = fmt_p(`1X_p`),
    `2xWeekly p`          = fmt_p(`2xWeekly_p`)
  ) %>%
  select(Outcome = outcome, Stratum = stratum, Timepoint = timepoint,
         `1X Δ (95% CI)`, `1X p`,
         `2xWeekly Δ (95% CI)`, `2xWeekly p`,
         `Retention %`, Concordant)

unadj_md <- comp_wide %>%
  filter(model == "unadj") %>%
  mutate(
    `1X Δ (95% CI)`       = fmt_ci(`1X_delta`,       `1X_lcl`,       `1X_ucl`),
    `2xWeekly Δ (95% CI)` = fmt_ci(`2xWeekly_delta`, `2xWeekly_lcl`, `2xWeekly_ucl`),
    `Retention %`         = ifelse(is.na(retention_pct), "—",
                                   sprintf("%.0f%%", retention_pct)),
    Concordant            = ifelse(concordant, "✓", "✗"),
    `1X p`                = fmt_p(`1X_p`),
    `2xWeekly p`          = fmt_p(`2xWeekly_p`)
  ) %>%
  select(Outcome = outcome, Stratum = stratum, Timepoint = timepoint,
         `1X Δ (95% CI)`, `1X p`,
         `2xWeekly Δ (95% CI)`, `2xWeekly p`,
         `Retention %`, Concordant)

md_lines <- c(
  "# §5.1 sensitivity — 1X primary vs 2xWeekly side-by-side",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "Per-cell Δ from baseline in mmHg with 95% CIs, side-by-side. Retention = |2xWeekly Δ| / |1X Δ| (suppressed when |1X Δ| < 0.5 mmHg, i.e., when the primary estimate is itself substantively zero). Concordance: same sign AND both p < 0.005, or both p ≥ 0.005.",
  "",
  "## Adjusted LMM (primary reporting frame)",
  "",
  knitr::kable(adj_md, format = "markdown", align = "lllllllrl"),
  "",
  "## Unadjusted LMM (sensitivity to covariate spec)",
  "",
  knitr::kable(unadj_md, format = "markdown", align = "lllllllrl"),
  "",
  "## Notes",
  "",
  "- 1X cohort: n = 20,021 users (§2.2 hypertensive filter). 2xWeekly cohort: n = 2,190 users (§2.2 filter applied to the 2xWeekly eligibility subset).",
  "- p-values shown are model-internal (from each cohort's respective LMM). At α = 0.005, both cohorts achieve significance on every non-null cell (the two Elevated DBP cells are NS in both, consistent with the no-headroom interpretation).",
  "- Retention % indicates how much of the 1X estimate's magnitude is preserved in the 2xWeekly cohort. Lower retention = greater attenuation = more RTM-related inflation suspected in 1X. The gradient (lowest retention at Elevated, highest at Stage 2) is consistent with RTM contributing proportionally more to smaller true effects.",
  "- This is a robustness check, not a magnitude correction. The 2xWeekly cohort is itself selection-biased (sustained logging engagement) and its estimates do not represent ground truth. The interpretive frame is 'direction-concordant; 1X magnitudes bounded below by 2xWeekly,' not '2xWeekly is the true effect.'"
)
writeLines(md_lines, file.path(results_dir, "5_1_comparison_paper.md"))

cat("Wrote:\n")
cat("  ", file.path(results_dir, "5_1_comparison_1X_vs_2xWeekly.csv"), "\n")
cat("  ", file.path(results_dir, "5_1_comparison_paper.md"), "\n")
cat("\n--- Adjusted comparison preview ---\n")
print(comp_wide %>% filter(model == "adj"), n = Inf, width = Inf)
