# 5_1_table_2.R
# SAP §5.1 Table 2: per-stratum Δ in mmHg (unadjusted + adjusted side-by-side),
# Cohen's d (baseline-SD standardized effect size, supplementary), and cross-
# stratum Δ-difference contrasts. For SBP and DBP.
#
# p-values reported from the ADJUSTED LMM (the fuller spec incl. baseline-HTN),
# Holm-adjusted within contrast family. The within-stratum Δ is led by its 95% CI
# (estimation-first). Raw + unadjusted-model + adjusted p-values are all exported
# in the CSVs for transparency.
#
# Cohen's d definition used here (harmonized with §5.5 Table S3):
#   d = adjusted LMM Δ (mmHg) / pooled baseline SD of the outcome
# Denominator = SD of baseline BP across the full analytic cohort (all strata
# pooled, one value per outcome), following Kazis et al. (Med Care 1989). Pooled
# rather than per-stratum to avoid range restriction from baseline-defined
# strata, and it uses fully-observed baseline BP (robust to follow-up attrition).
# Numerator = adjusted emmean Δ, matching the primary estimate. Signed
# (negative = BP reduction).
#
# Inputs:
#   Analyses/5_1/derived/emmeans_results.rds
#   Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/5_1/derived/results/):
#   5_1_table_2_within_stratum.csv   — per-stratum Δ rows (the primary table)
#   5_1_table_2_cross_stratum.csv    — cross-stratum Δ-difference rows
#   5_1_table_2_paper.md             — formatted, paper-ready markdown

library(tidyverse)

project_dir <- "./PRIMARY"
emm_path    <- file.path(project_dir, "Analyses", "5_1", "derived", "emmeans_results.rds")
uf_path     <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
results_dir <- file.path(project_dir, "Analyses", "5_1", "derived", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

emmeans_results <- read_rds(emm_path)
user_frame      <- read_rds(uf_path)

# --- Cohort baseline SD: denominator for Cohen's d (one value per outcome) ---
# Pooled across strata to avoid range restriction from baseline-defined strata;
# uses fully-observed baseline BP (robust to follow-up attrition). Kazis (1989).
baseline_sd <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female")) %>%
  summarise(SBP = sd(BASELINE_SBP, na.rm = TRUE),
            DBP = sd(BASELINE_DBP, na.rm = TRUE)) %>%
  unlist()
cat("Cohort baseline SD (Cohen's d denominator):\n")
cat("  SBP:", round(baseline_sd[["SBP"]], 4), " DBP:", round(baseline_sd[["DBP"]], 4), "\n\n")

# --- CI columns: detect if present, else derive from estimate + SE + df ---
# emmeans `contrast()` does not include CI columns by default; derive Wald CIs
# from estimate ± qt(.975, df) * SE (z-based when df = Inf).
ensure_ci <- function(dat) {
  nm <- names(dat)
  if (all(c("asymp.LCL", "asymp.UCL") %in% nm)) {
    dat$.lcl <- dat$asymp.LCL
    dat$.ucl <- dat$asymp.UCL
  } else if (all(c("lower.CL", "upper.CL") %in% nm)) {
    dat$.lcl <- dat$lower.CL
    dat$.ucl <- dat$upper.CL
  } else if (all(c("estimate", "SE") %in% nm)) {
    crit <- if ("df" %in% nm) ifelse(is.finite(dat$df), qt(0.975, dat$df), 1.959964) else 1.959964
    dat$.lcl <- dat$estimate - crit * dat$SE
    dat$.ucl <- dat$estimate + crit * dat$SE
  } else {
    stop("Could not find or derive CI columns; available: ", paste(nm, collapse = ", "))
  }
  dat
}

# --- Per-cell paired counts (n contributing a baseline+timepoint pair) ---
# One row per (outcome × stratum × timepoint): n_paired.
paired_n <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female")) %>%
  select(STUDY_ID, BASELINE_AHA,
         BASELINE_SBP, MO3_SBP, MO6_SBP,
         BASELINE_DBP, MO3_DBP, MO6_DBP) %>%
  mutate(
    SBP_MO3 = MO3_SBP - BASELINE_SBP,
    SBP_MO6 = MO6_SBP - BASELINE_SBP,
    DBP_MO3 = MO3_DBP - BASELINE_DBP,
    DBP_MO6 = MO6_DBP - BASELINE_DBP
  ) %>%
  select(STUDY_ID, BASELINE_AHA, SBP_MO3, SBP_MO6, DBP_MO3, DBP_MO6) %>%
  pivot_longer(-c(STUDY_ID, BASELINE_AHA),
               names_to = c("outcome", "timepoint"),
               names_sep = "_",
               values_to = "paired_diff") %>%
  filter(!is.na(paired_diff)) %>%
  group_by(outcome, BASELINE_AHA, timepoint) %>%
  summarise(n_paired = n(), .groups = "drop") %>%
  mutate(stratum = factor(BASELINE_AHA,
                          levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
                          labels = c("Elevated", "Stage 1", "Stage 2"))) %>%
  select(outcome, stratum, timepoint, n_paired)

# --- Pull within-stratum Δ for unadjusted + adjusted ---
extract_within <- function(tib, outcome_label, model_label) {
  tib %>%
    ensure_ci() %>%
    transmute(
      outcome   = outcome_label,
      stratum   = factor(BASELINE_AHA,
                         levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
                         labels = c("Elevated", "Stage 1", "Stage 2")),
      timepoint = sub(" - BASELINE", "", contrast),
      model     = model_label,
      delta     = estimate,
      lcl       = .lcl,
      ucl       = .ucl,
      p         = p.value,
      p_holm    = p.holm
    )
}

within_rows <- bind_rows(
  extract_within(emmeans_results$sbp_unadj$within_stratum, "SBP", "unadj"),
  extract_within(emmeans_results$sbp_adj$within_stratum,   "SBP", "adj"),
  extract_within(emmeans_results$dbp_unadj$within_stratum, "DBP", "unadj"),
  extract_within(emmeans_results$dbp_adj$within_stratum,   "DBP", "adj")
)

# Reshape to wide: one row per (outcome × stratum × timepoint), unadj/adj side-by-side
within_wide <- within_rows %>%
  pivot_wider(names_from = model,
              values_from = c(delta, lcl, ucl, p, p_holm),
              names_glue  = "{model}_{.value}") %>%
  left_join(paired_n, by = c("outcome", "stratum", "timepoint")) %>%
  mutate(
    d_baseline = adj_delta / unname(baseline_sd[outcome])   # adjusted Δ / pooled baseline SD (signed)
  ) %>%
  arrange(factor(outcome, levels = c("SBP", "DBP")), stratum, timepoint) %>%
  select(outcome, stratum, timepoint, n_paired,
         unadj_delta, unadj_lcl, unadj_ucl, unadj_p, unadj_p_holm,
         adj_delta,   adj_lcl,   adj_ucl,   adj_p,   adj_p_holm,
         d_baseline)

write_csv(within_wide, file.path(results_dir, "5_1_table_2_within_stratum.csv"))

# --- Cross-stratum Δ-difference contrasts ---
# emmeans `contrast(delta_emm, "pairwise", by = "contrast")` produces a tibble
# with columns: `contrast1` (stratum pair, e.g. "Elevated - Stage 1 Hypertension"),
# `contrast` (timepoint, e.g. "MO3 - BASELINE"), estimate, SE, df, z.ratio, p.value.
# The original by-variable becomes `contrast`; the new pairwise contrast becomes
# `contrast1` (because emmeans re-uses 'contrast' and disambiguates with a suffix).
extract_cross <- function(tib, outcome_label, model_label) {
  tib <- ensure_ci(tib)
  reserved <- c("contrast", "estimate", "SE", "df", "z.ratio", "t.ratio",
                "p.value", "asymp.LCL", "asymp.UCL", "lower.CL", "upper.CL",
                ".lcl", ".ucl")
  pair_col <- setdiff(names(tib), reserved)
  if (length(pair_col) == 0) {
    cat("DEBUG cross_stratum cols → ", paste(names(tib), collapse = ", "), "\n")
    stop("Could not identify pair column in cross_stratum output")
  }
  tibble(
    outcome   = outcome_label,
    timepoint = sub(" - BASELINE", "", as.character(tib$contrast)),
    pair      = as.character(tib[[ pair_col[1] ]]),
    model     = model_label,
    diff      = tib$estimate,
    lcl       = tib$.lcl,
    ucl       = tib$.ucl,
    p         = tib$p.value,
    p_holm    = tib$p.holm
  )
}

cross_rows <- bind_rows(
  extract_cross(emmeans_results$sbp_unadj$cross_stratum, "SBP", "unadj"),
  extract_cross(emmeans_results$sbp_adj$cross_stratum,   "SBP", "adj"),
  extract_cross(emmeans_results$dbp_unadj$cross_stratum, "DBP", "unadj"),
  extract_cross(emmeans_results$dbp_adj$cross_stratum,   "DBP", "adj")
)

cross_wide <- cross_rows %>%
  pivot_wider(names_from = model,
              values_from = c(diff, lcl, ucl, p, p_holm),
              names_glue  = "{model}_{.value}") %>%
  arrange(factor(outcome, levels = c("SBP", "DBP")), timepoint, pair) %>%
  select(outcome, timepoint, pair,
         unadj_diff, unadj_lcl, unadj_ucl, unadj_p, unadj_p_holm,
         adj_diff,   adj_lcl,   adj_ucl,   adj_p,   adj_p_holm)

write_csv(cross_wide, file.path(results_dir, "5_1_table_2_cross_stratum.csv"))

# --- Paper-ready markdown rendering ---
fmt_ci <- function(est, lcl, ucl, digits = 2) {
  sprintf("%.*f (%.*f, %.*f)", digits, est, digits, lcl, digits, ucl)
}
fmt_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}

within_md <- within_wide %>%
  mutate(
    `Δ unadj (95% CI)` = fmt_ci(unadj_delta, unadj_lcl, unadj_ucl),
    `Δ adj (95% CI)`   = fmt_ci(adj_delta,   adj_lcl,   adj_ucl),
    `p (adj, Holm)`    = fmt_p(adj_p_holm),
    `Cohen's d`        = sprintf("%.3f", d_baseline),
    n_paired           = format(n_paired, big.mark = ",")
  ) %>%
  select(Outcome = outcome, Stratum = stratum, Timepoint = timepoint,
         `n (paired)` = n_paired,
         `Δ unadj (95% CI)`, `Δ adj (95% CI)`, `p (adj, Holm)`, `Cohen's d`)

cross_md <- cross_wide %>%
  mutate(
    `Δ-diff unadj (95% CI)` = fmt_ci(unadj_diff, unadj_lcl, unadj_ucl),
    `Δ-diff adj (95% CI)`   = fmt_ci(adj_diff,   adj_lcl,   adj_ucl),
    `p (adj, Holm)`         = fmt_p(adj_p_holm)
  ) %>%
  select(Outcome = outcome, Timepoint = timepoint, `Stratum contrast` = pair,
         `Δ-diff unadj (95% CI)`, `Δ-diff adj (95% CI)`, `p (adj, Holm)`)

md_lines <- c(
  "# Table 2 — Per-stratum BP change and cross-stratum Δ-difference (SAP §5.1)",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "## Per-stratum Δ from baseline (mmHg)",
  "",
  knitr::kable(within_md, format = "markdown", align = "lllrlllr"),
  "",
  "## Cross-stratum Δ-difference (mmHg)",
  "",
  knitr::kable(cross_md, format = "markdown", align = "lllllr"),
  "",
  "### Notes",
  "",
  "- Unadjusted model: `BP ~ Time * AHA_baseline_stratum + (1 | STUDY_ID)`",
  "- Adjusted model adds: `age + sex + baseline_BMI + baseline_HTN_status + cyclic time_of_day + season + Time:age + Time:sex + Time:baseline_BMI + Time:baseline_HTN_status`",
  "- AHA reference level: Elevated. Time reference: BASELINE.",
  "- 95% CIs are asymptotic (z-based) Wald intervals from emmeans.",
  "- p-values shown are from the **adjusted** LMM (the fuller spec incl. baseline-HTN), **Holm-adjusted within contrast family**. The primary within-stratum Δ is led by its 95% CI (estimation-first); raw and unadjusted-model p-values are in the CSV exports.",
  "- Cohen's d = (adjusted LMM Δ) / (pooled baseline SD of the outcome, all strata combined; Kazis et al., Med Care 1989). Signed (negative = BP reduction). Harmonized with §5.5 Table S3. Supplementary effect-size summary."
)
writeLines(md_lines, file.path(results_dir, "5_1_table_2_paper.md"))

cat("Wrote:\n")
cat("  ", file.path(results_dir, "5_1_table_2_within_stratum.csv"), "\n")
cat("  ", file.path(results_dir, "5_1_table_2_cross_stratum.csv"), "\n")
cat("  ", file.path(results_dir, "5_1_table_2_paper.md"), "\n")
cat("\n--- Per-stratum Δ table preview ---\n")
print(within_wide, n = Inf, width = Inf)
cat("\n--- Cross-stratum Δ-difference preview ---\n")
print(cross_wide, n = Inf, width = Inf)
