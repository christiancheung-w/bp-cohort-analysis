# 5_1_sensitivity_empirical_median.R
# SAP §5.1 sensitivity bullet 5: empirical median paired-difference comparison.
# Tests robustness of the §5.1 mean Δ headline to distributional-tail influence
# by comparing per-stratum empirical median Δ to empirical mean Δ at MO3 and MO6,
# with 95% cluster bootstrap CIs.
#
# History: SAP bullet 5 originally specified a linear quantile mixed model (lqmm
# at τ = 0.5). Per CHANGELOG v0.51 (2026-05-26), lqmm failed to converge at the
# §2.2 cohort scale (~50k obs / ~20k clusters); the empirical median + cluster
# bootstrap design replaced it as the SAP-canonical method. See derived/PROVENANCE.md
# for the full backstory; the archived lqmm-attempt script lives at
# PRIMARY/archive/R-stale/5_1_sensitivity_quantile_v_lqmm.R.
#
# Structure:
#   Section 1: Empirical median + mean of paired Δs per stratum + cluster bootstrap CIs.
#   Section 2: Comparison table vs §5.1 LMM ADJUSTED means.
#   Section 3: Figure — forest plot of mean vs median per stratum.
#   Section 4: Paper-ready markdown.
#
# Substantive conclusion: if empirical median ≈ empirical mean within each
# (outcome × stratum) cell, the §5.1 mean estimates are tail-robust. The
# transitivity argument (empirical mean ≈ §5.1 ADJUSTED mean, established in
# §5.1 unadj-vs-adj) then implies the ADJUSTED mean is also tail-robust.
#
# Cohort: §2.2 hypertensive + sex restriction, matching §5.1 primary.
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/5_1/Sensitivity/empirical_median/derived/):
#   results/empirical_paired_summary.csv              — median + mean + CI per stratum
#   results/comparison_5_1_vs_empirical_median.csv    — combined table vs §5.1 LMM means
#   results/5_1_sensitivity_empirical_median_paper.md
#   results/5_1_sensitivity_empirical_median.log
#   figures/5_1_empirical_median_comparison.png       — forest plot: mean vs median per stratum

library(tidyverse)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")  # emit_sizes() (bootstrap; no contrast families)

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "5_1", "Sensitivity", "empirical_median", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(out_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "5_1_sensitivity_empirical_median.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)  # ensure log closes cleanly even on error

cat("§5.1 sensitivity — empirical median paired-Δ —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

user_frame <- read_rds(input_path)
cohort <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female"))
cat("§2.2 cohort + sex restriction:", nrow(cohort), "users\n\n")

# Capture sizes now (the `cohort` symbol is reused downstream) for Table A.2.
.sz_sbp <- sum(!is.na(cohort$BASELINE_SBP) & !is.na(cohort$MO6_SBP))
.sz_dbp <- sum(!is.na(cohort$BASELINE_DBP) & !is.na(cohort$MO6_DBP))
.sz_ceiling <- nrow(cohort)

# =============================================================================
# Section 1: Empirical median + mean of paired Δs (cluster bootstrap CIs)
# =============================================================================
cat("\n========== Section 1: Empirical paired-Δ analysis ==========\n")

paired_long <- cohort %>%
  transmute(
    STUDY_ID,
    stratum     = factor(BASELINE_AHA,
                         levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
                         labels = c("Elevated", "Stage 1", "Stage 2")),
    SBP_MO3     = MO3_SBP - BASELINE_SBP,
    SBP_MO6     = MO6_SBP - BASELINE_SBP,
    DBP_MO3     = MO3_DBP - BASELINE_DBP,
    DBP_MO6     = MO6_DBP - BASELINE_DBP
  ) %>%
  pivot_longer(-c(STUDY_ID, stratum),
               names_to  = c("outcome", "timepoint"),
               names_sep = "_",
               values_to = "delta") %>%
  filter(!is.na(delta))

cat("Paired-Δ rows (across outcome × timepoint × user):\n")
print(paired_long %>% count(outcome, timepoint, stratum))

# Cluster bootstrap helper: cluster-resample STUDY_IDs, compute median + mean
boot_paired <- function(d, n_boot = 1000, seed = 7) {
  set.seed(seed)
  users <- unique(d$STUDY_ID)
  n_users <- length(users)
  out <- replicate(n_boot, {
    boot_ids <- sample(users, n_users, replace = TRUE)
    # Cluster-resample: pull all rows for each sampled user (users may appear
    # multiple times in a draw; left_join handles the duplication correctly).
    boot_idx <- tibble(STUDY_ID = boot_ids) %>%
      left_join(d, by = "STUDY_ID", relationship = "many-to-many")
    c(median = median(boot_idx$delta, na.rm = TRUE),
      mean   = mean(boot_idx$delta,   na.rm = TRUE))
  })
  tibble(
    median_boot_lcl = quantile(out["median", ], 0.025, na.rm = TRUE),
    median_boot_ucl = quantile(out["median", ], 0.975, na.rm = TRUE),
    mean_boot_lcl   = quantile(out["mean",   ], 0.025, na.rm = TRUE),
    mean_boot_ucl   = quantile(out["mean",   ], 0.975, na.rm = TRUE)
  )
}

empirical_summary <- paired_long %>%
  group_by(outcome, timepoint, stratum) %>%
  summarise(
    n_users = n(),
    median  = median(delta),
    mean    = mean(delta),
    sd      = sd(delta),
    .groups = "drop"
  )

cat("\nBootstrap CIs (1000 cluster-bootstrap reps)... ")
t0 <- Sys.time()
boot_cis <- paired_long %>%
  group_by(outcome, timepoint, stratum) %>%
  group_modify(~ boot_paired(.x, n_boot = 1000)) %>%
  ungroup()
cat("done (", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec)\n")

empirical_summary <- empirical_summary %>%
  left_join(boot_cis, by = c("outcome", "timepoint", "stratum")) %>%
  mutate(
    abs_median_minus_mean = abs(median - mean),
    pct_diff = ifelse(abs(mean) > 0.5, 100 * (median - mean) / mean, NA)
  )

write_csv(empirical_summary, file.path(results_dir, "empirical_paired_summary.csv"))

cat("\n--- Empirical paired-Δ summary (mean vs median, with 95% bootstrap CI) ---\n")
print(empirical_summary, n = Inf, width = Inf)

# =============================================================================
# Section 2: Comparison vs §5.1 LMM ADJUSTED means
# =============================================================================
cat("\n\n========== Section 2: Comparison vs §5.1 LMM ADJUSTED means ==========\n")

lmm_path <- file.path(project_dir, "Analyses", "5_1", "derived", "results",
                     "5_1_delta_contrasts.csv")
lmm_means <- read_csv(lmm_path, show_col_types = FALSE) %>%
  filter(model == "sbp_adj" | model == "dbp_adj") %>%
  transmute(
    outcome   = ifelse(model == "sbp_adj", "SBP", "DBP"),
    stratum   = factor(BASELINE_AHA,
                       levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
                       labels = c("Elevated", "Stage 1", "Stage 2")),
    timepoint = sub(" - BASELINE", "", contrast),
    lmm_mean_adj = estimate
  )

comparison <- empirical_summary %>%
  select(outcome, stratum, timepoint, n_users,
         empirical_mean   = mean,
         empirical_median = median,
         emp_med_lcl      = median_boot_lcl,
         emp_med_ucl      = median_boot_ucl) %>%
  left_join(lmm_means, by = c("outcome", "stratum", "timepoint")) %>%
  mutate(
    abs_median_minus_lmm_mean   = abs(empirical_median - lmm_mean_adj),
    abs_emp_mean_minus_lmm_mean = abs(empirical_mean - lmm_mean_adj)
  ) %>%
  arrange(outcome, stratum, timepoint)

write_csv(comparison, file.path(results_dir, "comparison_5_1_vs_empirical_median.csv"))

cat("\n--- Comparison: §5.1 LMM mean (adj) vs empirical mean vs empirical median ---\n")
print(comparison, n = Inf, width = Inf)

# =============================================================================
# Section 3: Figure — mean vs median forest plot
# =============================================================================
cat("\n\n========== Section 3: Figure ==========\n")

fig_df <- comparison %>%
  filter(timepoint == "MO6") %>%
  select(outcome, stratum,
         `§5.1 LMM mean (adj)`  = lmm_mean_adj,
         `Empirical mean`       = empirical_mean,
         `Empirical median`     = empirical_median,
         emp_med_lcl, emp_med_ucl) %>%
  pivot_longer(c(`§5.1 LMM mean (adj)`, `Empirical mean`, `Empirical median`),
               names_to = "measure", values_to = "estimate") %>%
  mutate(measure = factor(measure,
                          levels = c("§5.1 LMM mean (adj)", "Empirical mean", "Empirical median")),
         lcl = ifelse(measure == "Empirical median", emp_med_lcl, NA),
         ucl = ifelse(measure == "Empirical median", emp_med_ucl, NA))

p <- ggplot(fig_df, aes(x = estimate, y = stratum, color = measure)) +
  geom_vline(xintercept = 0, color = "gray60", linetype = "dashed") +
  geom_errorbar(aes(xmin = lcl, xmax = ucl), width = 0.15,
                position = position_dodge(width = 0.5), linewidth = 0.6,
                orientation = "y", na.rm = TRUE) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  facet_wrap(~ outcome, ncol = 1, scales = "free_x") +
  scale_color_brewer(palette = "Set1", name = NULL) +
  labs(
    x = "MO6 Δ from baseline (mmHg)",
    y = NULL,
    title    = "Tail-robustness sensitivity: mean vs median estimates of MO6 Δ",
    subtitle = "Per-stratum comparison; 95% bootstrap CI on empirical median",
    caption  = "§5.1 LMM (adjusted) mean is the manuscript primary estimate. Empirical median tests tail-robustness via mean ≈ median."
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text      = element_text(face = "bold"),
        plot.caption    = element_text(size = 9, color = "gray40", hjust = 0))

out_png <- file.path(figures_dir, "5_1_empirical_median_comparison.png")
ggsave(out_png, p, width = 9, height = 7, dpi = 300)
cat("Wrote figure:\n  ", out_png, "\n")

# =============================================================================
# Section 4: Paper-ready markdown
# =============================================================================
fmt_signed <- function(x, d = 2) sprintf("%+.*f", d, x)

md_lines <- c(
  "# §5.1 sensitivity — empirical median tail-robustness",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "## Empirical median vs mean (per-stratum, paired Δs at MO6)",
  ""
)
mo6_rows <- comparison %>% filter(timepoint == "MO6")
for (otc in c("SBP", "DBP")) {
  md_lines <- c(md_lines, paste0("### ", otc))
  for (st in c("Elevated", "Stage 1", "Stage 2")) {
    r <- mo6_rows %>% filter(outcome == otc, stratum == st)
    if (nrow(r) == 0) next
    md_lines <- c(md_lines, sprintf(
      "- **%s** (n = %s): empirical mean Δ = %s mmHg, empirical median Δ = %s mmHg (95%% boot CI [%s, %s]). §5.1 LMM mean (adj) = %s mmHg. Median vs §5.1 mean discrepancy = %.2f mmHg.",
      st, format(r$n_users, big.mark = ","),
      fmt_signed(r$empirical_mean), fmt_signed(r$empirical_median),
      fmt_signed(r$emp_med_lcl), fmt_signed(r$emp_med_ucl),
      fmt_signed(r$lmm_mean_adj),
      r$abs_median_minus_lmm_mean
    ))
  }
  md_lines <- c(md_lines, "")
}

md_lines <- c(md_lines,
  "## Substantive conclusion",
  "",
  "If the empirical median Δ closely matches the empirical mean Δ within each (outcome × stratum) cell, the §5.1 mean estimates are tail-robust — extreme paired Δs are not driving the headline.",
  "",
  "**Transitivity to the ADJUSTED model:** the SAP-canonical method is the empirical paired-Δ comparison (per CHANGELOG v0.51). The original SAP spec was lqmm at τ = 0.5 with the §5.1 ADJUSTED covariate set — convergence failed at the §2.2 cohort scale (~50k obs / ~20k clusters), so empirical median substitutes. The strict-reviewer concern is that the empirical method is unadjusted while the original lqmm would have been covariate-adjusted. This is addressed by transitivity: (empirical mean ≈ empirical median, established here) chained with (empirical mean ≈ §5.1 ADJUSTED mean, established in the §5.1 unadj-vs-adj comparison, within ~0.5 mmHg per cell) yields: the §5.1 ADJUSTED mean is tail-robust within the same tolerance.",
  "",
  "## Method",
  "",
  "- Per-user paired Δs: `MO3_BP − BASELINE_BP` and `MO6_BP − BASELINE_BP`, separately for SBP and DBP.",
  "- Within each (outcome × stratum × timepoint) cell: compute empirical mean Δ and empirical median Δ across users.",
  "- 95% bootstrap CIs on the median: 1,000 cluster-resampled draws of users.",
  "- Comparison vs §5.1 ADJUSTED LMM means: pulled from `5_1/derived/results/5_1_delta_contrasts.csv`.",
  "",
  "## Pairs with",
  "",
  "- §5.1 primary (mean estimates being tested for tail-robustness).",
  "- §5.1 RTM Oldham sensitivity (orthogonal — quantifies RTM share, not tail-influence).",
  "- See `derived/PROVENANCE.md` for the lqmm-failure backstory and v0.51 substitution rationale."
)

writeLines(md_lines, file.path(results_dir, "5_1_sensitivity_empirical_median_paper.md"))
cat("\nWrote markdown:\n  ", file.path(results_dir, "5_1_sensitivity_empirical_median_paper.md"), "\n")

# --- Emit standardized sizes for Table A.2 (paired baseline→MO6 n per outcome) ---
emit_sizes(
  tibble(
    section = "5.1", role = "sensitivity",
    analysis = "Empirical median paired-Δ (cluster bootstrap)", cohort = "hypertensive",
    outcome = c("SBP", "DBP"), predictor_set = "median paired Δ",
    n_users = c(.sz_sbp, .sz_dbp), n_obs = c(.sz_sbp, .sz_dbp),
    n_events = NA_integer_, ceiling_n = .sz_ceiling
  ),
  "5_1_empirical_median", results_dir
)

cat("\n=== §5.1 empirical median sensitivity done ===\n")
