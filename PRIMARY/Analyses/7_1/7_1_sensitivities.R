# 7_1_sensitivities.R
# SAP §7.1 (v0.64+) sensitivity analysis: HTN-medication-naive subset.
#
# Single sensitivity (HTN-medication-naive) refit against BOTH §7.1 primary
# responder definitions: SBP-responder and DBP-responder. This is the
# "primary × sensitivity" cross-product that the dual-primary v0.64 framing
# requires (DBP-responder is no longer a sensitivity — it is co-primary).
#
# Cohort:
#   HTN-naive × SBP: §2.2 hypertensive + sex restriction + non-missing
#                    SBP_RESPONDER_6MO + HTN_MED_BASELINE == "never_logged"
#   HTN-naive × DBP: §2.2 hypertensive + sex restriction + non-missing
#                    DBP_RESPONDER_6MO + HTN_MED_BASELINE == "never_logged"
#
# Fits: 1 sensitivity × 2 responder_defs × 5 metrics = 10 LMMs total.
#
# Outputs (Analyses/7_1/derived/results/), all carrying `responder_def`:
#   7_1_sensitivity_anova.csv           — Type III ANOVA per (responder_def × metric)
#   7_1_sensitivity_emmeans.csv         — emmeans per cell
#   7_1_sensitivity_delta_contrasts.csv — per-group Δ from baseline
#   7_1_sensitivity_comparison.csv      — sensitivity-vs-primary p-value table
#   7_1_sensitivity_direction_check.csv — favorable-direction concordance per cell
#   7_1_sensitivities_paper.md          — paper-ready markdown

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")  # tidy_contrast(), emit_sizes()

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "7_1", "derived")
results_dir <- file.path(out_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "7_1_sensitivities.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)
cat("§7.1 sensitivities — HTN-medication-naive × {SBP, DBP} responder defs —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read frame, derive DBP-responder indicator inline ---
user_frame <- read_rds(input_path) %>%
  mutate(DBP_RESPONDER_6MO = as.integer((MO6_DBP - BASELINE_DBP) <= -5))
# NOTE: HTN_BASELINE_ON is NOT derived here — the HTN-naive subset is entirely
# never_logged, so HTN_BASELINE_ON would be a constant (not_on_at_baseline) for
# all rows. Including a constant covariate in lmer causes a rank-deficiency
# warning and the estimate is unidentifiable. Omit per §7.1 sensitivity spec.

# --- Five metrics ---
metrics <- c(
  RHR              = "RHR",
  SleepDuration    = "SLEEP_DURATION",
  SleepConsistency = "SLEEP_CONSISTENCY",
  Strain           = "STRAIN",
  HRV              = "HRV"
)

# Favorable direction at MO6 for responder group:
#   "lower"  = responders' MO6 Δ should be more negative than non-responders'
#   "higher" = responders' MO6 Δ should be more positive than non-responders'
metric_favorable_direction <- c(
  RHR              = "lower",
  SleepDuration    = "higher",
  SleepConsistency = "higher",
  Strain           = "higher",
  HRV              = "higher"
)

# --- Cohort builder: HTN-naive subset, parameterized by responder_def ---
build_htn_naive_cohort <- function(d, responder_def) {
  resp_col <- paste0(responder_def, "_RESPONDER_6MO")
  d %>%
    filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
           SEX %in% c("male", "female"),
           !is.na(.data[[resp_col]]),
           HTN_MED_BASELINE == "never_logged") %>%
    mutate(responder_group = factor(
      ifelse(.data[[resp_col]] == 1L, "responder", "non_responder"),
      levels = c("non_responder", "responder")
    ))
}

# --- Helper: build long-format dataframe for one metric ---
build_long <- function(metric_prefix, d) {
  cols <- paste0(metric_prefix, c("_BASELINE_MEAN", "_MO3_MEAN", "_MO6_MEAN"))
  d %>%
    select(STUDY_ID, responder_group, AGE_AT_START, SEX, BASELINE_BMI,
           all_of(cols)) %>%
    rename(
      baseline = !!sym(cols[1]),
      `3mo`    = !!sym(cols[2]),
      `6mo`    = !!sym(cols[3])
    ) %>%
    pivot_longer(c(baseline, `3mo`, `6mo`),
                 names_to = "Time", values_to = "value") %>%
    mutate(
      Time = factor(Time, levels = c("baseline", "3mo", "6mo")),
      SEX  = factor(SEX)
    ) %>%
    filter(!is.na(value), !is.na(AGE_AT_START), !is.na(BASELINE_BMI))
}

# --- Helper: fit one metric's LMM ---
fit_metric_lmm <- function(metric_name, metric_prefix, d) {
  long_d <- build_long(metric_prefix, d)
  cat(sprintf("\n  %s (%s): %d rows, %d users\n",
              metric_name, metric_prefix,
              nrow(long_d), n_distinct(long_d$STUDY_ID)))
  f <- as.formula(
    "value ~ Time * responder_group + AGE_AT_START + SEX + BASELINE_BMI + (1 | STUDY_ID)"
  )
  fit <- lmer(f, data = long_d, REML = TRUE)
  list(fit = fit, data = long_d, metric = metric_name)
}

# --- Build sensitivity cohorts: HTN-naive × {SBP, DBP} ---
responder_defs <- c("SBP", "DBP")
sensitivity_cohorts <- list()
cohort_sizes <- tibble(responder_def = character(),
                       n_users = integer(),
                       n_responder = integer(),
                       n_non_responder = integer())

cat("Sensitivity cohorts (HTN-medication-naive):\n")
for (rd in responder_defs) {
  scoh <- build_htn_naive_cohort(user_frame, rd)
  sensitivity_cohorts[[rd]] <- scoh
  resp_counts <- scoh %>% count(responder_group)
  cohort_sizes <- bind_rows(cohort_sizes, tibble(
    responder_def    = rd,
    n_users          = nrow(scoh),
    n_responder      = resp_counts$n[resp_counts$responder_group == "responder"],
    n_non_responder  = resp_counts$n[resp_counts$responder_group == "non_responder"]
  ))
  cat(sprintf("  HTN-naive × %s: n = %d users\n", rd, nrow(scoh)))
  print(resp_counts)
}
write_csv(cohort_sizes, file.path(results_dir, "7_1_sensitivity_cohort_sizes.csv"))

# --- Fit: HTN-naive × {SBP, DBP} × 5 metrics = 10 LMMs ---
all_anova   <- list()
all_emm     <- list()
all_delta   <- list()

for (rd in responder_defs) {
  scoh <- sensitivity_cohorts[[rd]]
  cat("\n\n========================================\n")
  cat(sprintf("Sensitivity: HTN_naive × %s-responder\n", rd))
  cat("========================================\n")

  for (m in names(metrics)) {
    fr <- fit_metric_lmm(m, metrics[[m]], scoh)
    a <- suppressMessages(anova(fr$fit, type = 3))
    anova_tbl <- as_tibble(a, rownames = "term") %>%
      mutate(sensitivity = "HTN_naive", responder_def = rd, metric = m)
    all_anova[[paste0(rd, "__", m)]] <- anova_tbl

    emm <- suppressMessages(emmeans(fr$fit, ~ Time | responder_group))
    delta <- suppressMessages(contrast(emm, method = "trt.vs.ctrl", adjust = "none"))
    all_emm[[paste0(rd, "__", m)]]   <- as_tibble(as.data.frame(emm)) %>%
      mutate(sensitivity = "HTN_naive", responder_def = rd, metric = m)
    all_delta[[paste0(rd, "__", m)]] <- tidy_contrast(delta) %>%  # adds p.holm within-family per v0.84
      mutate(sensitivity = "HTN_naive", responder_def = rd, metric = m)
  }
}

# --- Combine + write CSVs ---
anova_df <- bind_rows(all_anova) %>%
  select(sensitivity, responder_def, metric, term, everything())
write_csv(anova_df, file.path(results_dir, "7_1_sensitivity_anova.csv"))

emm_df <- bind_rows(all_emm) %>%
  select(sensitivity, responder_def, metric, everything())
write_csv(emm_df, file.path(results_dir, "7_1_sensitivity_emmeans.csv"))

delta_df <- bind_rows(all_delta) %>%
  select(sensitivity, responder_def, metric, everything())
write_csv(delta_df, file.path(results_dir, "7_1_sensitivity_delta_contrasts.csv"))

cat("\n\n--- Interaction p-values per (responder_def × metric) ---\n")
print(anova_df %>%
        filter(term == "Time:responder_group") %>%
        select(sensitivity, responder_def, metric, `F value`, `Pr(>F)`),
      n = Inf, width = Inf)

# --- Compare to primary §7.1 (matched on responder_def) ---
primary_anova <- read_csv(file.path(results_dir, "7_1_anova.csv"),
                          show_col_types = FALSE) %>%
  filter(term == "Time:responder_group") %>%
  transmute(responder_def, metric,
            primary_F = `F value`, primary_p = `Pr(>F)`,
            primary_sig = `Pr(>F)` < 0.005)

sens_interaction <- anova_df %>%
  filter(term == "Time:responder_group") %>%
  transmute(sensitivity, responder_def, metric,
            sens_F = `F value`, sens_p = `Pr(>F)`,
            sens_sig = `Pr(>F)` < 0.005)

comparison <- sens_interaction %>%
  left_join(primary_anova, by = c("responder_def", "metric")) %>%
  mutate(
    p_direction = ifelse(sens_p < primary_p, "tightens", "loosens"),
    sig_concordant = sens_sig == primary_sig
  ) %>%
  arrange(responder_def, metric)

write_csv(comparison, file.path(results_dir, "7_1_sensitivity_comparison.csv"))

cat("\n--- Comparison: HTN-naive vs §7.1 primary, by responder_def × metric ---\n")
print(comparison, n = Inf, width = Inf)

# --- Directionality check per (responder_def × metric) ---
direction_check <- delta_df %>%
  filter(contrast == "6mo - baseline") %>%
  select(sensitivity, responder_def, metric, responder_group, estimate) %>%
  pivot_wider(names_from = responder_group, values_from = estimate) %>%
  mutate(
    resp_minus_nonresp = responder - non_responder,
    favorable_direction = metric_favorable_direction[metric],
    favorable_concordant = case_when(
      favorable_direction == "lower"  ~ resp_minus_nonresp < 0,
      favorable_direction == "higher" ~ resp_minus_nonresp > 0,
      TRUE                            ~ NA
    )
  ) %>%
  select(sensitivity, responder_def, metric,
         non_responder, responder, resp_minus_nonresp,
         favorable_direction, favorable_concordant) %>%
  arrange(responder_def, metric)

write_csv(direction_check, file.path(results_dir, "7_1_sensitivity_direction_check.csv"))

cat("\n--- Directionality check at MO6 (responder − non_responder, vs favorable direction) ---\n")
print(direction_check, n = Inf, width = Inf)

direction_summary <- direction_check %>%
  group_by(responder_def) %>%
  summarise(n_metrics   = n(),
            n_favorable = sum(favorable_concordant),
            .groups     = "drop")
cat("\n--- Direction concordance summary by responder_def ---\n")
print(direction_summary)

# --- Convergence summary ---
convergence <- comparison %>%
  group_by(responder_def) %>%
  summarise(
    n_metrics              = n(),
    primary_sig_metrics    = sum(primary_sig),
    sens_sig_metrics       = sum(sens_sig),
    sig_concordant_metrics = sum(sig_concordant),
    .groups = "drop"
  )
cat("\n--- Convergence summary (HTN-naive vs primary) ---\n")
print(convergence)

# --- Emit standardized sizes for Table A.2 (sensitivity; per SAP review v0.84) ---
# Ceiling = primary cohort n for the matching responder_def (from primary CSVs).
primary_cohort_sizes <- tryCatch(
  readr::read_csv(file.path(results_dir, "7_1_cohort_sizes.csv"), show_col_types = FALSE),
  error = function(e) NULL
)
if (!is.null(primary_cohort_sizes)) {
  sizes_71s <- bind_rows(lapply(responder_defs, function(rd) {
    scoh <- sensitivity_cohorts[[rd]]
    ceiling_n <- primary_cohort_sizes$n_users[primary_cohort_sizes$responder_def == rd]
    bind_rows(lapply(names(metrics), function(m) {
      tibble(
        section       = "7.1",
        role          = "sensitivity",
        analysis      = "WHOOP trajectory LMM — HTN-naive",
        cohort        = "htn_naive",
        outcome       = m,
        predictor_set = rd,
        n_users       = as.integer(nrow(scoh)),   # before metric-specific NA drop
        n_obs         = NA_integer_,              # metric rows vary; use CSV if needed
        n_events      = NA_integer_,
        ceiling_n     = as.integer(ceiling_n)
      )
    }))
  }))
  emit_sizes(sizes_71s, "7_1_sensitivity", results_dir)
} else {
  message("emit_sizes skipped for sensitivity: primary 7_1_cohort_sizes.csv not found — run 7_1_whoop_trajectories.R first.")
}

# =============================================================================
# Paper-ready markdown
# =============================================================================
fmt_p <- function(p) {
  ifelse(p < 1e-10, formatC(p, format = "e", digits = 1),
  ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                    sprintf("%.3f", p)))
}

md_lines <- c(
  "# §7.1 sensitivities — HTN-medication-naive × {SBP, DBP} responder",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "Per SAP v0.64, §7.1 has two co-primary responder definitions (SBP-responder, DBP-responder). The HTN-medication-naive sensitivity is refit against both primary frames; this file reports both cross-product cells against the matching primary result. Bold = significant at α = 0.005.",
  "",
  "## Cohort sizes (HTN-naive × responder_def)",
  ""
)
for (i in seq_len(nrow(cohort_sizes))) {
  r <- cohort_sizes[i, ]
  md_lines <- c(md_lines, sprintf(
    "- **HTN-naive × %s**: n = %s users (responders: %s; non-responders: %s)",
    r$responder_def,
    format(r$n_users, big.mark = ","),
    format(r$n_responder, big.mark = ","),
    format(r$n_non_responder, big.mark = ",")
  ))
}

md_lines <- c(md_lines, "",
  "## Convergence summary (HTN-naive vs primary, per responder_def)",
  ""
)
for (i in seq_len(nrow(convergence))) {
  r <- convergence[i, ]
  md_lines <- c(md_lines, sprintf(
    "- **%s-responder**: HTN-naive sig metrics %d/%d (primary sig: %d/%d); significance concordant in %d/%d metrics.",
    r$responder_def,
    r$sens_sig_metrics, r$n_metrics,
    r$primary_sig_metrics, r$n_metrics,
    r$sig_concordant_metrics, r$n_metrics
  ))
}

md_lines <- c(md_lines, "",
  "## Cell-by-cell comparison (Time × responder_group p-values)",
  ""
)
for (rd in responder_defs) {
  md_lines <- c(md_lines, paste0("### ", rd, "-responder frame"))
  for (m in names(metrics)) {
    r <- comparison %>% filter(responder_def == rd, metric == m)
    if (nrow(r) == 0) next
    primary_bold <- ifelse(r$primary_sig, "**", "")
    sens_bold    <- ifelse(r$sens_sig, "**", "")
    md_lines <- c(md_lines, sprintf(
      "- **%s**: primary F = %.1f, %sp = %s%s   |   HTN-naive F = %.1f, %sp = %s%s",
      r$metric,
      r$primary_F, primary_bold, fmt_p(r$primary_p), primary_bold,
      r$sens_F,    sens_bold,    fmt_p(r$sens_p),    sens_bold
    ))
  }
  md_lines <- c(md_lines, "")
}

md_lines <- c(md_lines,
  "## Directionality at MO6 (responder − non_responder vs favorable direction)",
  ""
)
for (rd in responder_defs) {
  md_lines <- c(md_lines, paste0("### ", rd, "-responder frame"))
  dr <- direction_check %>% filter(responder_def == rd) %>% arrange(metric)
  for (i in seq_len(nrow(dr))) {
    r <- dr[i, ]
    mark <- ifelse(isTRUE(r$favorable_concordant), "✓", "✗")
    md_lines <- c(md_lines, sprintf(
      "- %s: responder − non_responder = %+0.3f, favorable = %s %s",
      r$metric, r$resp_minus_nonresp, r$favorable_direction, mark
    ))
  }
  md_lines <- c(md_lines, sprintf(
    "  → %d/%d metrics direction-concordant.",
    sum(dr$favorable_concordant), nrow(dr)
  ), "")
}

md_lines <- c(md_lines,
  "## Notes",
  "",
  "- **HTN-medication-naive subset** restricts to §2.2 users with `never_logged` HTN-medication status. Tests whether the trajectory divergences in §7.1 primary are driven by users on antihypertensive medications.",
  "- **Cross-product against both primary frames**: per SAP v0.64, the sensitivity is run against SBP-responder and DBP-responder primary definitions separately. Both cross-product cells are reported above to support the dual-primary framing without selection bias.",
  "- **Convergence criteria** (per SAP §7.1): rough magnitude on the interaction estimate and significance concordance with the matched-responder_def primary cell. A metric that flips significance class flags possible composition dependence.",
  "- **DBP-responder is NOT in this file as a separate sensitivity** — DBP-responder grouping is co-primary per v0.64 and is reported in `7_1_paper_table4.md`. The previous DBP_responder sensitivity cell has been folded into the primary table.",
  "- Pairs with: §7.1 primary (`7_1_whoop_trajectories.R`, `7_1_paper_table4.md`)."
)

writeLines(md_lines, file.path(results_dir, "7_1_sensitivities_paper.md"))
cat("\nWrote markdown:\n  ", file.path(results_dir, "7_1_sensitivities_paper.md"), "\n")

cat("\n=== §7.1 sensitivities done ===\n")
