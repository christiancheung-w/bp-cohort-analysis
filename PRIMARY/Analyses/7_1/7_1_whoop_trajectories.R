# 7_1_whoop_trajectories.R
# SAP §7.1 (v0.64+): WHOOP-metric trajectories by responder status, run as a
# DUAL-PRIMARY analysis across SBP-responder and DBP-responder groupings.
#
# For each responder definition, five LMMs are fit (one per WHOOP metric: RHR,
# sleep duration, sleep consistency, strain, HRV). The `Time × responder_group`
# interaction is the formal test of trajectory divergence (one p-value per
# metric × responder_def cell). Total fits: 5 metrics × 2 responder defs = 10.
#
# Outcome (per metric): user-mean computed within each 30-day window:
#   baseline = 30 days prior to STUDY_START_DATE → BASELINE_MEAN column
#   3mo      = 30-day window centered on day 90  → MO3_MEAN column
#   6mo      = 30-day window centered on day 180 → MO6_MEAN column
#
# Fixed effects: Time + responder_group + Time × responder_group
# Covariates:    age, sex, baseline_BMI, baseline_HTN_status (main-effect;
#               added per SAP v0.84; omitted in HTN-naive sensitivity — constant)
# Random:        (1 | STUDY_ID)
# Time coding:   baseline (reference), 3mo, 6mo
# Group coding:  non_responder (reference), responder = ΔBP ≤ −5 mmHg
#                (SBP or DBP per `responder_def`)
#
# Cohort: §2.2 hypertensive + sex restriction + non-missing responder indicator
# for the active responder_def. SBP and DBP cohorts differ only by which
# MO6 measurement was non-missing.
#
# Window coverage: window means come pre-filtered by `05_whoop_window_means.R`
# per the SAP §11 25% rule — any 30-day window with <8 days of WHOOP data has
# its 5 metric means set to NA upstream. §7.1 relies on this build contract;
# NA means are auto-dropped by lme4 during fitting.
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/7_1/derived/), all carrying a `responder_def` column:
#   fitted_models.rds                   — list keyed by responder_def, each a list of 5 fits
#   emmeans_results.rds                 — list keyed by responder_def
#   results/7_1_anova.csv               — Type III ANOVA per (responder_def × metric)
#   results/7_1_emmeans.csv             — per-(responder_def × Time × responder × metric) EMMs
#   results/7_1_delta_contrasts.csv     — per-group Δ from baseline at 3mo, 6mo
#   results/7_1_r2.csv                  — Nakagawa marginal + conditional R²
#   results/7_1_direction_check.csv     — favorable-direction concordance per cell
#   results/7_1_paper_table4.md         — paper-ready Table 4 (SBP + DBP sections)
#   results/7_1.log                     — sink log (Caveat: when sourced via
#                                         `source(..., echo = TRUE)` the
#                                         on.exit + sink interaction captures
#                                         only the initial echo line; analytic
#                                         outputs above are produced correctly
#                                         regardless. Same caveat as §5.1
#                                         empirical-median + §6.2 logistic. To
#                                         capture the full transcript, run via
#                                         `Rscript` or capture R console output
#                                         separately.)
#   figures/7_1_figure_3a_sbp.png       — 5-panel trajectory plot, SBP-responder
#   figures/7_1_figure_3b_dbp.png       — 5-panel trajectory plot, DBP-responder

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(ggh4x)         # facetted_pos_scales(): per-panel y-axis limits + decimals
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")  # htn_baseline_on(), tidy_contrast(), emit_sizes()

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "7_1", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(out_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "7_1.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)
cat("§7.1 WHOOP-metric trajectories — dual-primary (SBP + DBP responder) —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read frame, derive DBP-responder indicator inline (parallel to SBP) ---
user_frame <- read_rds(input_path) %>%
  mutate(DBP_RESPONDER_6MO = as.integer((MO6_DBP - BASELINE_DBP) <= -5),
         HTN_BASELINE_ON   = htn_baseline_on(HTN_MED_BASELINE))

# --- Five WHOOP metrics and their column patterns ---
metrics <- c(
  RHR              = "RHR",
  SleepDuration    = "SLEEP_DURATION",
  SleepConsistency = "SLEEP_CONSISTENCY",
  Strain           = "STRAIN",
  HRV              = "HRV"
)
metric_units <- c(
  RHR              = "bpm",
  SleepDuration    = "hours",
  SleepConsistency = "% (0-100)",
  Strain           = "0-21 scale",
  HRV              = "ms"
)

# Favorable direction for responder group at MO6 — invariant to responder_def.
# Cardiovascular literature: lower RHR favorable; higher sleep duration,
# sleep consistency, strain (activity), HRV favorable.
metric_favorable_direction <- c(
  RHR              = "lower",
  SleepDuration    = "higher",
  SleepConsistency = "higher",
  Strain           = "higher",
  HRV              = "higher"
)

# --- Cohort builder: parameterized by responder_def ∈ {"SBP", "DBP"} ---
build_cohort <- function(d, responder_def) {
  resp_col <- paste0(responder_def, "_RESPONDER_6MO")
  d %>%
    filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
           SEX %in% c("male", "female"),
           !is.na(.data[[resp_col]])) %>%
    mutate(responder_group = factor(
      ifelse(.data[[resp_col]] == 1L, "responder", "non_responder"),
      levels = c("non_responder", "responder")
    ))
}

# --- Helper: build long-format dataframe for one metric ---
build_long <- function(metric_prefix, d) {
  cols <- paste0(metric_prefix, c("_BASELINE_MEAN", "_MO3_MEAN", "_MO6_MEAN"))
  d %>%
    select(STUDY_ID, responder_group, AGE_AT_START, SEX, BASELINE_BMI, HTN_BASELINE_ON,
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
  cat(sprintf("\n=== %s (%s) ===\n", metric_name, metric_prefix))
  cat("Long-format rows:", nrow(long_d), "; unique users:", n_distinct(long_d$STUDY_ID), "\n")
  cat("Rows per (Time × group):\n")
  print(long_d %>% count(Time, responder_group))
  f <- as.formula(
    "value ~ Time * responder_group + AGE_AT_START + SEX + BASELINE_BMI + HTN_BASELINE_ON + (1 | STUDY_ID)"
  )
  fit <- lmer(f, data = long_d, REML = TRUE)
  list(fit = fit, data = long_d, metric = metric_name)
}

# --- Extraction helpers (per-metric → tibble rows) ---
extract_emm <- function(fr) {
  emm    <- suppressMessages(emmeans(fr$fit, ~ Time | responder_group))
  delta  <- suppressMessages(contrast(emm, method = "trt.vs.ctrl", adjust = "none"))
  consec <- suppressMessages(contrast(emm, method = "consec",      adjust = "none"))
  bg_emm  <- suppressMessages(emmeans(fr$fit, ~ responder_group | Time))
  between <- suppressMessages(contrast(bg_emm, method = "revpairwise", adjust = "none"))
  # tidy_contrast() (from cohort_defs.R) attaches 95% CIs + within-family Holm p
  # (p.holm = p.adjust(p.value, method = "holm") within each contrast object).
  # Families: delta = 4 contrasts (2 timepoints × 2 groups);
  #           consec = 4; between = 3 (one per timepoint). Per SAP v0.84.
  list(
    emmeans = as_tibble(as.data.frame(emm)) %>% mutate(metric = fr$metric),
    delta   = tidy_contrast(delta)  %>% mutate(metric = fr$metric),
    consec  = tidy_contrast(consec) %>% mutate(metric = fr$metric),
    between = tidy_contrast(between) %>% mutate(metric = fr$metric)
  )
}

# =============================================================================
# Outer loop: SBP-responder then DBP-responder (both primary per SAP v0.64)
# =============================================================================
responder_defs <- c("SBP", "DBP")

fits_by_def       <- list()
emm_results_by_def <- list()

anova_all_rows     <- list()
emmeans_all_rows   <- list()
delta_all_rows     <- list()
consec_all_rows    <- list()
between_all_rows   <- list()
r2_all_rows        <- list()
direction_all_rows <- list()

cohort_sizes <- tibble(responder_def = character(),
                       n_users = integer(),
                       n_responder = integer(),
                       n_non_responder = integer())

for (rd in responder_defs) {
  cat("\n\n############################################################\n")
  cat(sprintf("##  PRIMARY FRAME: %s-responder grouping\n", rd))
  cat("############################################################\n")

  cohort <- build_cohort(user_frame, rd)
  cat(sprintf("§2.2 cohort + %s-responder defined: %d users\n", rd, nrow(cohort)))
  cat("Responder distribution:\n")
  resp_counts <- cohort %>% count(responder_group)
  print(resp_counts)

  cohort_sizes <- bind_rows(cohort_sizes, tibble(
    responder_def    = rd,
    n_users          = nrow(cohort),
    n_responder      = resp_counts$n[resp_counts$responder_group == "responder"],
    n_non_responder  = resp_counts$n[resp_counts$responder_group == "non_responder"]
  ))

  fits <- list()
  for (m in names(metrics)) {
    fits[[m]] <- fit_metric_lmm(m, metrics[[m]], cohort)
  }
  fits_by_def[[rd]] <- fits

  # ANOVA per metric
  anova_rd <- map_dfr(fits, function(fr) {
    a <- suppressMessages(anova(fr$fit, type = 3))
    as_tibble(a, rownames = "term") %>% mutate(metric = fr$metric)
  }) %>%
    mutate(responder_def = rd) %>%
    select(responder_def, metric, term, everything())
  anova_all_rows[[rd]] <- anova_rd

  cat("\n--- Type III ANOVA: Time × responder_group p-values per metric ---\n")
  print(anova_rd %>%
          filter(term == "Time:responder_group") %>%
          select(responder_def, metric, `F value`, `Pr(>F)`))

  # R²
  r2_rd <- map_dfr(fits, function(fr) {
    r2 <- suppressMessages(r2_nakagawa(fr$fit))
    tibble(metric         = fr$metric,
           marginal_r2    = as.numeric(r2$R2_marginal),
           conditional_r2 = as.numeric(r2$R2_conditional))
  }) %>% mutate(responder_def = rd) %>%
    select(responder_def, everything())
  r2_all_rows[[rd]] <- r2_rd

  cat("\n--- Nakagawa R² ---\n")
  print(r2_rd)

  # emmeans + delta
  emm_results <- map(fits, extract_emm)
  emm_results_by_def[[rd]] <- emm_results

  emmeans_rd <- map_dfr(emm_results, ~ .x$emmeans) %>%
    mutate(responder_def = rd) %>% select(responder_def, everything())
  delta_rd <- map_dfr(emm_results, ~ .x$delta) %>%
    mutate(responder_def = rd) %>% select(responder_def, everything())
  consec_rd <- map_dfr(emm_results, ~ .x$consec) %>%
    mutate(responder_def = rd) %>% select(responder_def, everything())
  between_rd <- map_dfr(emm_results, ~ .x$between) %>%
    mutate(responder_def = rd) %>% select(responder_def, everything())
  emmeans_all_rows[[rd]] <- emmeans_rd
  delta_all_rows[[rd]]   <- delta_rd
  consec_all_rows[[rd]]  <- consec_rd
  between_all_rows[[rd]] <- between_rd

  cat("\n--- Per-group Δ from baseline (MO3, MO6) per metric ---\n")
  print(delta_rd %>%
          select(responder_def, metric, responder_group, contrast, estimate, SE, p.value),
        n = Inf, width = Inf)

  # Directionality at MO6
  direction_rd <- delta_rd %>%
    filter(contrast == "6mo - baseline") %>%
    select(responder_def, metric, responder_group, estimate) %>%
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
    select(responder_def, metric, non_responder, responder, resp_minus_nonresp,
           favorable_direction, favorable_concordant)
  direction_all_rows[[rd]] <- direction_rd

  cat("\n--- Directionality check at MO6 (responder − non_responder, vs favorable direction) ---\n")
  print(direction_rd, n = Inf, width = Inf)
  cat(sprintf("\n  → %d/%d metrics direction-concordant for %s-responder.\n",
              sum(direction_rd$favorable_concordant), nrow(direction_rd), rd))
}

# --- Combine + write CSVs (with responder_def column) ---
anova_rows     <- bind_rows(anova_all_rows)
emmeans_tbl    <- bind_rows(emmeans_all_rows)
delta_tbl      <- bind_rows(delta_all_rows)
consec_tbl     <- bind_rows(consec_all_rows)
between_tbl    <- bind_rows(between_all_rows)
r2_rows        <- bind_rows(r2_all_rows)
direction_check <- bind_rows(direction_all_rows)

write_csv(anova_rows,      file.path(results_dir, "7_1_anova.csv"))
write_csv(emmeans_tbl,     file.path(results_dir, "7_1_emmeans.csv"))
write_csv(delta_tbl,       file.path(results_dir, "7_1_delta_contrasts.csv"))
write_csv(consec_tbl,      file.path(results_dir, "7_1_consec_contrasts.csv"))        # adjacent-timepoint (line solid/dashed)
write_csv(between_tbl,     file.path(results_dir, "7_1_between_group_contrasts.csv")) # responder vs non at each timepoint (asterisk)
write_csv(r2_rows,         file.path(results_dir, "7_1_r2.csv"))
write_csv(direction_check, file.path(results_dir, "7_1_direction_check.csv"))
write_csv(cohort_sizes,    file.path(results_dir, "7_1_cohort_sizes.csv"))

# --- Save fits + emmeans (lists keyed by responder_def) ---
write_rds(fits_by_def,        file.path(out_dir, "fitted_models.rds"))
write_rds(emm_results_by_def, file.path(out_dir, "emmeans_results.rds"))

# --- Emit standardized sizes for Table A (per SAP review v0.84) ---
sizes_71 <- bind_rows(lapply(responder_defs, function(rd) {
  bind_rows(lapply(names(metrics), function(m) {
    fr  <- fits_by_def[[rd]][[m]]
    cs  <- cohort_sizes %>% filter(responder_def == rd)
    tibble(
      section       = "7.1",
      role          = "primary",
      analysis      = "WHOOP trajectory LMM",
      cohort        = "hypertensive",
      outcome       = m,
      predictor_set = rd,
      n_users       = as.integer(n_distinct(fr$data$STUDY_ID)),
      n_obs         = as.integer(nrow(fr$data)),
      n_events      = NA_integer_,
      ceiling_n     = as.integer(cs$n_users)
    )
  }))
}))
emit_sizes(sizes_71, "7_1", results_dir)

# =============================================================================
# Figure 3 — three variants
# =============================================================================

# Build the long-format figure dataframe with CIs and shared cosmetics.
fig_df_all <- emmeans_tbl %>%
  mutate(
    lcl = if ("asymp.LCL" %in% names(.)) asymp.LCL else lower.CL,
    ucl = if ("asymp.UCL" %in% names(.)) asymp.UCL else upper.CL
  )
if (!"lcl" %in% names(fig_df_all)) {
  fig_df_all <- emmeans_tbl %>%
    mutate(
      crit = if ("df" %in% names(.)) ifelse(is.finite(df), qt(0.975, df), 1.959964) else 1.959964,
      lcl  = emmean - crit * SE,
      ucl  = emmean + crit * SE
    )
}

fig_df_all <- fig_df_all %>%
  mutate(
    Time = factor(Time, levels = c("baseline", "3mo", "6mo")),
    time_num = case_when(
      Time == "baseline" ~ 0,
      Time == "3mo"      ~ 3,
      Time == "6mo"      ~ 6
    ),
    metric_label = factor(metric,
                          levels = c("RHR", "SleepConsistency", "SleepDuration", "Strain", "HRV"),
                          labels = c("RHR (bpm)", "Sleep Consistency (%)",
                                     "Sleep Duration (hours)", "Physical Activity (a.u.)", "HRV (ms)"))
  )

# Per-panel y scales — matched by metric_label factor-level order
# (RHR, Sleep Consistency, Sleep Duration, Physical Activity, HRV).
# 3a (SBP) uses the tight requested ranges; the SBP CI fits all but Sleep
# Consistency (dips to ~69.46), so that one limit is nudged out minimally to
# avoid clipping. 3b (DBP) + combined use slightly wider limits (union of both
# responder cohorts) so no CI clips there either.
y_scales_sbp <- list(
  scale_y_continuous(limits = c(62.75, 64.25), breaks = c(63, 64)),                 # RHR
  scale_y_continuous(limits = c(69.45, 70.55), breaks = c(69.5, 70, 70.5)),         # Sleep Consistency (nudged from 69.5–70.5)
  scale_y_continuous(limits = c(7.0, 7.15),    breaks = c(7.00, 7.05, 7.10, 7.15)), # Sleep Duration
  scale_y_continuous(limits = c(8, 9.5),       breaks = c(8, 8.5, 9, 9.5)),         # Physical Activity
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14)))                      # HRV: auto + asterisk headroom
)
y_scales_wide <- list(
  scale_y_continuous(limits = c(62.5, 64.4), breaks = c(63, 64)),                   # RHR
  scale_y_continuous(limits = c(69.2, 70.6), breaks = c(69.5, 70, 70.5)),           # Sleep Consistency
  scale_y_continuous(limits = c(7.0, 7.15),  breaks = c(7.00, 7.05, 7.10, 7.15)),   # Sleep Duration
  scale_y_continuous(limits = c(8, 9.6),     breaks = c(8, 8.5, 9, 9.5)),           # Physical Activity
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14)))                      # HRV: auto + asterisk headroom
)

# --- Figure 3a: SBP-responder (parallel to legacy figure) ---
make_parallel_figure <- function(df_def, responder_def, out_png, y_scales) {
  rdef <- responder_def
  resp_label <- sprintf("%s-Responder", responder_def)
  resp_color <- if (responder_def == "DBP") "#33a02c" else "#1f78b4"  # DBP green matches §6.2
  alpha_sig  <- 0.005
  time_of    <- c("baseline" = 0, "3mo" = 3, "6mo" = 6)

  # Point fill: significant vs BASELINE (delta / trt.vs.ctrl); baseline = empty (reference).
  # Uses Holm-adjusted p (p.holm) per SAP v0.84.
  pt_sig <- delta_tbl %>% filter(responder_def == rdef) %>%
    transmute(metric, responder_group,
              time_num = unname(time_of[sub(" .*$", "", contrast)]),
              sig_base = p.holm < alpha_sig)
  df_plot <- df_def %>%
    mutate(responder_label = factor(responder_group,
                                    levels = c("non_responder", "responder"),
                                    labels = c("Non-Responder", resp_label))) %>%
    left_join(pt_sig, by = c("metric", "responder_group", "time_num")) %>%
    mutate(pt_filled = dplyr::coalesce(sig_base, FALSE))   # baseline (no row) -> empty
  group_cols <- setNames(c("#999999", resp_color),
                         c("Non-Responder", resp_label))

  # Line segments: solid if significant change vs previous timepoint (consec), else dashed.
  em_lk <- df_plot %>% distinct(metric, metric_label, responder_group, responder_label, time_num, emmean)
  seg <- consec_tbl %>% filter(responder_def == rdef) %>%
    tidyr::separate(contrast, into = c("to", "from"), sep = " - ", remove = TRUE) %>%
    transmute(metric, responder_group,
              x = unname(time_of[from]), xend = unname(time_of[to]),
              seg_sig = p.holm < alpha_sig) %>%
    left_join(em_lk, by = c("metric", "responder_group", "x" = "time_num")) %>% rename(y = emmean) %>%
    left_join(em_lk %>% select(metric, responder_group, time_num, emmean),
              by = c("metric", "responder_group", "xend" = "time_num")) %>% rename(yend = emmean)

  # Asterisk: significant Responder vs Non-Responder difference at a timepoint,
  # placed just above the higher of the two CIs in that (metric × timepoint).
  ast <- between_tbl %>% filter(responder_def == rdef, p.holm < alpha_sig) %>%
    transmute(metric, time_num = unname(time_of[Time])) %>%
    left_join(df_def %>% group_by(metric, metric_label, time_num) %>%
                summarise(y = max(ucl), .groups = "drop"),
              by = c("metric", "time_num"))

  p <- ggplot(df_plot, aes(x = time_num, y = emmean,
                           color = responder_label, fill = responder_label,
                           group = responder_label)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.18, color = NA) +
    geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend,
                                 color = responder_label, linetype = seg_sig),
                 linewidth = 0.9, inherit.aes = FALSE) +
    geom_point(aes(shape = pt_filled), size = 2.5) +
    geom_text(data = ast, aes(x = time_num, y = y, label = "*"),
              vjust = 0, size = 6, fontface = "bold", color = "black", inherit.aes = FALSE) +
    facet_manual(vars(metric_label),
                 # 4-row grid so every panel spans 2 rows (equal size); Physical
                 # Activity (right) occupies the middle two rows → vertically centred.
                 design = "AE#\nAED\nCBD\nCB#",
                 scales = "free_y", strip.position = "left", axes = "x") +
    facetted_pos_scales(y = y_scales) +
    scale_x_continuous(breaks = c(0, 3, 6),
                       labels = c("Baseline", "3 Months", "6 Months")) +
    scale_color_manual(values = group_cols, name = NULL,
                       breaks = c(resp_label, "Non-Responder")) +
    scale_fill_manual(values  = group_cols, name = NULL,
                      breaks = c(resp_label, "Non-Responder")) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed"), guide = "none") +
    labs(
      x = NULL, y = NULL,
      title    = "Six-Month Wearable Derived Biometric Trajectories by Responder Status",
      subtitle = "Adjusted LMM estimated marginal means with 95% CI"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid        = element_blank(),
      panel.spacing     = grid::unit(1.5, "lines"),
      plot.margin       = margin(t = 6, r = 16, b = 6, l = 6),   # room for the last "6 Months" label
      axis.line         = element_line(color = "black", linewidth = 0.4),
      axis.ticks        = element_line(color = "black", linewidth = 0.4),
      axis.text.x       = element_text(face = "bold"),
      strip.placement   = "outside",
      strip.background  = element_blank(),
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 11),
      legend.position   = "bottom",
      plot.title        = element_text(face = "bold")
    )
  ggsave(out_png, p, width = 11, height = 7.5, dpi = 300)
  cat("Wrote figure:\n  ", out_png, "\n")
}

make_parallel_figure(
  fig_df_all %>% filter(responder_def == "SBP"),
  "SBP",
  file.path(figures_dir, "7_1_figure_3a_sbp.png"),
  y_scales = y_scales_sbp
)
make_parallel_figure(
  fig_df_all %>% filter(responder_def == "DBP"),
  "DBP",
  file.path(figures_dir, "7_1_figure_3b_dbp.png"),
  y_scales = y_scales_wide
)

# (Combined 4-line figure removed 2026-06-01 — not used; 3a/3b are the keepers.)

# =============================================================================
# Paper-ready Table 4 — parallel SBP / DBP sections
# =============================================================================
fmt_signed <- function(x, d = 3) sprintf("%+.*f", d, x)
fmt_p <- function(p) {
  ifelse(p < 1e-10, formatC(p, format = "e", digits = 1),
  ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                    sprintf("%.3f", p)))
}

md_lines <- c(
  "# §7.1 Table 4 — WHOOP-metric trajectories by responder status",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "Dual-primary analysis per SAP §7.1 v0.64: SBP-responder and DBP-responder are both primary groupings, each with its own §2.2 cohort filtered to non-missing responder indicator."
)

for (rd in responder_defs) {
  cs <- cohort_sizes %>% filter(responder_def == rd)
  md_lines <- c(md_lines, "",
                paste0("## ", rd, "-responder primary frame"),
                "",
                paste0("§2.2 cohort, n = ", format(cs$n_users, big.mark = ","),
                       " users with ", rd, "-responder defined ",
                       "(responders: ", format(cs$n_responder, big.mark = ","),
                       "; non-responders: ", format(cs$n_non_responder, big.mark = ","), ")."),
                "",
                "### Time × responder_group interaction (formal test of trajectory divergence)",
                "",
                "Bold = significant at α = 0.005 per SAP §9."
  )
  interaction_rows <- anova_rows %>%
    filter(responder_def == rd, term == "Time:responder_group") %>%
    arrange(metric)
  for (i in seq_len(nrow(interaction_rows))) {
    r <- interaction_rows[i, ]
    bold <- ifelse(r$`Pr(>F)` < 0.005, "**", "")
    md_lines <- c(md_lines, sprintf(
      "- **%s**: F = %.2f, %sp = %s%s",
      r$metric, r$`F value`, bold, fmt_p(r$`Pr(>F)`), bold
    ))
  }

  md_lines <- c(md_lines, "",
                "### Per-group Δ from baseline (MO3, MO6) per metric")
  for (m in names(metrics)) {
    md_lines <- c(md_lines, "", paste0("#### ", m, " (", metric_units[m], ")"))
    for (g in c("non_responder", "responder")) {
      rows <- delta_tbl %>%
        filter(responder_def == rd, metric == m, responder_group == g)
      md_lines <- c(md_lines, paste0("**", g, "**"))
      for (i in seq_len(nrow(rows))) {
        r <- rows[i, ]
        bold <- ifelse(r$p.holm < 0.005, "**", "")
        md_lines <- c(md_lines, sprintf(
          "- %s: %sΔ = %s (SE %.4f, p = %s, p.holm = %s)%s",
          r$contrast, bold, fmt_signed(r$estimate),
          r$SE, fmt_p(r$p.value), fmt_p(r$p.holm), bold
        ))
      }
    }
  }

  md_lines <- c(md_lines, "",
                "### Directionality at MO6 (responder − non_responder vs favorable direction)")
  dir_rows <- direction_check %>% filter(responder_def == rd) %>% arrange(metric)
  for (i in seq_len(nrow(dir_rows))) {
    r <- dir_rows[i, ]
    mark <- ifelse(isTRUE(r$favorable_concordant), "✓", "✗")
    md_lines <- c(md_lines, sprintf(
      "- %s: responder − non_responder = %s, favorable = %s %s",
      r$metric, fmt_signed(r$resp_minus_nonresp), r$favorable_direction, mark
    ))
  }
  md_lines <- c(md_lines, "",
                sprintf("  → %d/%d metrics direction-concordant.",
                        sum(dir_rows$favorable_concordant), nrow(dir_rows)))
}

md_lines <- c(md_lines, "",
  "## Notes",
  "",
  "- **Formal test**: `Time × responder_group` interaction p-value (Type III ANOVA, Satterthwaite df). Significant interaction = trajectory shape differs between responder groups.",
  "- **Cohort**: §2.2 hypertensive + sex restriction + responder defined (ΔBP ≤ −5 mmHg at MO6, SBP or DBP per frame). Users contribute up to 3 rows per metric (baseline / 3mo / 6mo); unbalanced data handled by REML.",
  "- **Time windows**: baseline = 30 days prior to STUDY_START_DATE; 3mo = 30-day window centered on day 90 (days 75–105); 6mo = 30-day window centered on day 180 (days 165–195). Window means already constructed in user frame.",
  "- **Significance**: α = 0.005 per SAP §9. Post-gate contrasts (Δ from baseline, consecutive Δ, between-group) use **within-family Holm-adjusted p** (`p.holm`) for significance flags and encoding on Figure 3 — consistent with §5.1/§5.5 universal Holm convention on contrast families (per SAP v0.84). Raw p retained alongside `p.holm` in CSVs and this table.",
  "- **Dual-primary framing**: SBP- and DBP-responder are both primary per SAP v0.64 (post-audit promotion; see SAP_CHANGELOG entry v0.64 and DEVIATIONS.md for framing context). HTN-medication-naive sensitivity is run against both primary frames in `7_1_sensitivities.R`.",
  "- Pairs with: §6.1 (LE8 predictors of BP change); §6.2 (WHOOP-signal predictors of SBP-responder); §7.2 (W/B-decomposed daily LMM)."
)

writeLines(md_lines, file.path(results_dir, "7_1_paper_table4.md"))
cat("\nWrote markdown:\n  ", file.path(results_dir, "7_1_paper_table4.md"), "\n")

cat("\nNext: 7_1_diagnostics.R (Q-Q residuals, fitted/residual + LOESS, BLUP Q-Q per fit — now 10 fits).\n")
cat("\n=== §7.1 WHOOP trajectories done (dual-primary) ===\n")
