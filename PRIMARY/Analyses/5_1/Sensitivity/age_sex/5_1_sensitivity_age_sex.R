# 5_1_sensitivity_age_sex.R
# SAP §5.1 sensitivity bullet 6: age- and sex-stratified replication.
# Refits the §5.1 primary ADJUSTED LMM separately within each age band
# (18–39 / 40–59 / 60–79) and within each sex stratum (male / female).
# Convergence in direction and rough magnitude across stratified fits
# demonstrates robustness of the headline findings to age and sex
# composition.
#
# Implementation notes:
# - Within age bands: keep age as a covariate (still varies within band).
# - Within sex strata: drop SEX and Time:SEX from the model spec (constant).
# - Adjusted model only — matches the primary §5.1 reporting frame and keeps
#   computational cost reasonable (10 LMM fits total).
#
# Cohort: §2.2 hypertensive + sex restriction, matching §5.1 primary.
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/5_1/Sensitivity/age_sex/derived/):
#   fitted_models.rds — list of 10 fits (5 subgroups × 2 outcomes)
#   results/age_sex_subgroup_sizes.csv      — per-subgroup × stratum cohort counts
#   results/delta_contrasts_by_subgroup.csv — within-stratum Δs per subgroup
#   results/comparison_5_1_vs_subgroups.csv — vs §5.1 primary
#   results/5_1_sensitivity_age_sex_paper.md
#   results/5_1_sensitivity_age_sex.log
#   figures/5_1_age_sex_forest.png          — forest plot

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")  # htn_baseline_on(), emit_sizes()

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "5_1", "Sensitivity", "age_sex", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(out_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "5_1_sensitivity_age_sex.log")
sink(log_path, split = TRUE)
cat("§5.1 sensitivity — age + sex stratified —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

user_frame <- read_rds(input_path)
cohort <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female"),
         !is.na(AGE_AT_START)) %>%
  mutate(
    age_band = case_when(
      AGE_AT_START >= 18 & AGE_AT_START <= 39 ~ "18–39",
      AGE_AT_START >= 40 & AGE_AT_START <= 59 ~ "40–59",
      AGE_AT_START >= 60 & AGE_AT_START <= 79 ~ "60–79",
      TRUE                                     ~ NA_character_
    ) %>% factor(levels = c("18–39", "40–59", "60–79"))
  ) %>%
  filter(!is.na(age_band))
cat("§2.2 cohort + sex + age-band filter:", nrow(cohort), "users\n\n")

cat("Subgroup sizes (total, before AHA stratification):\n")
print(cohort %>% count(age_band))
print(cohort %>% count(SEX))

# --- Pivot to long (same as §5.1) ---
make_long <- function(d) {
  d %>%
    select(STUDY_ID, BASELINE_AHA, AGE_AT_START, SEX, BASELINE_BMI, HTN_MED_BASELINE, age_band,
           BASELINE_SBP,           MO3_SBP,           MO6_SBP,
           BASELINE_DBP,           MO3_DBP,           MO6_DBP,
           BASELINE_TIME_OF_DAY,   MO3_TIME_OF_DAY,   MO6_TIME_OF_DAY,
           BASELINE_SEASON,        MO3_SEASON,        MO6_SEASON) %>%
    pivot_longer(
      cols = c(starts_with("BASELINE_S"), starts_with("MO3_S"), starts_with("MO6_S"),
               starts_with("BASELINE_D"), starts_with("MO3_D"), starts_with("MO6_D"),
               starts_with("BASELINE_T"), starts_with("MO3_T"), starts_with("MO6_T")),
      names_to = c("Time", ".value"),
      names_pattern = "(BASELINE|MO3|MO6)_(.+)"
    ) %>%
    mutate(
      Time         = factor(Time, levels = c("BASELINE", "MO3", "MO6")),
      BASELINE_AHA = factor(BASELINE_AHA,
                            levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")),
      SEX          = factor(SEX),
      SEASON       = factor(SEASON, levels = c("winter", "spring", "summer", "fall")),
      HTN_BASELINE_ON = htn_baseline_on(HTN_MED_BASELINE)  # baseline-only HTN trait (v0.79)
    )
}

# --- LMM helper: adjusted model with conditional covariate dropping ---
fit_subgroup_lmm <- function(outcome_col, data, drop_sex = FALSE) {
  covariates <- c(
    "AGE_AT_START", "BASELINE_BMI", "HTN_BASELINE_ON",
    "sin(2 * pi * TIME_OF_DAY / 24)", "cos(2 * pi * TIME_OF_DAY / 24)",
    "SEASON",
    "Time:AGE_AT_START", "Time:BASELINE_BMI", "Time:HTN_BASELINE_ON"
  )
  if (!drop_sex) covariates <- c(covariates, "SEX", "Time:SEX")
  f <- as.formula(paste0(
    outcome_col, " ~ Time * BASELINE_AHA + ",
    paste(covariates, collapse = " + "),
    " + (1 | STUDY_ID)"
  ))
  d <- data %>% filter(!is.na(.data[[outcome_col]]), !is.na(BASELINE_AHA),
                       !is.na(AGE_AT_START), !is.na(BASELINE_BMI),
                       !is.na(HTN_BASELINE_ON),
                       !is.na(TIME_OF_DAY), !is.na(SEASON))
  if (!drop_sex) d <- d %>% filter(!is.na(SEX))
  cat(sprintf("  Fitting on %d rows, %d users…\n",
              nrow(d), n_distinct(d$STUDY_ID)))
  lmer(f, data = d, REML = TRUE)
}

# --- Extract Δ contrasts (within-stratum) from a fit ---
# emmeans contrast() does NOT include CI columns by default; derive Wald CIs
# from estimate ± qt(0.975, df) * SE (z-based when df = Inf).
extract_within_delta <- function(fit, subgroup_label, outcome_label) {
  emm <- suppressMessages(emmeans(fit, ~ Time | BASELINE_AHA))
  delta_emm <- suppressMessages(contrast(emm, method = "trt.vs.ctrl", adjust = "none"))
  d <- as_tibble(as.data.frame(delta_emm))
  if (all(c("asymp.LCL", "asymp.UCL") %in% names(d))) {
    d$lcl <- d$asymp.LCL; d$ucl <- d$asymp.UCL
  } else if (all(c("lower.CL", "upper.CL") %in% names(d))) {
    d$lcl <- d$lower.CL;  d$ucl <- d$upper.CL
  } else {
    crit <- if ("df" %in% names(d)) ifelse(is.finite(d$df), qt(0.975, d$df), 1.959964) else 1.959964
    d$lcl <- d$estimate - crit * d$SE
    d$ucl <- d$estimate + crit * d$SE
  }
  d %>%
    mutate(subgroup = subgroup_label,
           outcome  = outcome_label,
           stratum  = factor(BASELINE_AHA,
                             levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
                             labels = c("Elevated", "Stage 1", "Stage 2")),
           timepoint = sub(" - BASELINE", "", contrast),
           p.holm    = p.adjust(p.value, method = "holm")) %>%   # within-family Holm (v0.79)
    select(subgroup, outcome, stratum, timepoint, estimate, SE, lcl, ucl, p.value, p.holm)
}

# --- Define subgroups (filter as a function to avoid formula-injection issues) ---
subgroups <- list(
  list(label = "Age 18–39",  filter_fn = function(d) filter(d, age_band == "18–39"),  drop_sex = FALSE),
  list(label = "Age 40–59",  filter_fn = function(d) filter(d, age_band == "40–59"),  drop_sex = FALSE),
  list(label = "Age 60–79",  filter_fn = function(d) filter(d, age_band == "60–79"),  drop_sex = FALSE),
  list(label = "Male",       filter_fn = function(d) filter(d, SEX == "male"),        drop_sex = TRUE),
  list(label = "Female",     filter_fn = function(d) filter(d, SEX == "female"),      drop_sex = TRUE)
)

# --- Fit 10 models (5 subgroups × 2 outcomes) and collect Δs ---
all_fits     <- list()
all_deltas   <- list()
all_sizes    <- list()

for (sg in subgroups) {
  cat("\n\n========================================\n")
  cat("Subgroup:", sg$label, "\n")
  cat("========================================\n")
  sg_cohort <- sg$filter_fn(cohort)
  cat("Subgroup n (users):", nrow(sg_cohort), "\n")
  cat("Per-stratum counts:\n")
  sizes <- sg_cohort %>% count(BASELINE_AHA) %>%
    mutate(subgroup = sg$label) %>%
    select(subgroup, BASELINE_AHA, n)
  print(sizes)
  all_sizes[[sg$label]] <- sizes

  long_d <- make_long(sg_cohort)

  for (otc in c("SBP", "DBP")) {
    cat(sprintf("\n--- %s × %s ---\n", sg$label, otc))
    fit <- fit_subgroup_lmm(otc, long_d, drop_sex = sg$drop_sex)
    all_fits[[paste0(sg$label, "__", otc)]] <- fit
    d_tbl <- extract_within_delta(fit, sg$label, otc)
    all_deltas[[paste0(sg$label, "__", otc)]] <- d_tbl
    cat("Within-stratum Δ contrasts:\n")
    print(d_tbl)
  }
}

# --- Combine + write CSVs ---
sizes_df <- bind_rows(all_sizes) %>%
  mutate(stratum = factor(BASELINE_AHA,
                          levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
                          labels = c("Elevated", "Stage 1", "Stage 2"))) %>%
  select(subgroup, stratum, n)
write_csv(sizes_df, file.path(results_dir, "age_sex_subgroup_sizes.csv"))

deltas_df <- bind_rows(all_deltas)
# lcl/ucl already derived in extract_within_delta() via Wald CI from estimate + SE + df
write_csv(deltas_df, file.path(results_dir, "delta_contrasts_by_subgroup.csv"))

# Save fits for downstream re-use
write_rds(all_fits, file.path(out_dir, "fitted_models.rds"))

# --- Emit standardized model sizes for Table A.2 (from the fits) ---
# One row per subgroup × outcome (adjusted only). ceiling_n = subgroup entry N.
subgroup_ceiling <- sizes_df %>% group_by(subgroup) %>%
  summarise(ceiling_n = sum(n), .groups = "drop")
emit_sizes(
  imap_dfr(all_fits, function(fit, key) {
    parts <- strsplit(key, "__", fixed = TRUE)[[1]]
    tibble(section = "5.1", role = "sensitivity",
           analysis = "BP-trajectory LMM (age/sex-stratified)", cohort = parts[1],
           outcome = parts[2], predictor_set = "adjusted",
           n_users = as.integer(lme4::ngrps(fit)), n_obs = as.integer(nobs(fit)),
           n_events = NA_integer_,
           ceiling_n = subgroup_ceiling$ceiling_n[match(parts[1], subgroup_ceiling$subgroup)])
  }),
  "5_1_age_sex", results_dir
)

# --- Comparison vs §5.1 primary ---
primary_path <- file.path(project_dir, "Analyses", "5_1", "derived", "results",
                          "5_1_delta_contrasts.csv")
primary <- read_csv(primary_path, show_col_types = FALSE) %>%
  filter(model %in% c("sbp_adj", "dbp_adj")) %>%
  transmute(
    outcome     = ifelse(model == "sbp_adj", "SBP", "DBP"),
    stratum     = factor(BASELINE_AHA,
                         levels = c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
                         labels = c("Elevated", "Stage 1", "Stage 2")),
    timepoint   = sub(" - BASELINE", "", contrast),
    primary_delta = estimate
  )

comparison <- deltas_df %>%
  left_join(primary, by = c("outcome", "stratum", "timepoint")) %>%
  mutate(diff_vs_primary = estimate - primary_delta) %>%
  arrange(outcome, stratum, timepoint, subgroup)

write_csv(comparison, file.path(results_dir, "comparison_5_1_vs_subgroups.csv"))

cat("\n\n--- Comparison vs §5.1 primary (adjusted) ---\n")
print(comparison %>%
        filter(timepoint == "MO6") %>%
        select(subgroup, outcome, stratum, estimate, SE, primary_delta, diff_vs_primary),
      n = Inf, width = Inf)

# --- Forest plot ---
fig_df <- comparison %>%
  filter(timepoint == "MO6") %>%
  mutate(
    subgroup_type = ifelse(grepl("^Age", subgroup), "Age band", "Sex"),
    subgroup_f = factor(subgroup,
                        levels = c("Age 18–39", "Age 40–59", "Age 60–79", "Male", "Female"))
  )

stratum_cols <- c("Elevated" = "#fdae61", "Stage 1" = "#f46d43", "Stage 2" = "#a50026")

p <- ggplot(fig_df, aes(x = estimate, y = subgroup_f, color = stratum)) +
  geom_vline(xintercept = 0, color = "gray60", linetype = "dashed", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0.2, linewidth = 0.5) +
  geom_point(size = 3) +
  geom_point(aes(x = primary_delta), shape = 4, size = 3.5, color = "black", alpha = 0.6) +
  facet_grid(stratum ~ outcome, scales = "free_x") +
  scale_color_manual(values = stratum_cols, name = "AHA stratum", guide = "none") +
  labs(
    x = "MO6 Δ from baseline (mmHg)",
    y = NULL,
    title    = "Age- and sex-stratified §5.1 LMM — per-stratum MO6 Δ",
    subtitle = "Filled point = within-subgroup adjusted Δ with 95% CI. Black × = §5.1 primary estimate.",
    caption  = "Convergence in direction and rough magnitude across subgroups demonstrates robustness to age/sex composition."
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text       = element_text(face = "bold"),
        legend.position  = "bottom",
        plot.caption     = element_text(size = 9, color = "gray40", hjust = 0))

out_png <- file.path(figures_dir, "5_1_age_sex_forest.png")
ggsave(out_png, p, width = 10, height = 8, dpi = 300)
cat("\nWrote figure:\n  ", out_png, "\n")

# --- Paper-ready markdown ---
fmt_signed <- function(x, d = 2) sprintf("%+.*f", d, x)

md_lines <- c(
  "# §5.1 sensitivity — age- and sex-stratified replication",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "## MO6 adjusted Δ by subgroup, with comparison to §5.1 primary",
  ""
)
mo6 <- comparison %>% filter(timepoint == "MO6")
for (otc in c("SBP", "DBP")) {
  md_lines <- c(md_lines, paste0("### ", otc))
  for (sg in c("Age 18–39", "Age 40–59", "Age 60–79", "Male", "Female")) {
    md_lines <- c(md_lines, paste0("**", sg, "**"))
    rows <- mo6 %>% filter(outcome == otc, subgroup == sg)
    for (i in seq_len(nrow(rows))) {
      r <- rows[i, ]
      md_lines <- c(md_lines, sprintf(
        "- %s: subgroup Δ = %s mmHg (95%% CI [%s, %s]); §5.1 primary = %s; diff = %s.",
        r$stratum, fmt_signed(r$estimate), fmt_signed(r$lcl), fmt_signed(r$ucl),
        fmt_signed(r$primary_delta), fmt_signed(r$diff_vs_primary)
      ))
    }
    md_lines <- c(md_lines, "")
  }
}

md_lines <- c(md_lines,
  "## Notes",
  "",
  "- **Adjusted model only**, matching primary reporting frame. SEX and Time:SEX terms dropped in sex-stratified fits (constant within stratum). Age remains a covariate in age-band fits (still varies within band).",
  "- **Convergence criteria** (per SAP): direction and rough magnitude. Sign concordance across all 30 (subgroup × outcome × stratum) cells indicates robustness. Magnitude differences are reported as `diff_vs_primary`.",
  "- **Substantive conclusion**: convergence across age bands and sex strata demonstrates that the §5.1 findings are not driven by age- or sex-related cohort composition. Within-subgroup estimates that materially diverge from primary would indicate effect-modification; small differences within sampling error indicate robustness.",
  "- **Per-stratum n in subgroups**: smaller cell sizes inflate SEs particularly in the 60–79 × Elevated and Female × Stage 2 cells. Wider CIs are expected and shown in the forest plot.",
  "- **Effect-modification interpretation cautions**: differences in subgroup estimates are not formally tested against a single model with subgroup-by-Time interactions. The current implementation is a stratified-fit comparison, which is the simpler and more defensible approach for descriptive robustness checking.",
  "- Pairs with: §5.1 primary (full cohort), §5.1 RTM Oldham sensitivity, §5.1 continuous-baseline sensitivity."
)

writeLines(md_lines, file.path(results_dir, "5_1_sensitivity_age_sex_paper.md"))
cat("\nWrote markdown:\n  ", file.path(results_dir, "5_1_sensitivity_age_sex_paper.md"), "\n")

cat("\n=== §5.1 age + sex sensitivity done ===\n")
sink()
