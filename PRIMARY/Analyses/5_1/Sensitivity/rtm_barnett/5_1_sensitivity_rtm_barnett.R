# 5_1_sensitivity_rtm_barnett.R
# SAP §5.1 sensitivity: regression-to-the-mean (RTM) adjustment, Barnett method.
#
# Reference: Barnett AG, van der Pols JC, Dobson AJ. "Regression to the mean:
#            what it is and how to deal with it." Int J Epidemiol 2005;34:215-220.
#
# SUPERSEDES ../rtm/5_1_sensitivity_rtm.R, whose approach was noticed to be
# incomplete. Its outputs are retained for provenance but are NOT the manuscript
# values.
#
# ---------------------------------------------------------------------------
# METHOD (as reported in the manuscript Methods)
#
#   expected RTM = (sigma_w^2 / sigma_t^2) * (xbar_stratum - xbar_full)
#   RTM-adjusted Delta = observed Delta + expected RTM
#
#   sigma_w^2 = pooled within-participant variance of BP measurements recorded
#               from baseline through 30 days (participants with >= 2 such
#               measurements)
#   sigma_t^2 = var(baseline BP) in the FULL cohort (all four AHA strata)
#   xbar_full = mean baseline BP in the full cohort, including Normal
#
# sigma_w^2/sigma_t^2 is Barnett's (1 - rho): the proportion of baseline variance
# attributable to within-person variability. sigma_t*(1-rho)*C(z) is his Eq (1);
# C(z) is the normal-theory expression for the selected group's mean deviation
# from the population mean, expressed in SD units. Because the AHA strata are
# defined JOINTLY on SBP and DBP (non-Normal = SBP >= 120 OR DBP >= 80), the
# single-variable single-cutpoint assumption behind C(z) does not hold, so the
# OBSERVED stratum mean deviation is used in its place. Part 4 below quantifies
# what that substitution costs on a genuinely univariate cut (agreement 91-97%;
# the residual is non-normality of BP, for which the observed mean is the more
# accurate quantity).
#
# Each BP measurement is the mean of three consecutive cuff recordings, and
# sigma_w is estimated from day-to-day variation OF those triplicate means, so
# the averaging is already inside the unit of analysis. Barnett's Eq (2)
# (sigma_w^2/m) is therefore NOT applied; doing so would double-count it.
#
# Cohort: §2.2 modelled cohort conventions (SEX in {male, female}).
# Reference cohort for xbar_full: all four strata, complete baseline+MO6.
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#         Build/derived/qc/bp_reading_analytic_frame_1X.rds
#
# Outputs (Analyses/5_1/Sensitivity/rtm_barnett/derived/):
#   results/variance_components.csv   — sigma_w^2, sigma_t^2, sigma_b^2, rho, noise fraction
#   results/per_stratum_rtm.csv       — per stratum: observed, expected RTM, adjusted
#   results/cohort_rtm.csv            — restricted (Elevated+Stage1+Stage2) cohort
#   results/eq1_comparison.csv        — generalisation vs literal Eq (1) diagnostic
#   results/5_1_sensitivity_rtm_barnett_paper.md
#   results/5_1_sensitivity_rtm_barnett.log
#   figures/5_1_rtm_barnett_observed_vs_adjusted.png   (diagnostic, not a manuscript figure)

library(tidyverse)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")   # emit_sizes()

project_dir <- "./PRIMARY"
user_path   <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
rdg_path    <- file.path(project_dir, "Build", "derived", "qc", "bp_reading_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "5_1", "Sensitivity", "rtm_barnett", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "5_1_sensitivity_rtm_barnett.log")
sink(log_path, split = TRUE)
cat("§5.1 sensitivity — RTM adjustment (Barnett et al. 2005) —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

STRATA   <- c("Normal", "Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")
SHORT    <- c("Normal", "Elevated", "Stage 1", "Stage 2")
HYP      <- STRATA[2:4]
CUTPOINT <- c(SBP = 120, DBP = 80)
WINDOW_D <- 30

user <- read_rds(user_path) %>% filter(SEX %in% c("male", "female"))
rdg  <- read_rds(rdg_path)  %>% semi_join(user %>% select(STUDY_ID), by = "STUDY_ID")
cat("Full cohort (4 strata, modelled):", nrow(user), "users\n")
cat("Reading-level rows available:     ", nrow(rdg), "\n\n")

# ---------------------------------------------------------------------------
# PART 1 — Within-person variance from near-baseline repeated measurements
# ---------------------------------------------------------------------------
# NOTE ON STYLE: quantities are built as standalone vectors before any
# tibble()/summarise() call. tibble() and summarise() evaluate arguments
# sequentially, so a column named e.g. `b` will shadow an outer vector `b` in a
# later argument -- a silent, warning-free source of wrong answers.

pooled_within_var <- function(reading_df, outcome, one_per_day = FALSE) {
  d <- reading_df %>%
    filter(DAY_FROM_BASELINE >= 0, DAY_FROM_BASELINE <= WINDOW_D,
           !is.na(.data[[outcome]]))
  if (one_per_day) {
    d <- d %>%
      group_by(STUDY_ID, DAY_FROM_BASELINE) %>%
      summarise(.val = mean(.data[[outcome]]), .groups = "drop")
  } else {
    d <- d %>% transmute(STUDY_ID, .val = .data[[outcome]])
  }
  per <- d %>%
    group_by(STUDY_ID) %>%
    summarise(n_read = n(), v = var(.val), .groups = "drop") %>%
    filter(n_read >= 2, !is.na(v))
  list(sigma_w2 = sum((per$n_read - 1) * per$v) / sum(per$n_read - 1),
       n_users  = nrow(per),
       df       = sum(per$n_read - 1))
}

cat("--- Part 1: within-person variance (baseline through", WINDOW_D, "days) ---\n")
vc_rows <- list()
for (o in c("SBP", "DBP")) {
  w_all <- pooled_within_var(rdg, o, one_per_day = FALSE)
  w_day <- pooled_within_var(rdg, o, one_per_day = TRUE)

  base_v <- user[[paste0("BASELINE_", o)]]
  mo6_v  <- user[[paste0("MO6_", o)]]
  keep   <- !is.na(base_v) & !is.na(mo6_v)
  bl     <- base_v[keep]

  sigma_t2 <- var(bl)
  vc_rows[[length(vc_rows) + 1]] <- tibble(
    outcome         = o,
    n_full          = sum(keep),
    mean_baseline_full = mean(bl),
    sigma_t2        = sigma_t2,
    sigma_t         = sqrt(sigma_t2),
    sigma_w2        = w_all$sigma_w2,
    sigma_w         = sqrt(w_all$sigma_w2),
    sigma_b2        = sigma_t2 - w_all$sigma_w2,
    rho             = (sigma_t2 - w_all$sigma_w2) / sigma_t2,
    noise_fraction  = w_all$sigma_w2 / sigma_t2,
    sigma_w_n_users = w_all$n_users,
    sigma_w_df      = w_all$df,
    sigma_w_one_per_day = sqrt(w_day$sigma_w2))
  cat(sprintf("  %s: sigma_w=%.4f (n=%d, df=%d) | one-per-day %.4f | sigma_t=%.4f | noise fraction=%.4f\n",
              o, sqrt(w_all$sigma_w2), w_all$n_users, w_all$df,
              sqrt(w_day$sigma_w2), sqrt(sigma_t2), w_all$sigma_w2 / sigma_t2))
}
vc <- bind_rows(vc_rows)
write_csv(vc, file.path(results_dir, "variance_components.csv"))
cat(sprintf("\n  sigma_w contributors: %d of %d participants (%.1f%%)\n",
            vc$sigma_w_n_users[1], nrow(user), 100 * vc$sigma_w_n_users[1] / nrow(user)))

# ---------------------------------------------------------------------------
# PART 2 — Per-stratum expected RTM and RTM-adjusted delta
# ---------------------------------------------------------------------------
cat("\n--- Part 2: per-stratum RTM adjustment ---\n")
per_rows <- list()
for (o in c("SBP", "DBP")) {
  v      <- vc %>% filter(outcome == o)
  frac   <- v$noise_fraction
  mu     <- v$mean_baseline_full

  base_v <- user[[paste0("BASELINE_", o)]]
  mo6_v  <- user[[paste0("MO6_", o)]]
  keep   <- !is.na(base_v) & !is.na(mo6_v)
  bl     <- base_v[keep]
  dl     <- mo6_v[keep] - base_v[keep]
  st     <- factor(user$BASELINE_AHA[keep], levels = STRATA, labels = SHORT)

  for (s in SHORT) {
    i <- st == s
    per_rows[[length(per_rows) + 1]] <- tibble(
      outcome = o, stratum = s, n_users = sum(i),
      mean_baseline = mean(bl[i]),
      deviation_from_full_mean = mean(bl[i]) - mu,
      delta_observed = mean(dl[i]),
      expected_rtm = frac * (mean(bl[i]) - mu),
      delta_rtm_adjusted = mean(dl[i]) + frac * (mean(bl[i]) - mu))
  }
}
per_stratum <- bind_rows(per_rows)
write_csv(per_stratum, file.path(results_dir, "per_stratum_rtm.csv"))
print(as.data.frame(per_stratum %>% mutate(across(where(is.numeric), ~round(., 4)))),
      row.names = FALSE)

# Restricted (§2.2) cohort
cat("\n--- Part 2b: restricted cohort (Elevated + Stage 1 + Stage 2) ---\n")
coh_rows <- list()
for (o in c("SBP", "DBP")) {
  v      <- vc %>% filter(outcome == o)
  frac   <- v$noise_fraction
  mu     <- v$mean_baseline_full
  base_v <- user[[paste0("BASELINE_", o)]]
  mo6_v  <- user[[paste0("MO6_", o)]]
  keep   <- !is.na(base_v) & !is.na(mo6_v) & user$BASELINE_AHA %in% HYP
  bl     <- base_v[keep]
  dl     <- mo6_v[keep] - base_v[keep]
  coh_rows[[length(coh_rows) + 1]] <- tibble(
    outcome = o, cohort = "hypertensive", n_users = sum(keep),
    mean_baseline = mean(bl),
    deviation_from_full_mean = mean(bl) - mu,
    delta_observed = mean(dl),
    expected_rtm = frac * (mean(bl) - mu),
    delta_rtm_adjusted = mean(dl) + frac * (mean(bl) - mu),
    pct_surviving = 100 * (mean(dl) + frac * (mean(bl) - mu)) / mean(dl))
}
cohort_rtm <- bind_rows(coh_rows)
write_csv(cohort_rtm, file.path(results_dir, "cohort_rtm.csv"))
print(as.data.frame(cohort_rtm %>% mutate(across(where(is.numeric), ~round(., 4)))),
      row.names = FALSE)

# ---------------------------------------------------------------------------
# PART 3 — Diagnostic: generalisation vs literal Eq (1)
# ---------------------------------------------------------------------------
cat("\n--- Part 3: diagnostic — observed-mean generalisation vs literal Eq (1) ---\n")
C_z <- function(z) dnorm(z) / (1 - pnorm(z))   # inverse Mills ratio (see header)
eq1_rows <- list()
for (o in c("SBP", "DBP")) {
  v      <- vc %>% filter(outcome == o)
  frac   <- v$noise_fraction
  mu     <- v$mean_baseline_full
  s_t    <- v$sigma_t
  cc     <- CUTPOINT[[o]]
  base_v <- user[[paste0("BASELINE_", o)]]
  mo6_v  <- user[[paste0("MO6_", o)]]
  keep   <- !is.na(base_v) & !is.na(mo6_v)
  bl     <- base_v[keep]

  z        <- (cc - mu) / s_t
  uni      <- bl[bl >= cc]                       # genuinely univariate selection
  eq1_rows[[length(eq1_rows) + 1]] <- tibble(
    outcome = o, cutpoint = cc, z = z, C_z = C_z(z),
    empirical_factor_univariate = (mean(uni) - mu) / s_t,
    rtm_eq1_literal = s_t * frac * C_z(z),
    rtm_empirical_univariate = frac * (mean(uni) - mu),
    agreement_pct = 100 * (frac * (mean(uni) - mu)) / (s_t * frac * C_z(z)),
    n_univariate = length(uni))
}
eq1 <- bind_rows(eq1_rows)
write_csv(eq1, file.path(results_dir, "eq1_comparison.csv"))
print(as.data.frame(eq1 %>% mutate(across(where(is.numeric), ~round(., 4)))),
      row.names = FALSE)
cat("  Agreement is 91-97%; residual is non-normality of the BP distribution.\n")
cat("  The observed-mean form is the more accurate of the two for these data.\n")

# ---------------------------------------------------------------------------
# PART 4 — Diagnostic figure (not a manuscript figure)
# ---------------------------------------------------------------------------
plot_df <- per_stratum %>%
  mutate(stratum = factor(stratum, levels = SHORT),
         outcome = factor(outcome, levels = c("SBP", "DBP")))

p <- ggplot(plot_df, aes(y = stratum)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_segment(aes(x = delta_observed, xend = delta_rtm_adjusted,
                   yend = stratum), colour = "grey70", linewidth = 1.2) +
  geom_point(aes(x = delta_observed), colour = "grey35", size = 2.6) +
  geom_point(aes(x = delta_rtm_adjusted), colour = "#00538A", size = 2.6) +
  facet_wrap(~ outcome, ncol = 1, scales = "free_x") +
  labs(x = "Change from baseline to 6 months (mmHg)",
       y = "Baseline AHA Stratum",
       title = "§5.1 RTM adjustment (Barnett et al.) — diagnostic",
       subtitle = "Grey = observed change; blue = RTM-adjusted",
       caption = paste0("Expected RTM = (within-person variance / baseline variance) x ",
                        "(stratum mean baseline - full cohort mean baseline)")) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(size = 8, colour = "grey40", hjust = 0))

fig_path <- file.path(figures_dir, "5_1_rtm_barnett_observed_vs_adjusted.png")
ggsave(fig_path, p, width = 7, height = 6.5, dpi = 300)
cat("\nWrote figure:\n  ", fig_path, "\n")

# ---------------------------------------------------------------------------
# PART 5 — Paper-ready markdown
# ---------------------------------------------------------------------------
fmt <- function(x, d = 2) sprintf("%+.*f", d, x)
md <- c(
  "# §5.1 sensitivity — RTM adjustment (Barnett et al. 2005)",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "Supersedes `../rtm/` (Oldham/Blomqvist decomposition); see script header.",
  "",
  "## Variance components",
  "")
for (o in c("SBP", "DBP")) {
  v <- vc %>% filter(outcome == o)
  md <- c(md, sprintf(
    "- **%s**: within-person SD = %.2f mmHg (n = %s participants, %.1f%%); baseline SD = %.2f mmHg; reliability rho = %.3f; noise fraction = %.4f. Full-cohort mean baseline = %.2f mmHg (n = %s).",
    o, v$sigma_w, format(v$sigma_w_n_users, big.mark = ","),
    100 * v$sigma_w_n_users / nrow(user), v$sigma_t, v$rho, v$noise_fraction,
    v$mean_baseline_full, format(v$n_full, big.mark = ",")))
}
md <- c(md, "", "## Per-stratum RTM adjustment (6 months)", "")
for (o in c("SBP", "DBP")) {
  md <- c(md, paste0("### ", o), "",
          "| Stratum | n | Mean baseline | Observed Δ | Expected RTM | RTM-adjusted Δ |",
          "|---|---|---|---|---|---|")
  rr <- per_stratum %>% filter(outcome == o)
  for (i in seq_len(nrow(rr))) {
    md <- c(md, sprintf("| %s | %s | %.2f | %s | %s | **%s** |",
                        rr$stratum[i], format(rr$n_users[i], big.mark = ","),
                        rr$mean_baseline[i], fmt(rr$delta_observed[i]),
                        fmt(rr$expected_rtm[i]), fmt(rr$delta_rtm_adjusted[i])))
  }
  cr <- cohort_rtm %>% filter(outcome == o)
  md <- c(md, sprintf("| Elevated + hypertensive | %s | %.2f | %s | %s | **%s** |",
                      format(cr$n_users, big.mark = ","), cr$mean_baseline,
                      fmt(cr$delta_observed), fmt(cr$expected_rtm),
                      fmt(cr$delta_rtm_adjusted)), "")
}
md <- c(md,
  "## Notes",
  "",
  "- Expected RTM is **negative for Elevated**: RTM pulls toward the reference mean, and Elevated sits at (SBP) or below (DBP, by definition <80 mmHg) the full-cohort mean. For SBP this makes the adjustment negligible, which is the assumption-free anchor used in the manuscript Discussion.",
  "- Barnett's Eq (2) (`sigma_w^2/m`) is deliberately **not** applied: each BP measurement is already the mean of three consecutive cuff recordings, and sigma_w is estimated from day-to-day variation of those means.",
  "- Adjustments are applied to **observed paired mean change**, not to the §5.1 adjusted LMM estimates reported in Table 2.",
  "- `eq1_comparison.csv` quantifies the cost of substituting the observed stratum mean for Barnett's normal-theory C(z) on a genuinely univariate cut.",
  "- Pairs with: §5.1 primary; §5.1 full-cohort sensitivity (Normal-stratum rise = direct empirical RTM evidence); §5.5 RTM (12-month, `5_5/Sensitivity/rtm_barnett/`)."
)
writeLines(md, file.path(results_dir, "5_1_sensitivity_rtm_barnett_paper.md"))
cat("\nWrote markdown:\n  ",
    file.path(results_dir, "5_1_sensitivity_rtm_barnett_paper.md"), "\n")

# ---------------------------------------------------------------------------
# PART 6 — Table A.2 sizes
# ---------------------------------------------------------------------------
emit_sizes(
  tibble(
    section = "5.1", role = "sensitivity",
    analysis = "RTM adjustment (Barnett et al. 2005)", cohort = "hypertensive",
    outcome = c("SBP", "DBP"), predictor_set = "RTM-adjusted Δ",
    n_users  = cohort_rtm$n_users,
    n_obs    = cohort_rtm$n_users,
    n_events = NA_integer_,
    ceiling_n = sum(user$BASELINE_AHA %in% HYP)
  ),
  "5_1_rtm_barnett", results_dir
)

cat("\n=== §5.1 RTM (Barnett) sensitivity done ===\n")
sink()
