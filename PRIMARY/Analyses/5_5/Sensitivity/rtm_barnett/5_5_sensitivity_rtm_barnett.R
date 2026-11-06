# 5_5_sensitivity_rtm_barnett.R
# SAP §5.5 sensitivity: regression-to-the-mean (RTM) adjustment at 12 months,
# Barnett method. Parallels 5_1/Sensitivity/rtm_barnett/.
#
# Reference: Barnett AG, van der Pols JC, Dobson AJ. "Regression to the mean:
#            what it is and how to deal with it." Int J Epidemiol 2005;34:215-220.
#
# SUPERSEDES ../rtm/5_5_sensitivity_rtm.R, whose approach was noticed to be
# incomplete (it applied the same formula as the §5.1 version it was based on).
# Its outputs are retained for provenance but are NOT the manuscript values.
#
# ---------------------------------------------------------------------------
# METHOD — identical to §5.1, applied to the 12-month endpoint:
#
#   expected RTM = (sigma_w^2 / sigma_t^2) * (xbar_stratum - xbar_full)
#   RTM-adjusted Delta = observed Delta + expected RTM
#
# ONE DELIBERATE DIFFERENCE FROM §5.1 — the reference cohort.
# Both sigma_t^2 and xbar_full are computed within the MO12-COMPLETER subset
# (all four AHA strata, complete baseline+MO12), not the full 6-month cohort.
# Rationale: mu must be the mean of the population the restricted group was
# selected from, and the 12-month analysis population IS the completers. Using
# the 6-month cohort mean would mix populations. Consequence: the noise fraction
# is LARGER at 12 months (baseline variance is narrower among completers), so the
# same sigma_w implies a larger proportional correction.
#
# Each BP measurement is the mean of three consecutive cuff recordings and
# sigma_w is estimated from day-to-day variation OF those means, so Barnett's
# Eq (2) (sigma_w^2/m) is NOT applied -- it would double-count the averaging.
#
# CAVEAT to carry into any text: MO12 completers are ~7% of the cohort and are
# self-selected. These estimates describe the completers.
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#         Build/derived/qc/bp_reading_analytic_frame_1X.rds
#
# Outputs (Analyses/5_5/Sensitivity/rtm_barnett/derived/):
#   results/variance_components.csv
#   results/per_stratum_rtm.csv
#   results/cohort_rtm.csv
#   results/5_5_sensitivity_rtm_barnett_paper.md
#   results/5_5_sensitivity_rtm_barnett.log
#   figures/5_5_rtm_barnett_observed_vs_adjusted.png   (diagnostic)

library(tidyverse)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")   # emit_sizes()

project_dir <- "./PRIMARY"
user_path   <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
rdg_path    <- file.path(project_dir, "Build", "derived", "qc", "bp_reading_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "5_5", "Sensitivity", "rtm_barnett", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "5_5_sensitivity_rtm_barnett.log")
sink(log_path, split = TRUE)
cat("§5.5 sensitivity — 12-month RTM adjustment (Barnett et al. 2005) —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

STRATA   <- c("Normal", "Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")
SHORT    <- c("Normal", "Elevated", "Stage 1", "Stage 2")
HYP      <- STRATA[2:4]
WINDOW_D <- 30

user <- read_rds(user_path) %>% filter(SEX %in% c("male", "female"))
rdg  <- read_rds(rdg_path)  %>% semi_join(user %>% select(STUDY_ID), by = "STUDY_ID")

n_mo12 <- sum(!is.na(user$BASELINE_SBP) & !is.na(user$MO12_SBP))
n_mo12_hyp <- sum(!is.na(user$BASELINE_SBP) & !is.na(user$MO12_SBP) &
                    user$BASELINE_AHA %in% HYP)
cat("Full modelled cohort:          ", nrow(user), "users\n")
cat("MO12 completers (all strata):  ", n_mo12,
    sprintf("(%.1f%%)\n", 100 * n_mo12 / nrow(user)))
cat("MO12 completers (hypertensive):", n_mo12_hyp, "\n\n")

# ---------------------------------------------------------------------------
# PART 1 — Within-person variance (same estimator as §5.1)
# ---------------------------------------------------------------------------
# Standalone vectors are built before any tibble()/summarise() call: both
# evaluate arguments sequentially, so a new column name silently shadows an
# outer vector of the same name in later arguments.

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

cat("--- Part 1: variance components (reference cohort = MO12 completers) ---\n")
vc_rows <- list()
for (o in c("SBP", "DBP")) {
  w_all <- pooled_within_var(rdg, o, one_per_day = FALSE)
  w_day <- pooled_within_var(rdg, o, one_per_day = TRUE)

  base_v <- user[[paste0("BASELINE_", o)]]
  mo12_v <- user[[paste0("MO12_", o)]]
  keep   <- !is.na(base_v) & !is.na(mo12_v)
  bl     <- base_v[keep]

  sigma_t2 <- var(bl)
  vc_rows[[length(vc_rows) + 1]] <- tibble(
    outcome            = o,
    n_reference        = sum(keep),
    mean_baseline_full = mean(bl),
    sigma_t2           = sigma_t2,
    sigma_t            = sqrt(sigma_t2),
    sigma_w2           = w_all$sigma_w2,
    sigma_w            = sqrt(w_all$sigma_w2),
    sigma_b2           = sigma_t2 - w_all$sigma_w2,
    rho                = (sigma_t2 - w_all$sigma_w2) / sigma_t2,
    noise_fraction     = w_all$sigma_w2 / sigma_t2,
    sigma_w_n_users    = w_all$n_users,
    sigma_w_df         = w_all$df,
    sigma_w_one_per_day = sqrt(w_day$sigma_w2))
  cat(sprintf("  %s: sigma_w=%.4f | sigma_t=%.4f (MO12 completers) | noise fraction=%.4f | mean baseline=%.2f\n",
              o, sqrt(w_all$sigma_w2), sqrt(sigma_t2),
              w_all$sigma_w2 / sigma_t2, mean(bl)))
}
vc <- bind_rows(vc_rows)
write_csv(vc, file.path(results_dir, "variance_components.csv"))

# ---------------------------------------------------------------------------
# PART 2 — Per-stratum and restricted-cohort RTM adjustment at 12 months
# ---------------------------------------------------------------------------
cat("\n--- Part 2: per-stratum RTM adjustment (12 months) ---\n")
per_rows <- list()
for (o in c("SBP", "DBP")) {
  v      <- vc %>% filter(outcome == o)
  frac   <- v$noise_fraction
  mu     <- v$mean_baseline_full

  base_v <- user[[paste0("BASELINE_", o)]]
  mo12_v <- user[[paste0("MO12_", o)]]
  keep   <- !is.na(base_v) & !is.na(mo12_v)
  bl     <- base_v[keep]
  dl     <- mo12_v[keep] - base_v[keep]
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

cat("\n--- Part 2b: restricted cohort (Elevated + Stage 1 + Stage 2) ---\n")
coh_rows <- list()
for (o in c("SBP", "DBP")) {
  v      <- vc %>% filter(outcome == o)
  frac   <- v$noise_fraction
  mu     <- v$mean_baseline_full
  base_v <- user[[paste0("BASELINE_", o)]]
  mo12_v <- user[[paste0("MO12_", o)]]
  keep   <- !is.na(base_v) & !is.na(mo12_v) & user$BASELINE_AHA %in% HYP
  bl     <- base_v[keep]
  dl     <- mo12_v[keep] - base_v[keep]
  coh_rows[[length(coh_rows) + 1]] <- tibble(
    outcome = o, cohort = "hypertensive MO12 completers", n_users = sum(keep),
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
# PART 3 — Diagnostic figure
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
  labs(x = "Change from baseline to 12 months (mmHg)",
       y = "Baseline AHA Stratum",
       title = "§5.5 RTM adjustment at 12 months (Barnett et al.) — diagnostic",
       subtitle = "Grey = observed change; blue = RTM-adjusted",
       caption = "Reference cohort is the 12-month completer subset; estimates describe completers.") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(size = 8, colour = "grey40", hjust = 0))

fig_path <- file.path(figures_dir, "5_5_rtm_barnett_observed_vs_adjusted.png")
ggsave(fig_path, p, width = 7, height = 6.5, dpi = 300)
cat("\nWrote figure:\n  ", fig_path, "\n")

# ---------------------------------------------------------------------------
# PART 4 — Paper-ready markdown
# ---------------------------------------------------------------------------
fmt <- function(x, d = 2) sprintf("%+.*f", d, x)
md <- c(
  "# §5.5 sensitivity — 12-month RTM adjustment (Barnett et al. 2005)",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "Supersedes `../rtm/` (Oldham/Blomqvist decomposition); see script header.",
  "Parallels `5_1/Sensitivity/rtm_barnett/`, with the reference cohort restricted",
  "to 12-month completers.",
  "",
  "## Variance components",
  "")
for (o in c("SBP", "DBP")) {
  v <- vc %>% filter(outcome == o)
  md <- c(md, sprintf(
    "- **%s**: within-person SD = %.2f mmHg; baseline SD among completers = %.2f mmHg; reliability rho = %.3f; noise fraction = %.4f. Completer-cohort mean baseline = %.2f mmHg (n = %s).",
    o, v$sigma_w, v$sigma_t, v$rho, v$noise_fraction,
    v$mean_baseline_full, format(v$n_reference, big.mark = ",")))
}
md <- c(md, "", "## Per-stratum RTM adjustment (12 months)", "")
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
  "- **Completer selection**: fit on ~7% of the cohort, self-selected by returning at 12 months. Estimates describe completers and are read alongside the better-powered §5.1 6-month RTM.",
  "- The reference cohort (for both `sigma_t^2` and `xbar_full`) is the MO12 completer subset, so the noise fraction is larger than at 6 months and the proportional correction is correspondingly larger.",
  "- 6- and 12-month estimates are computed on **different participant sets**, so they do not constitute a within-person comparison. Any 'deepened from 6 to 12 months' claim requires restricting the 6-month estimate to MO12 completers.",
  "- Barnett's Eq (2) (`sigma_w^2/m`) is deliberately not applied; see script header.",
  "- Adjustments are applied to observed paired mean change, not to the §5.5 adjusted LMM estimates.",
  "- Pairs with: §5.5 primary (12-month LMM); §5.1 RTM (`5_1/Sensitivity/rtm_barnett/`)."
)
writeLines(md, file.path(results_dir, "5_5_sensitivity_rtm_barnett_paper.md"))
cat("\nWrote markdown:\n  ",
    file.path(results_dir, "5_5_sensitivity_rtm_barnett_paper.md"), "\n")

# ---------------------------------------------------------------------------
# PART 5 — Table A.2 sizes
# ---------------------------------------------------------------------------
emit_sizes(
  tibble(
    section = "5.5", role = "sensitivity",
    analysis = "12-month RTM adjustment (Barnett et al. 2005)",
    cohort = "hypertensive MO12 completers",
    outcome = c("SBP", "DBP"), predictor_set = "RTM-adjusted Δ",
    n_users  = cohort_rtm$n_users,
    n_obs    = cohort_rtm$n_users,
    n_events = NA_integer_,
    ceiling_n = n_mo12_hyp
  ),
  "5_5_rtm_barnett", results_dir
)

cat("\n=== §5.5 RTM (Barnett) sensitivity done ===\n")
sink()
