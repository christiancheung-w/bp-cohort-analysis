# rtm_barnett_bootstrap.R
# RTM-corrected BP change at 6 and 12 months, with participant-bootstrap 95% CIs.
#
#   RTM-corrected change = adjusted LMM change from baseline (per stratum)
#                          + expected RTM
#   expected RTM         = (sigma_w^2 / sigma_t^2) * (stratum mean baseline - reference mean baseline)
#
# Expected RTM follows ../rtm_barnett/ (Barnett et al. 2005, observed stratum mean
# deviation in place of C(z)); here it is added to the adjusted model estimate
# instead of the raw mean change.
#
# Horizons:
#   MO6  - adjusted LMM as 5_1_primary_lmm.R; reference = all four strata with baseline + MO6
#   MO12 - adjusted LMM as 5_5_primary_lmm.R; reference = all four strata with baseline + MO12
#   sigma_w^2 (both) = pooled within-participant variance of readings on days 0-30,
#                      all modelled participants
#
# Bootstrap: participants resampled with replacement within baseline AHA stratum.
# Each resample refits both LMMs and recomputes sigma_w^2, sigma_t^2 and all means.
# CIs = 2.5th / 97.5th percentiles. Per-resample seeds are drawn from SEED, so draws
# do not depend on the number of cores.
#
# Usage (project root): Rscript PRIMARY/Analyses/5_1/Sensitivity/rtm_barnett_bootstrap/rtm_barnett_bootstrap.R [B=2000] [cores=16]
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#         Build/derived/qc/bp_reading_analytic_frame_1X.rds
# Check:  the original sample must reproduce the adjusted MO6/MO12 - BASELINE
#         contrasts in 5_1 / 5_5 delta_contrasts.csv and expected_rtm in both
#         rtm_barnett per_stratum_rtm.csv files; the script stops otherwise.
# Outputs (Analyses/5_1/Sensitivity/rtm_barnett_bootstrap/derived/results/):
#   rtm_bootstrap_summary.csv - horizon x outcome x stratum: point estimates + CIs
#   rtm_bootstrap_draws.csv   - every resample
#   rtm_bootstrap.log

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(emmeans)
  library(parallel)
})
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")   # htn_baseline_on()
emm_options(lmer.df = "asymptotic")

args  <- commandArgs(trailingOnly = TRUE)
B     <- if (length(args) >= 1) as.integer(args[1]) else 2000L
CORES <- if (length(args) >= 2) as.integer(args[2]) else 16L
SEED     <- 20261102L
WINDOW_D <- 30
TOL      <- 1e-7

project_dir <- "./PRIMARY"
user_path   <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
rdg_path    <- file.path(project_dir, "Build", "derived", "qc", "bp_reading_analytic_frame_1X.rds")
ref_paths   <- list(
  MO6  = list(delta = file.path(project_dir, "Analyses", "5_1", "derived", "results", "5_1_delta_contrasts.csv"),
              rtm   = file.path(project_dir, "Analyses", "5_1", "Sensitivity", "rtm_barnett", "derived", "results", "per_stratum_rtm.csv")),
  MO12 = list(delta = file.path(project_dir, "Analyses", "5_5", "derived", "results", "5_5_delta_contrasts.csv"),
              rtm   = file.path(project_dir, "Analyses", "5_5", "Sensitivity", "rtm_barnett", "derived", "results", "per_stratum_rtm.csv")))
out_dir     <- file.path(project_dir, "Analyses", "5_1", "Sensitivity", "rtm_barnett_bootstrap", "derived")
results_dir <- file.path(out_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

sink(file.path(results_dir, "rtm_bootstrap.log"), split = TRUE)
cat("RTM-corrected change, participant bootstrap —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("B =", B, "| cores =", CORES, "| seed =", SEED, "\n\n")

STRATA <- c("Normal", "Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension")
SHORT  <- c("Normal", "Elevated", "Stage 1", "Stage 2")
HYP    <- STRATA[2:4]
TIMES  <- list(MO6 = c("BASELINE", "MO3", "MO6"), MO12 = c("BASELINE", "MO3", "MO6", "MO12"))

user <- read_rds(user_path) %>% filter(SEX %in% c("male", "female"))
rdg  <- read_rds(rdg_path)  %>% semi_join(user %>% select(STUDY_ID), by = "STUDY_ID")
cat("Modelled participants (all four strata):", nrow(user), "\n")

# Per-participant within-variance pieces, days 0-30: (n - 1) * var and (n - 1).
# A resample's pooled sigma_w^2 = sum(ss * w) / sum(df * w), w = times drawn.
within_pieces <- function(o) {
  per <- rdg %>%
    filter(DAY_FROM_BASELINE >= 0, DAY_FROM_BASELINE <= WINDOW_D, !is.na(.data[[o]])) %>%
    group_by(STUDY_ID) %>%
    summarise(n_read = n(), v = var(.data[[o]]), .groups = "drop") %>%
    filter(n_read >= 2, !is.na(v))
  m  <- match(user$STUDY_ID, per$STUDY_ID)
  list(ss = ifelse(is.na(m), 0, (per$n_read[m] - 1) * per$v[m]),
       df = ifelse(is.na(m), 0, per$n_read[m] - 1))
}
WP <- list(SBP = within_pieces("SBP"), DBP = within_pieces("DBP"))
rm(rdg); invisible(gc())

expected_rtm <- function(w, o, horizon) {
  base <- user[[paste0("BASELINE_", o)]]
  fu   <- user[[paste0(horizon, "_", o)]]
  keep <- !is.na(base) & !is.na(fu) & w > 0
  bl   <- rep(base[keep], w[keep])
  st   <- factor(rep(user$BASELINE_AHA[keep], w[keep]), levels = STRATA)
  sw2  <- sum(WP[[o]]$ss * w) / sum(WP[[o]]$df * w)
  st2  <- var(bl)
  mu   <- mean(bl)
  m    <- as.numeric(tapply(bl, st, mean))
  tibble(stratum = STRATA, mean_baseline = m, ref_mean_baseline = mu,
         sigma_w2 = sw2, sigma_t2 = st2, noise_fraction = sw2 / st2,
         expected_rtm = (sw2 / st2) * (m - mu))
}

to_long <- function(u, times) {
  cols <- c(outer(times, c("SBP", "DBP", "TIME_OF_DAY", "SEASON"), paste, sep = "_"))
  u %>%
    select(STUDY_ID, BASELINE_AHA, AGE_AT_START, SEX, BASELINE_BMI, HTN_MED_BASELINE, all_of(cols)) %>%
    pivot_longer(all_of(cols), names_to = c("Time", ".value"),
                 names_pattern = paste0("(", paste(times, collapse = "|"), ")_(.+)")) %>%
    mutate(Time            = factor(Time, levels = times),
           BASELINE_AHA    = factor(BASELINE_AHA, levels = HYP),
           SEX             = factor(SEX),
           SEASON          = factor(SEASON, levels = c("winter", "spring", "summer", "fall")),
           HTN_BASELINE_ON = htn_baseline_on(HTN_MED_BASELINE))
}

fit_delta <- function(long, o) {
  f <- as.formula(paste0(
    o, " ~ Time * BASELINE_AHA",
    " + AGE_AT_START + SEX + BASELINE_BMI + HTN_BASELINE_ON",
    " + sin(2 * pi * TIME_OF_DAY / 24) + cos(2 * pi * TIME_OF_DAY / 24)",
    " + SEASON",
    " + Time:AGE_AT_START + Time:SEX + Time:BASELINE_BMI + Time:HTN_BASELINE_ON",
    " + (1 | STUDY_ID)"))
  d <- long %>% filter(!is.na(.data[[o]]), !is.na(BASELINE_AHA),
                       !is.na(AGE_AT_START), !is.na(SEX), !is.na(BASELINE_BMI),
                       !is.na(HTN_BASELINE_ON), !is.na(TIME_OF_DAY), !is.na(SEASON))
  fit <- lmer(f, data = d, REML = TRUE)
  k   <- nlevels(d$Time)
  ct  <- summary(contrast(suppressMessages(emmeans(fit, ~ Time | BASELINE_AHA)),
                          list(delta = c(-1, rep(0, k - 2), 1))), infer = TRUE)
  n   <- d %>% distinct(STUDY_ID, BASELINE_AHA) %>% count(BASELINE_AHA)
  tibble(stratum   = as.character(ct$BASELINE_AHA),
         n_users   = n$n[match(as.character(ct$BASELINE_AHA), as.character(n$BASELINE_AHA))],
         delta     = ct$estimate,
         delta_lcl = ct$asymp.LCL,
         delta_ucl = ct$asymp.UCL,
         converged = length(fit@optinfo$conv$lme4$messages) == 0)
}

one_rep <- function(idx) {
  w  <- tabulate(idx, nbins = nrow(user))
  ub <- user[idx, ]
  if (!identical(idx, seq_len(nrow(user)))) ub$STUDY_ID <- paste0("b", seq_len(nrow(ub)))
  hyp <- ub[ub$BASELINE_AHA %in% HYP, ]
  out <- list()
  for (h in names(TIMES)) {
    long <- to_long(hyp, TIMES[[h]])
    for (o in c("SBP", "DBP")) {
      out[[length(out) + 1]] <- fit_delta(long, o) %>%
        left_join(expected_rtm(w, o, h), by = "stratum") %>%
        mutate(horizon = h, outcome = o, rtm_corrected = delta + expected_rtm)
    }
  }
  bind_rows(out)
}

strata_idx <- split(seq_len(nrow(user)), factor(user$BASELINE_AHA, levels = STRATA))
draw_idx   <- function() unlist(lapply(strata_idx, function(ix) ix[sample.int(length(ix), length(ix), replace = TRUE)]),
                                use.names = FALSE)

# --- Original sample + reproduction check ---
cat("--- Original sample ---\n")
point <- one_rep(seq_len(nrow(user)))
for (h in names(TIMES)) {
  ref_d <- read_csv(ref_paths[[h]]$delta, show_col_types = FALSE) %>%
    filter(contrast == paste(h, "- BASELINE"), model %in% c("sbp_adj", "dbp_adj")) %>%
    transmute(outcome = toupper(sub("_adj", "", model)), stratum = BASELINE_AHA, ref_delta = estimate)
  ref_r <- read_csv(ref_paths[[h]]$rtm, show_col_types = FALSE) %>%
    transmute(outcome, stratum = STRATA[match(stratum, SHORT)], ref_rtm = expected_rtm)
  chk <- point %>% filter(horizon == h) %>%
    inner_join(ref_d, by = c("outcome", "stratum")) %>%
    inner_join(ref_r, by = c("outcome", "stratum"))
  dd <- max(abs(chk$delta - chk$ref_delta)); dr <- max(abs(chk$expected_rtm - chk$ref_rtm))
  cat(sprintf("  %s: %d cells matched | max |delta - frozen| = %.2e | max |RTM - frozen| = %.2e\n",
              h, nrow(chk), dd, dr))
  if (nrow(chk) != 6 || dd > TOL || dr > TOL) stop("Original sample does not reproduce the frozen ", h, " estimates.")
}

# --- Bootstrap ---
cat("\n--- Bootstrap ---\n")
set.seed(SEED)
seeds <- sample.int(.Machine$integer.max, B)
run_b <- function(b) {
  set.seed(seeds[b])
  r <- tryCatch(withCallingHandlers(one_rep(draw_idx()),
                                    warning = function(w) invokeRestart("muffleWarning")),
                error = function(e) NULL)
  if (!is.null(r)) r$b <- b
  r
}
t0    <- Sys.time()
draws <- mclapply(seq_len(B), run_b, mc.cores = CORES, mc.preschedule = TRUE)
n_failed <- sum(vapply(draws, is.null, logical(1)))
draws <- bind_rows(draws)
cat(sprintf("  %d resamples in %.1f min | failed: %d | non-converged fits: %d of %d\n",
            B, as.numeric(difftime(Sys.time(), t0, units = "mins")), n_failed,
            sum(!draws$converged), nrow(draws)))

# --- Summary ---
q <- function(x, p) unname(quantile(x, p, type = 7))
boot <- draws %>%
  group_by(horizon, outcome, stratum) %>%
  summarise(B_ok              = n(),
            rtm_corrected_se  = sd(rtm_corrected),
            rtm_corrected_lcl = q(rtm_corrected, 0.025),
            rtm_corrected_ucl = q(rtm_corrected, 0.975),
            delta_boot_lcl    = q(delta, 0.025),
            delta_boot_ucl    = q(delta, 0.975),
            expected_rtm_lcl  = q(expected_rtm, 0.025),
            expected_rtm_ucl  = q(expected_rtm, 0.975),
            .groups = "drop")
summary_tbl <- point %>%
  select(horizon, outcome, stratum, n_users, delta, delta_lcl, delta_ucl,
         mean_baseline, ref_mean_baseline, sigma_w2, sigma_t2, noise_fraction,
         expected_rtm, rtm_corrected) %>%
  left_join(boot, by = c("horizon", "outcome", "stratum")) %>%
  mutate(horizon = factor(horizon, levels = names(TIMES)),
         outcome = factor(outcome, levels = c("SBP", "DBP")),
         stratum = factor(stratum, levels = HYP)) %>%
  arrange(horizon, outcome, stratum) %>%
  mutate(across(c(horizon, outcome, stratum), as.character))

write_csv(summary_tbl, file.path(results_dir, "rtm_bootstrap_summary.csv"))
write_csv(draws,       file.path(results_dir, "rtm_bootstrap_draws.csv"))

cat("\n")
print(as.data.frame(summary_tbl %>%
        transmute(horizon, outcome, stratum, n_users,
                  adjusted = sprintf("%.2f (%.2f, %.2f)", delta, delta_lcl, delta_ucl),
                  expected_rtm = sprintf("%.2f", expected_rtm),
                  rtm_corrected = sprintf("%.2f (%.2f, %.2f)", rtm_corrected, rtm_corrected_lcl, rtm_corrected_ucl))),
      row.names = FALSE)
cat("\nWrote:", file.path(results_dir, "rtm_bootstrap_summary.csv"), "\n")
cat("       ", file.path(results_dir, "rtm_bootstrap_draws.csv"), "\n")
sink()
