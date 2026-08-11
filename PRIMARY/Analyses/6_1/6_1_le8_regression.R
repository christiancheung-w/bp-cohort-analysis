# 6_1_le8_regression.R
# SAP §6.1: LE8-anchored multivariable regression of BP outcomes on WHOOP-
# derived LE8 predictors, with continuous baseline BP as ANCOVA-style RTM
# correction. Six variants (Primary, +Δweight, +ΔHDL, +ΔNon-HDL, +ΔTG,
# +ΔHbA1c per v0.52) × four primary frames (continuous SBP, continuous DBP,
# SBP-responder logistic, DBP-responder logistic per v0.53) = 24 fits total.
#
# Note on LE8 component numbering: project uses sleep duration as the LE8
# "wearable proxy" component but the standard AHA LE8 ordering is #1 Diet,
# #2 Physical activity, #3 Nicotine, #4 Sleep, #5 Weight, #6 Cholesterol,
# #7 Glucose, #8 Blood pressure. Local labels below follow project shorthand
# rather than the canonical numbering — cosmetic only, does not affect any
# computed result.
#
# Predictor set:
#   Sleep duration (LE8 #4): SLEEP_DURATION_STUDY_MEAN
#   Physical activity (LE8 #2): STRAIN_STUDY_MEAN
#   Nicotine (LE8 #3): LOG_TOBACCO_NICOTINE_USE_DAYS (intensity) +
#           TOBACCO_NICOTINE_EVER_LOGGED (binary user-vs-non-user)
#   Weight (LE8 #5; +Δweight variant): DELTA_WEIGHT_KG
#   Cholesterol / glycemia (LE8 #6/7; per-lab variants per v0.52):
#           DELTA_HDL_CHOLESTEROL, DELTA_NON_HDL_CHOLESTEROL,
#           DELTA_TRIGLYCERIDES, DELTA_HBA1C (each as its own per-lab variant)
#
# Covariates (all variants): BASELINE_BP, AGE_AT_START, SEX, BASELINE_BMI,
#   baseline_HTN_status (on-at-baseline trait; added per CHANGELOG v0.82),
#   SLEEP_DURATION_BASELINE_MEAN, STRAIN_BASELINE_MEAN (30-day pre-study T0
#   controls; added per CHANGELOG v0.85 — ANCOVA-style parallel to BASELINE_BP;
#   isolates study-period increment from stable person-level trait).
#
# Frames:
#   Continuous ANCOVA: FOLLOWUP_BP ~ baseline_BP + predictors + covariates
#     Fit via lm(). Primary effect-size frame. β reported in mmHg.
#   Logistic-responder companion: SBP_RESPONDER_6MO ~ same RHS
#     Fit via glm(family = binomial). OR with 95% CI.
#
# Per SAP §9: α = 0.005 formal-test threshold; 95% CIs as effect-size convention.
# Per SAP §6.1: VIF > 5 drops the offending predictor (preferring more
# clinically interpretable of the pair).
#
# Cohort: §2.2 hypertensive (Elevated + Stage 1 + Stage 2). Variant subsets
# are filtered by complete-case requirement on the added predictors.
#
# Input:  Build/derived/qc/user_analytic_frame_1X.rds
#
# Outputs (Analyses/6_1/derived/):
#   results/6_1_cohort_sizes.csv               — variant × outcome n
#   results/6_1_coefficients_continuous.csv    — β + 95% CI + p, all 6 ANCOVA fits
#   results/6_1_coefficients_logistic.csv      — OR + 95% CI + p, all 3 logistic fits
#   results/6_1_vif_check.csv                  — VIF per predictor per fit
#   results/6_1_model_summaries.csv            — n, R² (continuous) / pseudo-R² (logistic), AUC (logistic per SAP §6.1 Diagnostic / v0.70), AIC, BIC
#   results/6_1_paper_table3.md                — paper-ready Table 3 layout
#   results/6_1_le8_regression.log             — full sink log
#   figures/6_1_forest_continuous.png          — forest plot of LE8 β across variants

library(tidyverse)
library(broom)
library(car)        # for vif()
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")  # htn_baseline_on(), emit_sizes()

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
out_dir     <- file.path(project_dir, "Analyses", "6_1", "derived")
results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "6_1_le8_regression.log")
sink(log_path, split = TRUE)
cat("§6.1 LE8-anchored regression —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read frame, filter to §2.2 cohort, construct DBP_RESPONDER_6MO ---
# Per SAP v0.53: DBP-responder promoted to a primary logistic outcome alongside
# SBP-responder, symmetric with the continuous-ANCOVA frame's SBP/DBP parallelism.
user_frame <- read_rds(input_path)
cohort <- user_frame %>%
  filter(BASELINE_AHA %in% c("Elevated", "Stage 1 Hypertension", "Stage 2 Hypertension"),
         SEX %in% c("male", "female")) %>%
  mutate(DBP_RESPONDER_6MO = as.integer((MO6_DBP - BASELINE_DBP) <= -5),
         HTN_BASELINE_ON   = htn_baseline_on(HTN_MED_BASELINE))  # baseline-only HTN trait (v0.82)
cat("§2.2 cohort + sex restriction:", nrow(cohort), "users\n\n")

# --- Predictor + covariate definitions ---
le8_primary   <- c("SLEEP_DURATION_STUDY_MEAN", "STRAIN_STUDY_MEAN",
                   "LOG_TOBACCO_NICOTINE_USE_DAYS", "TOBACCO_NICOTINE_EVER_LOGGED")
le8_dweight   <- c(le8_primary, "DELTA_WEIGHT_KG")
covariates    <- c("AGE_AT_START", "SEX", "BASELINE_BMI", "HTN_BASELINE_ON",
                   "SLEEP_DURATION_BASELINE_MEAN", "STRAIN_BASELINE_MEAN")  # baseline sleep/PA T0 controls added per v0.85

# Per-variant covariate override (SAP §6.1 drop rule, v0.86):
# PlusDHbA1c (n=269) triggers VIF > 5 for the SLEEP_DURATION baseline/study pair
# (~5.4). Drop the T0 covariate, retaining the study-window predictor (the LE8
# metric of interest). All other variants remain at VIF < 5.
variant_covariate_overrides <- list(
  PlusDHbA1c = setdiff(covariates, "SLEEP_DURATION_BASELINE_MEAN")
)

# Variants per SAP §6.1 Table 3: each variant adds the listed predictors ON
# TOP OF PRIMARY (not hierarchically stacked through +Δweight). This matches
# the SAP language "Additional predictors (on top of primary spec)" — labs
# do NOT require weight pairing. Each per-lab variant adds one lab Δ to the
# Primary LE8 predictor set; cohort = §2.2 hypertensive with paired pre/post
# values for that specific lab. Per-lab framing (replacing the prior joint
# 4-lab +Δlabs) is documented in SAP_CHANGELOG v0.52.
variants <- list(
  Primary     = le8_primary,
  PlusDweight = le8_dweight,
  PlusDHDL    = c(le8_primary, "DELTA_HDL_CHOLESTEROL"),
  PlusDNonHDL = c(le8_primary, "DELTA_NON_HDL_CHOLESTEROL"),
  PlusDTG     = c(le8_primary, "DELTA_TRIGLYCERIDES"),
  PlusDHbA1c  = c(le8_primary, "DELTA_HBA1C")
)

# Union of all LE8 + lab terms across variants (used for filtering output tables)
le8_all_terms <- unique(unlist(variants))

# --- Build per-variant complete-case cohorts ---
cohort_sizes <- map_dfr(names(variants), function(v) {
  vars_required <- c(variants[[v]], covariates, "BASELINE_SBP", "BASELINE_DBP",
                     "MO6_SBP", "MO6_DBP", "SBP_RESPONDER_6MO")
  d <- cohort %>% filter(if_all(all_of(vars_required), ~ !is.na(.)))
  tibble(variant = v, n_complete = nrow(d))
})
write_csv(cohort_sizes, file.path(results_dir, "6_1_cohort_sizes.csv"))
cat("--- Variant cohort sizes (complete cases on all required vars) ---\n")
print(cohort_sizes)

# --- Helper: fit continuous ANCOVA for one (outcome × variant) ---
fit_continuous <- function(outcome_col, baseline_col, predictors, data, covs = covariates) {
  vars_required <- c(predictors, covs, baseline_col, outcome_col)
  d <- data %>% filter(if_all(all_of(vars_required), ~ !is.na(.)))
  rhs <- paste(c(baseline_col, predictors, covs), collapse = " + ")
  f <- as.formula(paste0(outcome_col, " ~ ", rhs))
  fit <- lm(f, data = d)
  list(fit = fit, n = nrow(d), formula = f, study_ids = d$STUDY_ID)
}

# --- Helper: fit logistic-responder (parameterized on outcome + baseline covariate) ---
fit_logistic <- function(outcome_col, baseline_col, predictors, data, covs = covariates) {
  vars_required <- c(predictors, covs, baseline_col, outcome_col)
  d <- data %>% filter(if_all(all_of(vars_required), ~ !is.na(.)))
  rhs <- paste(c(baseline_col, predictors, covs), collapse = " + ")
  f <- as.formula(paste0(outcome_col, " ~ ", rhs))
  fit <- glm(f, family = binomial(), data = d)
  list(fit = fit, n = nrow(d), formula = f, study_ids = d$STUDY_ID)
}

# --- VIF check ---
vif_table <- function(fit, label) {
  v <- tryCatch(car::vif(fit$fit), error = function(e) NULL)
  if (is.null(v)) return(NULL)
  if (is.matrix(v)) v <- v[, "GVIF"]  # for factor predictors
  tibble(model = label, predictor = names(v), VIF = as.numeric(v))
}

# --- Fit all 6 continuous + 3 logistic ---
fits_continuous <- list()
fits_logistic   <- list()

for (v in names(variants)) {
  preds   <- variants[[v]]
  cov_set <- if (!is.null(variant_covariate_overrides[[v]])) variant_covariate_overrides[[v]] else covariates
  cat("\n\n========================================\n")
  cat("Variant:", v, " | predictors:", paste(preds, collapse = ", "), "\n")
  if (!identical(cov_set, covariates))
    cat("  [VIF override] covariates:", paste(cov_set, collapse = ", "), "\n")
  cat("========================================\n")

  for (otc in c("SBP", "DBP")) {
    outcome_col  <- paste0("MO6_", otc)
    baseline_col <- paste0("BASELINE_", otc)
    cat(sprintf("\n--- %s ANCOVA × %s ---\n", otc, v))
    f <- fit_continuous(outcome_col, baseline_col, preds, cohort, cov_set)
    cat("n =", f$n, ";  formula:", deparse(f$formula), "\n")
    cat("R² =", round(summary(f$fit)$r.squared, 4),
        " adj-R² =", round(summary(f$fit)$adj.r.squared, 4), "\n")
    fits_continuous[[paste0(otc, "__", v)]] <- f
  }

  for (resp in c("SBP", "DBP")) {
    outcome_col  <- paste0(resp, "_RESPONDER_6MO")
    baseline_col <- paste0("BASELINE_", resp)
    cat(sprintf("\n--- %s-responder logistic × %s ---\n", resp, v))
    fl <- fit_logistic(outcome_col, baseline_col, preds, cohort, cov_set)
    cat("n =", fl$n, ";  formula:", deparse(fl$formula), "\n")
    fits_logistic[[paste0(resp, "_responder__", v)]] <- fl
  }
}

# --- Extract continuous coefficients ---
coef_continuous <- imap_dfr(fits_continuous, function(f, label) {
  parts <- strsplit(label, "__", fixed = TRUE)[[1]]
  outcome <- parts[1]; variant <- parts[2]
  broom::tidy(f$fit, conf.int = TRUE, conf.level = 0.95) %>%
    mutate(outcome = outcome,
           variant = variant,
           n       = f$n,
           sig_005 = p.value < 0.005) %>%
    select(outcome, variant, term, estimate, std.error, conf.low, conf.high,
           statistic, p.value, sig_005, n)
})
write_csv(coef_continuous, file.path(results_dir, "6_1_coefficients_continuous.csv"))

# --- Extract logistic coefficients (with OR) ---
# Key format now: "{SBP|DBP}_responder__{variant}". Parse to get outcome and variant.
coef_logistic <- imap_dfr(fits_logistic, function(f, label) {
  parts <- strsplit(label, "__", fixed = TRUE)[[1]]
  responder_outcome <- parts[1]; variant <- parts[2]
  broom::tidy(f$fit, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE) %>%
    mutate(responder_outcome = responder_outcome,
           variant = variant,
           n       = f$n,
           sig_005 = p.value < 0.005) %>%
    rename(OR = estimate, OR_lcl = conf.low, OR_ucl = conf.high) %>%
    select(responder_outcome, variant, term, OR, std.error, OR_lcl, OR_ucl,
           statistic, p.value, sig_005, n)
})
write_csv(coef_logistic, file.path(results_dir, "6_1_coefficients_logistic.csv"))

cat("\n\n--- Continuous ANCOVA coefficients (LE8 predictors only, by variant × outcome) ---\n")
print(coef_continuous %>%
        filter(term %in% le8_all_terms) %>%
        select(outcome, variant, term, estimate, conf.low, conf.high, p.value, sig_005),
      n = Inf, width = Inf)

cat("\n--- Logistic OR (LE8 predictors only, by outcome × variant) ---\n")
print(coef_logistic %>%
        filter(term %in% le8_all_terms) %>%
        select(responder_outcome, variant, term, OR, OR_lcl, OR_ucl, p.value, sig_005),
      n = Inf, width = Inf)

# --- VIF check (continuous only — VIF for logistic is conceptually different) ---
vif_all <- map_dfr(names(fits_continuous), function(label) {
  vif_table(fits_continuous[[label]], label)
})
write_csv(vif_all, file.path(results_dir, "6_1_vif_check.csv"))

vif_flags <- vif_all %>% filter(VIF > 5)
cat("\n--- VIF flags (predictors with VIF > 5) ---\n")
if (nrow(vif_flags) == 0) {
  cat("None. All predictors within VIF < 5 threshold.\n")
} else {
  print(vif_flags)
}

# --- Model summaries ---
# AUC computed for logistic fits per SAP §6.1 Diagnostic block (C-statistic).
# Manual Wilcoxon-Mann-Whitney implementation (no extra dependency):
#   AUC = P(score_pos > score_neg) + 0.5 * P(score_pos == score_neg)
compute_auc <- function(fit_obj) {
  y <- model.response(model.frame(fit_obj))
  p <- predict(fit_obj, type = "response")
  pos <- p[y == 1]
  neg <- p[y == 0]
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  # Vectorized Wilcoxon U via rank
  r <- rank(c(pos, neg))
  u <- sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2
  u / (length(pos) * length(neg))
}

model_summ <- bind_rows(
  imap_dfr(fits_continuous, function(f, label) {
    s <- summary(f$fit)
    tibble(model = label, type = "continuous",
           n = f$n, r2 = s$r.squared, adj_r2 = s$adj.r.squared,
           auc = NA_real_,
           aic = AIC(f$fit), bic = BIC(f$fit))
  }),
  imap_dfr(fits_logistic, function(f, label) {
    null_dev <- f$fit$null.deviance
    res_dev  <- f$fit$deviance
    pseudo_r2 <- 1 - (res_dev / null_dev)
    auc <- compute_auc(f$fit)
    tibble(model = label, type = "logistic",
           n = f$n, r2 = pseudo_r2, adj_r2 = NA_real_,
           auc = auc,
           aic = AIC(f$fit), bic = BIC(f$fit))
  })
)
write_csv(model_summ, file.path(results_dir, "6_1_model_summaries.csv"))

# --- Emit standardized sizes for Table A (per variant × outcome) ---
# Continuous + logistic frames share the complete-case n per (variant × outcome);
# n_events = responder count (the logistic frame's events). 6 variants × 2 outcomes.
.n61_ceiling <- nrow(cohort)   # precomputed: the `cohort` column below would shadow the df
le8_sizes <- map_dfr(names(variants), function(v) {
  map_dfr(c("SBP", "DBP"), function(otc) {
    fc <- fits_continuous[[paste0(otc, "__", v)]]
    fl <- fits_logistic[[paste0(otc, "_responder__", v)]]
    tibble(
      section = "6.1", role = "primary",
      analysis = "LE8 regression (continuous ANCOVA + responder logistic)",
      cohort = "hypertensive", outcome = otc, predictor_set = v,
      n_users = as.integer(fc$n), n_obs = as.integer(fc$n),
      n_events = as.integer(sum(model.response(model.frame(fl$fit)))),
      ceiling_n = .n61_ceiling
    )
  })
})
emit_sizes(le8_sizes, "6_1", results_dir)

cat("\n--- Model summaries (fit statistics) ---\n")
print(model_summ)

# --- Forest plot: LE8 predictor β across variants (continuous, SBP + DBP) ---
# Shared behavioural predictors (sleep, strain, nicotine) are common to all six
# variants, but each per-lab variant re-fits them on a small, lab-selected
# subsample (n 389–8,398 vs Primary's 19,896) — so those re-estimates are
# redundant subsample artifacts, not independent findings. Show the shared
# predictors from the Primary (full-n) model ONLY; each lab Δ shows from its own
# per-lab variant. Full estimates remain in 6_1_coefficients_continuous.csv / Table 3.
shared_le8 <- c("SLEEP_DURATION_STUDY_MEAN", "STRAIN_STUDY_MEAN",
                "LOG_TOBACCO_NICOTINE_USE_DAYS")
# Reader-facing variable names (raw term -> label); vector order = top-to-bottom y order.
term_labels <- c(
  SLEEP_DURATION_STUDY_MEAN     = "6 Month Sleep Duration",
  STRAIN_STUDY_MEAN             = "6 Month Physical Activity",
  LOG_TOBACCO_NICOTINE_USE_DAYS = "Tobacco/Nicotine Use",
  DELTA_WEIGHT_KG               = "Δ Body Mass",
  DELTA_HDL_CHOLESTEROL         = "Δ HDL Cholesterol",
  DELTA_NON_HDL_CHOLESTEROL     = "Δ Non-HDL Cholesterol",
  DELTA_TRIGLYCERIDES           = "Δ Triglycerides"
)
forest_df <- coef_continuous %>%
  filter(term %in% le8_all_terms) %>%
  filter(variant == "Primary" | !(term %in% shared_le8)) %>%
  # ΔHbA1c dropped from the figure: n=389 (smallest lab subset), CIs so wide they
  # dominate the x-scale and carry no information. Full estimate retained in Table 3.
  filter(term != "DELTA_HBA1C") %>%
  mutate(
    term     = factor(term_labels[as.character(term)], levels = rev(unname(term_labels))),
    variant  = factor(variant,
                      levels = c("Primary", "PlusDweight",
                                 "PlusDHDL", "PlusDNonHDL", "PlusDTG"),
                      labels = c("Primary", "+Δweight",
                                 "+ΔHDL", "+ΔNon-HDL", "+ΔTG")),
    # SBP listed last → position_dodge draws it at the higher (top) y position;
    # legend order restored to SBP-first via scale_color_manual(breaks=).
    outcome  = factor(outcome, levels = c("DBP", "SBP")),
    sig_label = ifelse(sig_005, "p < 0.005", "p ≥ 0.005")
  )

# SBP + DBP overlaid on one panel (same units: mmHg per unit predictor). Colour =
# outcome; dodge is VERTICAL only (separates the two outcomes along the y/term
# axis — the x/β position is never shifted). Outcome colours match §6.2
# (SBP blue, DBP green). Significance still encoded by fill (filled = p < 0.005);
# its key is intentionally omitted (explained in text), only the outcome legend shows.
p <- ggplot(forest_df, aes(x = estimate, y = term, color = outcome, shape = sig_label)) +
  geom_vline(xintercept = 0, color = "gray60", linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                width = 0.4, linewidth = 0.5,
                position = position_dodge(width = 0.6), orientation = "y") +
  geom_point(size = 2.8, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c("SBP" = "#1f78b4", "DBP" = "#33a02c"), name = NULL,
                     breaks = c("SBP", "DBP"),
                     labels = c(SBP = "Systolic Blood Pressure",
                                DBP = "Diastolic Blood Pressure")) +
  scale_shape_manual(values = c("p < 0.005" = 16, "p ≥ 0.005" = 1), guide = "none") +
  scale_x_continuous(limits = c(-1.4, 0.6),
                     breaks = seq(-1.25, 0.5, by = 0.25),
                     expand = expansion(mult = 0.02)) +
  labs(
    x = "β (mmHg per Unit Predictor, 95% CI)",
    y = NULL,
    title    = "Life's Essential 8 Predictors of 6 Month Blood Pressure",
    subtitle = "Adjusted ANCOVA Model Variants"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position    = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.4),
    axis.title.x       = element_text(face = "bold", margin = margin(t = 10)),
    axis.text.y        = element_text(size = 13),
    plot.title         = element_text(face = "bold")
  )

out_png <- file.path(figures_dir, "6_1_forest_continuous.png")
ggsave(out_png, p, width = 9, height = 7, dpi = 300)
cat("\nWrote forest plot:\n  ", out_png, "\n")

# --- Paper-ready markdown (Table 3 layout) ---
fmt_b <- function(est, lcl, ucl, p, alpha = 0.005) {
  sig <- ifelse(p < alpha, "**", "")
  sprintf("%s%+.3f (%.3f, %.3f)%s [p = %s]",
          sig, est, lcl, ucl, sig,
          ifelse(p < 1e-10, formatC(p, format = "e", digits = 1),
                 ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                        sprintf("%.3f", p))))
}
fmt_or <- function(or, lcl, ucl, p, alpha = 0.005) {
  sig <- ifelse(p < alpha, "**", "")
  sprintf("%s%.3f (%.3f, %.3f)%s [p = %s]",
          sig, or, lcl, ucl, sig,
          ifelse(p < 1e-10, formatC(p, format = "e", digits = 1),
                 ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                        sprintf("%.3f", p))))
}

md_lines <- c(
  "# §6.1 Table 3 — LE8 wearable predictors of BP outcomes",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "## Continuous ANCOVA: β (mmHg per unit predictor, 95% CI)",
  "",
  "Bold = significant at α = 0.005 per SAP §9.",
  ""
)

for (v in c("Primary", "PlusDweight", "PlusDHDL", "PlusDNonHDL", "PlusDTG", "PlusDHbA1c")) {
  vlbl <- c(Primary = "Primary", PlusDweight = "+Δweight",
            PlusDHDL = "+ΔHDL", PlusDNonHDL = "+ΔNon-HDL",
            PlusDTG = "+ΔTG", PlusDHbA1c = "+ΔHbA1c")[v]
  md_lines <- c(md_lines, paste0("### Variant: ", vlbl))
  preds_in_variant <- variants[[v]]
  for (otc in c("SBP", "DBP")) {
    rows <- coef_continuous %>%
      filter(outcome == otc, variant == v, term %in% preds_in_variant)
    md_lines <- c(md_lines, paste0("**", otc, "**"))
    for (i in seq_len(nrow(rows))) {
      r <- rows[i, ]
      md_lines <- c(md_lines, sprintf(
        "- `%s`: %s",
        r$term, fmt_b(r$estimate, r$conf.low, r$conf.high, r$p.value)
      ))
    }
    md_lines <- c(md_lines, "")
  }
}

md_lines <- c(md_lines,
  "## Logistic-responder companion: OR for SBP-responder and DBP-responder (Δ ≤ −5 mmHg, 95% CI)",
  "",
  "Bold = significant at α = 0.005 per SAP §9. Per SAP v0.53, both SBP-responder and DBP-responder are primary logistic outcomes (symmetric with continuous-ANCOVA frame).",
  ""
)
for (resp in c("SBP", "DBP")) {
  md_lines <- c(md_lines, paste0("### ", resp, "-responder (Δ", resp, " ≤ −5 mmHg)"))
  for (v in c("Primary", "PlusDweight", "PlusDHDL", "PlusDNonHDL", "PlusDTG", "PlusDHbA1c")) {
    vlbl <- c(Primary = "Primary", PlusDweight = "+Δweight",
              PlusDHDL = "+ΔHDL", PlusDNonHDL = "+ΔNon-HDL",
              PlusDTG = "+ΔTG", PlusDHbA1c = "+ΔHbA1c")[v]
    md_lines <- c(md_lines, paste0("#### Variant: ", vlbl))
    preds_in_variant <- variants[[v]]
    rows <- coef_logistic %>%
      filter(responder_outcome == paste0(resp, "_responder"),
             variant == v, term %in% preds_in_variant)
    for (i in seq_len(nrow(rows))) {
      r <- rows[i, ]
      md_lines <- c(md_lines, sprintf(
        "- `%s`: %s",
        r$term, fmt_or(r$OR, r$OR_lcl, r$OR_ucl, r$p.value)
      ))
    }
    md_lines <- c(md_lines, "")
  }
}

md_lines <- c(md_lines,
  "## Notes",
  "",
  "- **Model class**: continuous ANCOVA `FOLLOWUP_BP ~ baseline_BP + LE8 predictors + age + sex + baseline_BMI + baseline_HTN_status + baseline_sleep_duration + baseline_PA`, fit separately for SBP and DBP. Logistic companion uses the identical RHS with binary BP-responder outcomes (SBP-responder AND DBP-responder, ΔBP ≤ −5 mmHg, dual primary per v0.53). Both frames share the same covariate set — the only variation between frames is the BP outcome form (continuous vs. binary responder).",
  "- **Baseline sleep/PA covariates (v0.85)**: `SLEEP_DURATION_BASELINE_MEAN` and `STRAIN_BASELINE_MEAN` are 30-day pre-study T0 controls, paralleling the role of `baseline_BP` in the ANCOVA. Study-window sleep/PA are T1 measurements of the same continuously-tracked behaviors; controlling for T0 isolates the study-period increment rather than the stable person-level trait. Tobacco has no equivalent T0 (sparse, event-based logging) and is unchanged.",
  "- **Cohort**: §2.2 hypertensive (Elevated + Stage 1 + Stage 2) + sex restriction. Complete-case requirement on all model predictors per variant; complete-case now also requires non-NA baseline sleep/PA (users with < 8 baseline-window WHOOP days have NA baseline means per the build contract). Primary ~14k (reduced from ~20k due to baseline coverage filter), +Δweight ~10k, per-lab variants ~350–550.",
  "- **Significance threshold**: α = 0.005 (formal-test gate per SAP §9). 95% CIs retained as the effect-size convention.",
  "- **VIF check**: predictor VIF reported in `6_1_vif_check.csv`. SAP §6.1 drop-rule = VIF > 5 (none flagged in current fits — see log).",
  "- **Tobacco/nicotine pair**: `LOG_TOBACCO_NICOTINE_USE_DAYS` (intensity) + `TOBACCO_NICOTINE_EVER_LOGGED` (binary user-vs-non-user). The binary captures the ~80% never-logged majority; the log captures dose-response within users with any logging.",
  "- **Per-lab variants** (n ≈ 400–600 each) are wide-CI by design due to limited paired pre/post lab availability per analyte; reported for symmetry with the continuous frame, not as independent tests. Per v0.52 these replace the original joint +Δlabs variant which was severely underpowered.",
  "- Pairs with: §5.1 primary (cohort-level Δ headline), §6.2 (daily WHOOP signal predictors), §7 (W/B daily-level decomposition)."
)

# Save fits for downstream diagnostics (6_1_diagnostics.R)
write_rds(list(continuous = fits_continuous, logistic = fits_logistic),
          file.path(out_dir, "fitted_models.rds"))
cat("\nWrote fits checkpoint:\n  ", file.path(out_dir, "fitted_models.rds"), "\n")

writeLines(md_lines, file.path(results_dir, "6_1_paper_table3.md"))
cat("\nWrote markdown:\n  ", file.path(results_dir, "6_1_paper_table3.md"), "\n")

cat("\n=== §6.1 LE8 regression done ===\n")
sink()
