# 6_1_diagnostics.R
# SAP §6.1 diagnostic outputs:
#   Continuous ANCOVA fits (12, lm):    Q-Q residuals + fitted/residual + LOESS overlay
#   Logistic-responder fits (12, glm):  decile-based calibration plot
#   Cook's distance:                    lab variants only (16 fits = 4 lab variants
#                                       × 2 outcomes × 2 frames) — at n = 389–592 the
#                                       4/n threshold is informative; Primary (~20k)
#                                       and +Δweight (~12k) skip Cook's because the
#                                       threshold is uninformative at large n (same
#                                       rationale as v0.58 §6.2/§8 Cook's removal).
#
# Nakagawa R² is not applicable (no random effects in lm/glm). Adjusted R² for
# continuous fits and AUC for logistic fits are already in 6_1_model_summaries.csv.
#
# LOESS on the fitted/residual scatter uses a 5,000-row subsample (the
# smoothScatter density is computed on full data; only the LOESS curve is
# on the subsample). At n ≈ 20k for Primary, full-data LOESS with span=0.5
# would run slowly; subsample is a faithful proxy of curvature.
#
# Input:  Analyses/6_1/derived/fitted_models.rds
#           list with $continuous = 12 fits, $logistic = 12 fits;
#           each fit entry: list(fit, n, formula, study_ids)
# Output: Analyses/6_1/derived/
#   diagnostics/6_1_continuous_qq.png            — 2×6 grid: outcome × variant
#   diagnostics/6_1_continuous_fitted_resid.png  — 2×6 grid + LOESS overlay
#   diagnostics/6_1_logistic_calibration.png     — 2×6 grid: outcome × variant
#   results/6_1_moments.csv                      — 12 rows, continuous residual moments
#   results/6_1_calibration_data.csv             — 120 rows, logistic deciles
#   results/6_1_cooks_summary.csv                — 16 rows, lab variants × frame × outcome
#   results/6_1_cooks_lab_full.csv               — per-observation Cook's (lab variants)
#   results/6_1_correlation_matrix.csv           — pairwise correlations among LE8 wearable predictors (Primary cohort, per SAP §6.1 Diagnostic / v0.70)
#   results/6_1_diagnostics.log                  — sink log

library(tidyverse)

project_dir  <- "./PRIMARY"
in_path      <- file.path(project_dir, "Analyses", "6_1", "derived", "fitted_models.rds")
diag_dir     <- file.path(project_dir, "Analyses", "6_1", "derived", "diagnostics")
results_dir  <- file.path(project_dir, "Analyses", "6_1", "derived", "results")
dir.create(diag_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "6_1_diagnostics.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("§6.1 diagnostics —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Load checkpoint ---
fitted_models   <- read_rds(in_path)
fits_continuous <- fitted_models$continuous
fits_logistic   <- fitted_models$logistic

cat("Loaded fits:\n")
cat("  continuous:", length(fits_continuous), "(expected 12)\n")
cat("  logistic:  ", length(fits_logistic),   "(expected 12)\n\n")

# --- Variant + outcome metadata ---
variant_order  <- c("Primary", "PlusDweight", "PlusDHDL", "PlusDNonHDL", "PlusDTG", "PlusDHbA1c")
variant_labels <- c(
  Primary     = "Primary",
  PlusDweight = "+Δweight",
  PlusDHDL    = "+ΔHDL",
  PlusDNonHDL = "+ΔNonHDL",
  PlusDTG     = "+ΔTG",
  PlusDHbA1c  = "+ΔHbA1c"
)
lab_variants  <- c("PlusDHDL", "PlusDNonHDL", "PlusDTG", "PlusDHbA1c")
outcome_order <- c("SBP", "DBP")

# Continuous fit keys: "<outcome>__<variant>" (e.g., "SBP__Primary")
# Logistic   fit keys: "<outcome>_responder__<variant>" (e.g., "SBP_responder__Primary")
parse_key_continuous <- function(key) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  list(outcome = parts[1], variant = parts[2])
}
parse_key_logistic <- function(key) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  list(responder = parts[1], variant = parts[2])
}

# --- Moments helper (reused pattern from §7.1) ---
moments_rows <- list()
moments_summary <- function(x, model_key, quantity) {
  x <- x[is.finite(x)]
  mu    <- mean(x)
  sigma <- sd(x)
  skew  <- mean((x - mu)^3) / sigma^3
  kurt  <- mean((x - mu)^4) / sigma^4 - 3  # excess kurtosis (normal = 0)
  cat(sprintf("  %-45s  n = %6d  mean = %+7.4f  sd = %6.3f  skew = %+6.3f  excess_kurt = %+6.3f\n",
              quantity, length(x), mu, sigma, skew, kurt))
  moments_rows[[length(moments_rows) + 1]] <<- tibble(
    model_key = model_key, quantity = quantity,
    n = length(x), mean = mu, sd = sigma,
    skewness = skew, excess_kurtosis = kurt
  )
  invisible(NULL)
}

# =============================================================================
# Section 1: Continuous diagnostics — residual moments + Q-Q + fitted/residual
# =============================================================================
cat("\n========== Section 1: Continuous ANCOVA diagnostics ==========\n")

continuous_diag <- lapply(names(fits_continuous), function(key) {
  fr   <- fits_continuous[[key]]
  meta <- parse_key_continuous(key)
  list(
    key       = key,
    label     = paste0(meta$outcome, " — ", variant_labels[[meta$variant]],
                       sprintf(" (n=%s)", format(fr$n, big.mark = ","))),
    outcome   = meta$outcome,
    variant   = meta$variant,
    residuals = residuals(fr$fit),
    fitted    = fitted(fr$fit),
    n         = fr$n
  )
})
names(continuous_diag) <- names(fits_continuous)

# Reorder for consistent panel grid: outcome first (rows), variant second (cols)
order_idx <- order(
  match(sapply(continuous_diag, `[[`, "outcome"), outcome_order),
  match(sapply(continuous_diag, `[[`, "variant"), variant_order)
)
continuous_diag <- continuous_diag[order_idx]

cat("\n=== Residual moments (continuous fits) ===\n")
cat("Flag values: |skew| > 1 = noticeable; |excess_kurt| > 1 = noticeable heavy/light tails.\n\n")
for (key in names(continuous_diag)) {
  d <- continuous_diag[[key]]
  cat("--- ", d$label, " ---\n", sep = "")
  moments_summary(d$residuals, key, "residuals")
}

moments_table <- bind_rows(moments_rows)
write_csv(moments_table, file.path(results_dir, "6_1_moments.csv"))
cat("\nWrote moments CSV:\n  ", file.path(results_dir, "6_1_moments.csv"), "\n")

# --- Q-Q residuals PNG: 2×6 grid (rows = outcome, cols = variant) ---
out_png_qq <- file.path(diag_dir, "6_1_continuous_qq.png")
png(out_png_qq, width = 2400, height = 1200, res = 150)
par(mfrow = c(2, 6), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

for (key in names(continuous_diag)) {
  d <- continuous_diag[[key]]
  qqnorm(d$residuals,
         main = d$label,
         pch = 16, cex = 0.3, col = adjustcolor("gray30", alpha = 0.4))
  qqline(d$residuals, col = "red", lty = 2, lwd = 1.5)
}
mtext("§6.1 continuous-ANCOVA Q-Q residuals — 6 variants × 2 outcomes",
      outer = TRUE, cex = 1.0, font = 2)
dev.off()
cat("\nWrote Q-Q grid:\n  ", out_png_qq, "\n")

# --- Fitted/residual PNG: 2×6 grid + LOESS overlay (subsampled to 5k) ---
out_png_fr <- file.path(diag_dir, "6_1_continuous_fitted_resid.png")
png(out_png_fr, width = 2400, height = 1200, res = 150)
par(mfrow = c(2, 6), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

set.seed(7)  # reproducible LOESS subsample
for (key in names(continuous_diag)) {
  d <- continuous_diag[[key]]
  smoothScatter(d$fitted, d$residuals,
                main = d$label,
                xlab = "Fitted value", ylab = "Residual",
                nrpoints = 0)
  abline(h = 0, col = "red", lty = 2, lwd = 1.5)

  # LOESS overlay — subsampled if n > 5000 to keep runtime manageable
  n_sample   <- min(5000, length(d$fitted))
  samp_idx   <- sample(seq_along(d$fitted), n_sample)
  loess_data <- data.frame(fitted    = d$fitted[samp_idx],
                           residuals = d$residuals[samp_idx])
  smooth_x   <- seq(min(loess_data$fitted), max(loess_data$fitted), length.out = 200)
  loess_fit  <- tryCatch(
    loess(residuals ~ fitted, data = loess_data, span = 0.5),
    error = function(e) NULL
  )
  if (!is.null(loess_fit)) {
    lines(smooth_x, predict(loess_fit, newdata = data.frame(fitted = smooth_x)),
          col = "blue", lwd = 2)
  }
}
mtext("§6.1 continuous-ANCOVA fitted vs. residual — 6 variants × 2 outcomes (LOESS subsampled to 5k)",
      outer = TRUE, cex = 1.0, font = 2)
dev.off()
cat("Wrote fitted/residual grid:\n  ", out_png_fr, "\n")

# =============================================================================
# Section 2: Logistic calibration data + plot
# =============================================================================
cat("\n========== Section 2: Logistic calibration ==========\n")

logistic_diag <- lapply(names(fits_logistic), function(key) {
  fr          <- fits_logistic[[key]]
  meta        <- parse_key_logistic(key)
  resp_short  <- sub("_responder$", "", meta$responder)  # "SBP" or "DBP"
  outcome_col <- paste0(resp_short, "_RESPONDER_6MO")

  d         <- fr$fit$model
  pred_prob <- predict(fr$fit, type = "response")
  obs       <- d[[outcome_col]]
  list(
    key       = key,
    label     = paste0(resp_short, "-responder — ", variant_labels[[meta$variant]],
                       sprintf(" (n=%s)", format(fr$n, big.mark = ","))),
    outcome   = resp_short,
    variant   = meta$variant,
    pred_prob = pred_prob,
    obs       = obs,
    n         = fr$n
  )
})
names(logistic_diag) <- names(fits_logistic)

# Reorder for consistent panel display (outcome × variant)
order_idx <- order(
  match(sapply(logistic_diag, `[[`, "outcome"), outcome_order),
  match(sapply(logistic_diag, `[[`, "variant"), variant_order)
)
logistic_diag <- logistic_diag[order_idx]

# Build calibration deciles per fit, long-format
calib_rows <- lapply(names(logistic_diag), function(key) {
  d <- logistic_diag[[key]]
  tibble(pred_prob = d$pred_prob, obs = d$obs) %>%
    mutate(pred_decile = ntile(pred_prob, 10)) %>%
    group_by(pred_decile) %>%
    summarise(n         = n(),
              mean_pred = mean(pred_prob),
              mean_obs  = mean(obs),
              .groups   = "drop") %>%
    mutate(model_key = key,
           outcome   = d$outcome,
           variant   = d$variant) %>%
    select(model_key, outcome, variant, pred_decile, n, mean_pred, mean_obs)
})
calib_data <- bind_rows(calib_rows)
write_csv(calib_data, file.path(results_dir, "6_1_calibration_data.csv"))

cat("\nCalibration deciles (long-format):", nrow(calib_data),
    "rows (expected 12 fits × 10 deciles = 120)\n")
cat("Wrote calibration CSV:\n  ", file.path(results_dir, "6_1_calibration_data.csv"), "\n")

# --- Calibration PNG: 2×6 grid, facet_grid(outcome ~ variant) ---
calib_plot_df <- calib_data %>%
  mutate(outcome     = factor(outcome, levels = outcome_order,
                              labels = c("SBP-responder", "DBP-responder")),
         variant_lbl = factor(variant_labels[variant], levels = variant_labels))

p_calib <- ggplot(calib_plot_df, aes(x = mean_pred, y = mean_obs)) +
  geom_abline(slope = 1, intercept = 0, color = "gray60", linetype = "dashed") +
  geom_point(size = 2, color = "#1f78b4") +
  geom_line(color = "#1f78b4", linewidth = 0.5) +
  facet_grid(outcome ~ variant_lbl, scales = "free") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Mean predicted probability (per decile)",
    y = "Observed proportion of responders",
    title    = "§6.1 logistic-responder calibration — 6 variants × 2 outcomes",
    subtitle = "10 deciles of predicted probability per fit; dashed = ideal calibration",
    caption  = "Per SAP §6.1 spec. Lab variants (HDL/NonHDL/TG/HbA1c) on small n (389–592)."
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title   = element_text(face = "bold"),
        strip.text   = element_text(face = "bold", size = 9),
        plot.caption = element_text(size = 9, color = "gray40", hjust = 0))

out_png_calib <- file.path(diag_dir, "6_1_logistic_calibration.png")
ggsave(out_png_calib, p_calib, width = 16, height = 7, dpi = 150)
cat("Wrote calibration grid:\n  ", out_png_calib, "\n")

# =============================================================================
# Section 3: Cook's distance — lab variants only (16 fits)
# =============================================================================
cat("\n========== Section 3: Cook's distance — lab variants only ==========\n")
cat("Lab variants:", paste(lab_variants, collapse = ", "), "\n")
cat("Skipping Primary + PlusDweight (n > 8k → 4/n threshold ≈ 0.0003, uninformative).\n\n")

cooks_summary_rows <- list()
cooks_full_rows    <- list()

process_cooks <- function(fits_list, frame_label, parse_fn) {
  for (key in names(fits_list)) {
    meta <- parse_fn(key)
    if (!(meta$variant %in% lab_variants)) next

    fr          <- fits_list[[key]]
    cooks_d_vec <- cooks.distance(fr$fit)
    threshold   <- 4 / fr$n
    above_idx   <- which(cooks_d_vec > threshold)

    outcome_short <- if (!is.null(meta$outcome)) meta$outcome else sub("_responder$", "", meta$responder)

    # Summary row
    cooks_summary_rows[[length(cooks_summary_rows) + 1]] <<- tibble(
      frame               = frame_label,
      variant             = meta$variant,
      outcome             = outcome_short,
      model_key           = key,
      n                   = fr$n,
      threshold_4_over_n  = threshold,
      max_cooks_d         = max(cooks_d_vec, na.rm = TRUE),
      n_above_threshold   = length(above_idx),
      pct_above_threshold = 100 * length(above_idx) / fr$n
    )

    # Per-observation rows — map Cook's values to STUDY_IDs via saved study_ids
    if (length(fr$study_ids) == length(cooks_d_vec)) {
      cooks_full_rows[[length(cooks_full_rows) + 1]] <<- tibble(
        frame              = frame_label,
        variant            = meta$variant,
        outcome            = outcome_short,
        model_key          = key,
        STUDY_ID            = fr$study_ids,
        cooks_d            = as.numeric(cooks_d_vec),
        threshold_4_over_n = threshold,
        above_threshold    = as.logical(cooks_d_vec > threshold)
      )
    } else {
      cat(sprintf("  WARNING: study_ids length (%d) ≠ cooks_d length (%d) for %s — skipping per-obs CSV rows.\n",
                  length(fr$study_ids), length(cooks_d_vec), key))
    }

    cat(sprintf("  %-50s n=%4d  threshold=%.4f  max_D=%.4f  above=%3d (%.1f%%)\n",
                key, fr$n, threshold, max(cooks_d_vec, na.rm = TRUE),
                length(above_idx), 100 * length(above_idx) / fr$n))
  }
}

cat("--- Continuous fits ---\n")
process_cooks(fits_continuous, "continuous", parse_key_continuous)

cat("\n--- Logistic fits ---\n")
process_cooks(fits_logistic, "logistic", parse_key_logistic)

cooks_summary <- bind_rows(cooks_summary_rows)
cooks_full    <- bind_rows(cooks_full_rows)

write_csv(cooks_summary, file.path(results_dir, "6_1_cooks_summary.csv"))
write_csv(cooks_full,    file.path(results_dir, "6_1_cooks_lab_full.csv"))

cat("\nWrote Cook's summary CSV (", nrow(cooks_summary), " rows):\n  ",
    file.path(results_dir, "6_1_cooks_summary.csv"), "\n", sep = "")
cat("Wrote Cook's per-observation CSV (", nrow(cooks_full), " rows, lab variants only):\n  ",
    file.path(results_dir, "6_1_cooks_lab_full.csv"), "\n", sep = "")

# =============================================================================
# Section 4: Pairwise correlation matrix among LE8 wearable predictors
# =============================================================================
# Per SAP §6.1 Diagnostic block — pairwise correlation matrix among the LE8
# wearable predictors as a collinearity-assessment diagnostic (complements the
# VIF check in `6_1_le8_regression.R` / `6_1_vif_check.csv`). Computed on the
# Primary-variant cohort.

cat("\n========== Section 4: Pairwise correlation matrix (LE8 wearable predictors) ==========\n")

primary_mf <- model.frame(fits_continuous[["SBP__Primary"]]$fit)
le8_predictors <- c("SLEEP_DURATION_STUDY_MEAN", "STRAIN_STUDY_MEAN",
                    "LOG_TOBACCO_NICOTINE_USE_DAYS", "TOBACCO_NICOTINE_EVER_LOGGED")

cor_mat <- primary_mf %>%
  dplyr::select(dplyr::any_of(le8_predictors)) %>%
  cor(use = "pairwise.complete.obs")

cor_long <- as.data.frame(as.table(cor_mat)) %>%
  dplyr::rename(predictor_1 = Var1, predictor_2 = Var2, correlation = Freq) %>%
  dplyr::mutate(dplyr::across(dplyr::starts_with("predictor"), as.character))

write_csv(cor_long, file.path(results_dir, "6_1_correlation_matrix.csv"))

cat("Pairwise correlations among LE8 wearable predictors (Primary cohort, n =",
    nrow(primary_mf), "):\n")
print(round(cor_mat, 3))
cat("\nWrote correlation matrix CSV:\n  ",
    file.path(results_dir, "6_1_correlation_matrix.csv"), "\n", sep = "")

# =============================================================================
# Section 5: Interpretation guide
# =============================================================================
cat("\n========== Diagnostic interpretation guide ==========\n")
cat("Q-Q residuals     — look for S-shape (heavy/light tails) or banana shape (skew).\n")
cat("                     Light deviation at extremes is OK; systematic deviation is a flag.\n")
cat("Fitted/residual   — blue LOESS line should be ~flat near y=0. Fan shape = heteroscedasticity.\n")
cat("Calibration       — points should fall near the dashed 1:1 line.\n")
cat("Cook's distance   — at lab-variant n ~400-600, 4/n threshold ≈ 0.007–0.010 is informative.\n")
cat("                     Observations above threshold flag potentially-influential users.\n")
cat("                     See 6_1_cooks_lab_full.csv for per-STUDY_ID values + filter on above_threshold = TRUE.\n")
cat("Correlation matrix — pairwise correlations among LE8 wearable predictors (Primary cohort).\n")
cat("                     |r| > 0.7 = strong collinearity concern; complement with VIF in 6_1_vif_check.csv.\n")

cat("\nAUC (logistic) and adjusted-R² (continuous) already in 6_1_model_summaries.csv\n")
cat("  (generated by 6_1_le8_regression.R).\n")

cat("\n=== §6.1 diagnostics done ===\n")
