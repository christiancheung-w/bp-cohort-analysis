# 7_1_diagnostics.R
# SAP §7.1 (v0.64+) diagnostic plots: Q-Q residuals, fitted/residual, Q-Q BLUPs.
# Now spans the DUAL-PRIMARY structure: 5 WHOOP metrics × 2 responder defs (SBP, DBP) = 10 fits.
#
# Per SAP §7.1 diagnostic spec — Flag criteria:
#   Q-Q residuals: clear S-shape (heavy/light tails), banana shape (skew),
#                  isolated extreme points
#   Fitted/residual: fan-shaped spread, systematic curvature
#   Q-Q BLUPs: same criteria as residual Q-Q
#
# At the §2.2 cohort scale, formal normality tests (Shapiro-Wilk, KS) reject
# any deviation from normality. Per the §5.1 v0.4 framing the diagnostic
# intent is descriptive flag-spotting, not formal hypothesis test. We use
# visual Q-Q + numeric moment summaries (mean, SD, skewness, excess kurtosis).
#
# Nakagawa R² (marginal + conditional) is already computed and saved by
# 7_1_whoop_trajectories.R at derived/results/7_1_r2.csv — not duplicated here.
#
# Input:  Analyses/7_1/derived/fitted_models.rds
#         (list keyed by responder_def, each value a list of 5 metric fits)
# Outputs:
#   Analyses/7_1/derived/diagnostics/7_1_diagnostics_sbp.png (5 × 3 grid, SBP frame)
#   Analyses/7_1/derived/diagnostics/7_1_diagnostics_dbp.png (5 × 3 grid, DBP frame)
#   Analyses/7_1/derived/results/7_1_moments.csv (with responder_def column)
#   Analyses/7_1/derived/results/7_1_diagnostics.log

library(tidyverse)
library(lme4)

project_dir  <- "./PRIMARY"
in_path      <- file.path(project_dir, "Analyses", "7_1", "derived", "fitted_models.rds")
diag_dir     <- file.path(project_dir, "Analyses", "7_1", "derived", "diagnostics")
results_dir  <- file.path(project_dir, "Analyses", "7_1", "derived", "results")
dir.create(diag_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "7_1_diagnostics.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("§7.1 diagnostics run (dual-primary, 10 fits) —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

fits_by_def <- read_rds(in_path)
if (!all(c("SBP", "DBP") %in% names(fits_by_def))) {
  stop("Expected `fits_by_def` keyed by 'SBP' and 'DBP'; got: ",
       paste(names(fits_by_def), collapse = ", "))
}

# --- Metric labels (parallels `metrics` + `metric_units` in 7_1_whoop_trajectories.R) ---
metric_labels <- c(
  RHR              = "RHR (bpm)",
  SleepDuration    = "Sleep duration (hours)",
  SleepConsistency = "Sleep consistency (% 0-100)",
  Strain           = "Strain (0-21)",
  HRV              = "HRV (ms)"
)

# --- Moment summary helper (carries responder_def) ---
moments_rows <- list()
moments_summary <- function(x, responder_def, model_name, quantity) {
  x <- x[is.finite(x)]
  mu    <- mean(x)
  sigma <- sd(x)
  skew  <- mean((x - mu)^3) / sigma^3
  kurt  <- mean((x - mu)^4) / sigma^4 - 3  # excess kurtosis (normal = 0)
  cat(sprintf("  %-30s  n = %6d  mean = %+7.4f  sd = %6.3f  skew = %+6.3f  excess_kurt = %+6.3f\n",
              quantity, length(x), mu, sigma, skew, kurt))
  moments_rows[[length(moments_rows) + 1]] <<- tibble(
    responder_def   = responder_def,
    model           = model_name,
    quantity        = quantity,
    n               = length(x),
    mean            = mu,
    sd              = sigma,
    skewness        = skew,
    excess_kurtosis = kurt
  )
  invisible(NULL)
}

# --- Per-fit: collect residuals, fitted values, BLUPs ---
diag_data <- function(fit, label) {
  list(
    label     = label,
    residuals = residuals(fit),
    fitted    = fitted(fit),
    blups     = ranef(fit)$STUDY_ID[, "(Intercept)"]
  )
}

# --- Storage-pattern helper: unwrap list(fit, data, metric) or bare merMod ---
get_lmer <- function(x) {
  if (inherits(x, "merMod")) {
    x
  } else if (is.list(x) && "fit" %in% names(x) && inherits(x$fit, "merMod")) {
    x$fit
  } else {
    stop("Unrecognized fitted-models entry: expected merMod or list with $fit.")
  }
}

# --- Build per-(responder_def × metric) diagnostic-data list ---
dlist_by_def <- list()
for (rd in c("SBP", "DBP")) {
  fits <- fits_by_def[[rd]]
  dlist_by_def[[rd]] <- lapply(names(fits), function(m) {
    diag_data(get_lmer(fits[[m]]), metric_labels[[m]])
  })
  names(dlist_by_def[[rd]]) <- names(fits)
}

# --- Console: numeric moment summaries (per responder_def × metric) ---
cat("\n=== Residual + BLUP moment summaries ===\n")
cat("(For reference: a normal distribution has skew = 0, excess_kurt = 0.\n")
cat(" Flag values: |skew| > 1 = noticeable skew; |excess_kurt| > 1 = noticeable heavy/light tails.)\n\n")

for (rd in c("SBP", "DBP")) {
  cat(sprintf("\n##### Responder def: %s\n\n", rd))
  dlist <- dlist_by_def[[rd]]
  for (name in names(dlist)) {
    d <- dlist[[name]]
    cat("--- ", d$label, " ---\n", sep = "")
    moments_summary(d$residuals, rd, name, "residuals")
    moments_summary(d$blups,     rd, name, "BLUPs")
    cat("\n")
  }
}

moments_table <- bind_rows(moments_rows)
write_csv(moments_table, file.path(results_dir, "7_1_moments.csv"))
cat("\nWrote moments CSV:\n  ", file.path(results_dir, "7_1_moments.csv"), "\n")

# --- PNG factory: one 5 × 3 grid per responder_def ---
write_diag_png <- function(rd, dlist, out_png) {
  png(out_png, width = 1800, height = 3000, res = 150)
  par(mfrow = c(5, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  for (name in names(dlist)) {
    d <- dlist[[name]]

    # Q-Q residuals
    qqnorm(d$residuals,
           main = paste0(d$label, "\nQ-Q residuals"),
           pch = 16, cex = 0.25, col = adjustcolor("gray30", alpha = 0.4))
    qqline(d$residuals, col = "red", lty = 2, lwd = 1.5)

    # Fitted vs residual — smoothScatter handles large n
    smoothScatter(d$fitted, d$residuals,
                  main = paste0(d$label, "\nFitted vs residual"),
                  xlab = "Fitted value", ylab = "Residual",
                  nrpoints = 0)
    abline(h = 0, col = "red", lty = 2, lwd = 1.5)
    smooth_x <- seq(min(d$fitted), max(d$fitted), length.out = 200)
    loess_fit <- tryCatch(
      loess(residuals ~ fitted,
            data = data.frame(fitted = d$fitted, residuals = d$residuals),
            span = 0.5, control = loess.control(surface = "direct")),
      error = function(e) NULL
    )
    if (!is.null(loess_fit)) {
      lines(smooth_x, predict(loess_fit, newdata = data.frame(fitted = smooth_x)),
            col = "blue", lwd = 2)
    }

    # Q-Q BLUPs
    qqnorm(d$blups,
           main = paste0(d$label, "\nQ-Q user random intercepts (BLUPs)"),
           pch = 16, cex = 0.4, col = adjustcolor("gray30", alpha = 0.5))
    qqline(d$blups, col = "red", lty = 2, lwd = 1.5)
  }

  mtext(sprintf("SAP §7.1 LMM diagnostics — %s-responder frame (5 metrics × {Q-Q residuals, fitted/residual + LOESS, Q-Q BLUPs})", rd),
        outer = TRUE, cex = 0.9, font = 2)
  dev.off()
  cat("\nWrote diagnostic grid PNG:\n  ", out_png, "\n")
}

write_diag_png("SBP", dlist_by_def[["SBP"]],
               file.path(diag_dir, "7_1_diagnostics_sbp.png"))
write_diag_png("DBP", dlist_by_def[["DBP"]],
               file.path(diag_dir, "7_1_diagnostics_dbp.png"))

cat("\nInterpretation guide:\n")
cat("  Q-Q residuals     — look for S-shape (heavy/light tails) or banana shape (skew).\n")
cat("                       Light deviation at extremes is OK; systematic deviation is a flag.\n")
cat("  Fitted/residual   — blue LOESS line should be ~flat near y=0. Fan shape = heteroscedasticity.\n")
cat("  Q-Q BLUPs         — checks normality of user-level random intercept distribution.\n")
cat("\nNakagawa R² (marginal + conditional) already in derived/results/7_1_r2.csv\n")
cat("  (generated by 7_1_whoop_trajectories.R; now has responder_def column).\n")
