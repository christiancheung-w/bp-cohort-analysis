# 5_1_diagnostics.R
# SAP §5.1 diagnostic plots: Q-Q residuals, fitted/residual, Q-Q BLUPs.
# Per SAP §5.1 diagnostic spec — Flag criteria:
#   Q-Q residuals: clear S-shape (heavy/light tails), banana shape (skew),
#                  isolated extreme points
#   Fitted/residual: fan-shaped spread, systematic curvature
#   Q-Q BLUPs: same criteria as residual Q-Q
#
# At n ≈ 50k observations, formal normality tests (Shapiro-Wilk, KS) reject
# any deviation from normality. Per SAP v0.4 the diagnostic intent is
# descriptive flag-spotting, not formal hypothesis test. We use visual Q-Q +
# numeric moment summaries (mean, SD, skewness, excess kurtosis).
#
# Input:  Analyses/5_1/derived/fitted_models.rds (4 lmerMod objects)
# Output: Analyses/5_1/derived/diagnostics/5_1_diagnostics.png (4 × 3 grid)
#         Console: numeric residual + BLUP moment summaries

library(tidyverse)
library(lme4)

project_dir  <- "./PRIMARY"
in_path      <- file.path(project_dir, "Analyses", "5_1", "derived", "fitted_models.rds")
diag_dir     <- file.path(project_dir, "Analyses", "5_1", "derived", "diagnostics")
results_dir  <- file.path(project_dir, "Analyses", "5_1", "derived", "results")
dir.create(diag_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Capture full console output to a log file (split = TRUE → also show on console)
log_path <- file.path(results_dir, "5_1_diagnostics.log")
sink(log_path, split = TRUE)
cat("§5.1 diagnostics run —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

fitted_models <- read_rds(in_path)

# --- Moment summary helper ---
moments_rows <- list()  # accumulate for CSV export
moments_summary <- function(x, name, quantity) {
  x <- x[is.finite(x)]
  mu    <- mean(x)
  sigma <- sd(x)
  skew  <- mean((x - mu)^3) / sigma^3
  kurt  <- mean((x - mu)^4) / sigma^4 - 3  # excess kurtosis (normal = 0)
  cat(sprintf("  %-30s  n = %6d  mean = %+7.4f  sd = %6.3f  skew = %+6.3f  excess_kurt = %+6.3f\n",
              quantity, length(x), mu, sigma, skew, kurt))
  moments_rows[[length(moments_rows) + 1]] <<- tibble(
    model        = name,
    quantity     = quantity,
    n            = length(x),
    mean         = mu,
    sd           = sigma,
    skewness     = skew,
    excess_kurtosis = kurt
  )
  invisible(NULL)
}

# --- Per-model: collect residuals, fitted values, BLUPs ---
diag_data <- function(fit, label) {
  list(
    label     = label,
    residuals = residuals(fit),
    fitted    = fitted(fit),
    blups     = ranef(fit)$STUDY_ID[, "(Intercept)"]
  )
}

dlist <- list(
  sbp_unadj = diag_data(fitted_models$sbp_unadj, "SBP — unadjusted"),
  sbp_adj   = diag_data(fitted_models$sbp_adj,   "SBP — adjusted"),
  dbp_unadj = diag_data(fitted_models$dbp_unadj, "DBP — unadjusted"),
  dbp_adj   = diag_data(fitted_models$dbp_adj,   "DBP — adjusted")
)

# --- Console: numeric moment summaries ---
cat("\n=== Residual + BLUP moment summaries ===\n")
cat("(For reference: a normal distribution has skew = 0, excess_kurt = 0.\n")
cat(" Flag values: |skew| > 1 = noticeable skew; |excess_kurt| > 1 = noticeable heavy/light tails.)\n\n")

for (name in names(dlist)) {
  d <- dlist[[name]]
  cat("--- ", d$label, " ---\n", sep = "")
  moments_summary(d$residuals, name, "residuals")
  moments_summary(d$blups,     name, "BLUPs")
  cat("\n")
}

# Save structured moments as CSV
moments_table <- bind_rows(moments_rows)
write_csv(moments_table, file.path(results_dir, "5_1_moments.csv"))
cat("\nWrote moments CSV:\n  ", file.path(results_dir, "5_1_moments.csv"), "\n")

# --- PNG: 4 × 3 grid of diagnostic plots ---
out_png <- file.path(diag_dir, "5_1_diagnostics.png")
png(out_png, width = 1800, height = 2400, res = 150)
par(mfrow = c(4, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

for (name in names(dlist)) {
  d <- dlist[[name]]

  # Q-Q residuals
  qqnorm(d$residuals,
         main = paste0(d$label, "\nQ-Q residuals"),
         pch = 16, cex = 0.25, col = adjustcolor("gray30", alpha = 0.4))
  qqline(d$residuals, col = "red", lty = 2, lwd = 1.5)

  # Fitted vs residual — smoothScatter handles 50k points
  smoothScatter(d$fitted, d$residuals,
                main = paste0(d$label, "\nFitted vs residual"),
                xlab = "Fitted value", ylab = "Residual",
                nrpoints = 0)
  abline(h = 0, col = "red", lty = 2, lwd = 1.5)
  # LOESS overlay for systematic-curvature inspection
  smooth_x <- seq(min(d$fitted), max(d$fitted), length.out = 200)
  loess_fit <- tryCatch(
    loess(residuals ~ fitted, data = data.frame(fitted = d$fitted, residuals = d$residuals),
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

mtext("SAP §5.1 LMM diagnostics — 4 models × {Q-Q residuals, fitted/residual + LOESS, Q-Q BLUPs}",
      outer = TRUE, cex = 0.9, font = 2)
dev.off()

cat("\nWrote diagnostic grid PNG:\n  ", out_png, "\n")
cat("\nInterpretation guide:\n")
cat("  Q-Q residuals     — look for S-shape (heavy/light tails) or banana shape (skew).\n")
cat("                       Light deviation at extremes is OK; systematic deviation is a flag.\n")
cat("  Fitted/residual   — blue LOESS line should be ~flat near y=0. Fan shape = heteroscedasticity.\n")
cat("  Q-Q BLUPs         — checks normality of user-level random intercept distribution.\n")

# Close the sink to flush log file
sink()
