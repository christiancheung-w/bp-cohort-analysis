# 5_5_diagnostics.R
# SAP §5.5 diagnostic plots: Q-Q residuals, fitted/residual, Q-Q BLUPs.
# Identical structure to 5_1_diagnostics.R, applied to the 4-timepoint LMM
# fits in Analyses/5_5/derived/fitted_models.rds.
#
# Input:  Analyses/5_5/derived/fitted_models.rds (4 lmerMod objects)
# Output: Analyses/5_5/derived/diagnostics/5_5_diagnostics.png  (4 × 3 grid)
#         Analyses/5_5/derived/results/5_5_moments.csv
#         Analyses/5_5/derived/results/5_5_diagnostics.log

library(tidyverse)
library(lme4)

project_dir  <- "./PRIMARY"
in_path      <- file.path(project_dir, "Analyses", "5_5", "derived", "fitted_models.rds")
diag_dir     <- file.path(project_dir, "Analyses", "5_5", "derived", "diagnostics")
results_dir  <- file.path(project_dir, "Analyses", "5_5", "derived", "results")
dir.create(diag_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "5_5_diagnostics.log")
sink(log_path, split = TRUE)
cat("§5.5 diagnostics run —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

fitted_models <- read_rds(in_path)

# --- Moment summary helper ---
moments_rows <- list()
moments_summary <- function(x, name, quantity) {
  x <- x[is.finite(x)]
  mu    <- mean(x)
  sigma <- sd(x)
  skew  <- mean((x - mu)^3) / sigma^3
  kurt  <- mean((x - mu)^4) / sigma^4 - 3
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

diag_data <- function(fit, label) {
  list(
    label     = label,
    residuals = residuals(fit),
    fitted    = fitted(fit),
    blups     = ranef(fit)$STUDY_ID[, "(Intercept)"]
  )
}

dlist <- list(
  sbp_unadj = diag_data(fitted_models$sbp_unadj, "SBP — unadjusted (4-tp)"),
  sbp_adj   = diag_data(fitted_models$sbp_adj,   "SBP — adjusted (4-tp)"),
  dbp_unadj = diag_data(fitted_models$dbp_unadj, "DBP — unadjusted (4-tp)"),
  dbp_adj   = diag_data(fitted_models$dbp_adj,   "DBP — adjusted (4-tp)")
)

cat("\n=== Residual + BLUP moment summaries ===\n")
cat("(Flag thresholds: |skew| > 1 or |excess_kurt| > 1 → noticeable departure. Per SAP §5.1 v0.4 + v0.69 project-wide diagnostic-threshold convention.)\n\n")

for (name in names(dlist)) {
  d <- dlist[[name]]
  cat("--- ", d$label, " ---\n", sep = "")
  moments_summary(d$residuals, name, "residuals")
  moments_summary(d$blups,     name, "BLUPs")
  cat("\n")
}

moments_table <- bind_rows(moments_rows)
write_csv(moments_table, file.path(results_dir, "5_5_moments.csv"))
cat("\nWrote moments CSV:\n  ", file.path(results_dir, "5_5_moments.csv"), "\n")

# --- PNG: 4 × 3 grid ---
out_png <- file.path(diag_dir, "5_5_diagnostics.png")
png(out_png, width = 1800, height = 2400, res = 150)
par(mfrow = c(4, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

for (name in names(dlist)) {
  d <- dlist[[name]]

  qqnorm(d$residuals,
         main = paste0(d$label, "\nQ-Q residuals"),
         pch = 16, cex = 0.25, col = adjustcolor("gray30", alpha = 0.4))
  qqline(d$residuals, col = "red", lty = 2, lwd = 1.5)

  smoothScatter(d$fitted, d$residuals,
                main = paste0(d$label, "\nFitted vs residual"),
                xlab = "Fitted value", ylab = "Residual",
                nrpoints = 0)
  abline(h = 0, col = "red", lty = 2, lwd = 1.5)
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

  qqnorm(d$blups,
         main = paste0(d$label, "\nQ-Q user random intercepts (BLUPs)"),
         pch = 16, cex = 0.4, col = adjustcolor("gray30", alpha = 0.5))
  qqline(d$blups, col = "red", lty = 2, lwd = 1.5)
}

mtext("SAP §5.5 LMM diagnostics — 4 models × {Q-Q residuals, fitted/residual + LOESS, Q-Q BLUPs}",
      outer = TRUE, cex = 0.9, font = 2)
dev.off()

cat("\nWrote diagnostic grid PNG:\n  ", out_png, "\n")
cat("\nInterpretation guide (same as §5.1):\n")
cat("  Q-Q residuals     — S-shape = heavy/light tails; banana = skew. Light extreme deviation OK.\n")
cat("  Fitted/residual   — flat blue LOESS near y=0; fan shape = heteroscedasticity.\n")
cat("  Q-Q BLUPs         — normality of user random intercept distribution.\n")
cat("  At n ≈ 50k, descriptive flag-spotting only; LMM inference robust to mild departures.\n")

sink()
