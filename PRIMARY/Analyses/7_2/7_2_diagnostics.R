# 7_2_diagnostics.R
# SAP §7.2 LMM diagnostics: residual Q-Q, fitted-vs-residual + LOESS, BLUP Q-Q.
# Per audit finding 7.2.2 (diagnostics missing from primary script).
#
# Inputs:  Analyses/7_2/derived/fitted_models.rds
#          Analyses/7_2/derived/cohort_wb_checkpoint.rds
# Outputs: Analyses/7_2/derived/figures/diagnostics/7_2_diagnostics.png
#            (2 × 3 grid: outcomes × {Q-Q resid, fitted/resid + LOESS, BLUP Q-Q})
#          Analyses/7_2/derived/results/7_2_diagnostics.log  (Caveat: when sourced
#                                                             via source(..., echo = TRUE)
#                                                             the on.exit + sink interaction
#                                                             captures only the initial
#                                                             echo line; analytic outputs
#                                                             are produced correctly
#                                                             regardless. Same caveat
#                                                             as §5.1 / §6.2 / §7.1.)
#          Analyses/7_2/derived/results/7_2_moments.csv

library(tidyverse)
library(lme4)

project_dir <- "./PRIMARY"
out_dir     <- file.path(project_dir, "Analyses", "7_2", "derived")
results_dir <- file.path(out_dir, "results")
diag_dir    <- file.path(out_dir, "figures", "diagnostics")
dir.create(diag_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "7_2_diagnostics.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)
cat("§7.2 LMM diagnostics —",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Load fits + cohort_wb checkpoint ---
fits      <- read_rds(file.path(out_dir, "fitted_models.rds"))
cohort_wb <- read_rds(file.path(out_dir, "cohort_wb_checkpoint.rds"))
cat("Loaded:", nrow(cohort_wb), "rows from",
    n_distinct(cohort_wb$STUDY_ID), "users\n\n")

# --- Build diagnostic list per outcome ---
outcomes <- list(
  SBP = list(fit = fits$SBP, label = "SBP − baseline_SBP"),
  DBP = list(fit = fits$DBP, label = "DBP − baseline_DBP")
)

dlist <- lapply(outcomes, function(o) {
  list(
    label     = o$label,
    residuals = residuals(o$fit),
    fitted    = fitted(o$fit),
    blups     = ranef(o$fit)$STUDY_ID[, "(Intercept)"]
  )
})

# --- Moments table ---
moments_rows <- list()
for (name in names(dlist)) {
  d <- dlist[[name]]
  std_resid <- (d$residuals - mean(d$residuals)) / sd(d$residuals)
  moments_rows[[name]] <- tibble(
    outcome      = d$label,
    n_obs        = length(d$residuals),
    n_users      = length(d$blups),
    resid_mean   = mean(d$residuals),
    resid_sd     = sd(d$residuals),
    resid_skew   = mean(std_resid^3),
    resid_xskurt = mean(std_resid^4) - 3,
    blup_mean    = mean(d$blups),
    blup_sd      = sd(d$blups)
  )
}
moments_table <- bind_rows(moments_rows)
write_csv(moments_table, file.path(results_dir, "7_2_moments.csv"))
cat("\nMoments table:\n")
print(moments_table, width = Inf)

# --- PNG: 2 × 3 grid (outcomes × {Q-Q resid, fitted/resid + LOESS, BLUP Q-Q}) ---
out_png <- file.path(diag_dir, "7_2_diagnostics.png")
png(out_png, width = 1800, height = 1200, res = 150)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

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
                xlab = "Fitted value (mmHg)", ylab = "Residual (mmHg)",
                nrpoints = 0)
  abline(h = 0, col = "red", lty = 2, lwd = 1.5)
  # LOESS overlay for systematic-curvature inspection (direct surface for large n)
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

mtext("SAP §7.2 LMM diagnostics — 2 outcomes × {Q-Q residuals, fitted/residual + LOESS, Q-Q BLUPs}",
      outer = TRUE, cex = 0.9, font = 2)
dev.off()

cat("\nWrote diagnostic grid PNG:\n  ", out_png, "\n")
cat("\nInterpretation guide:\n")
cat("  Q-Q residuals     — look for S-shape (heavy/light tails) or banana shape (skew).\n")
cat("                       Light deviation at extremes is OK; systematic deviation is a flag.\n")
cat("  Fitted/residual   — blue LOESS line should be ~flat near y=0. Fan shape = heteroscedasticity.\n")
cat("  Q-Q BLUPs         — checks normality of user-level random intercept distribution.\n")
cat("\nNakagawa R² (marginal + conditional) is in 7_2_model_summaries.csv\n")
cat("  (generated by 7_2_wb_lmm.R).\n")

cat("\n=== §7.2 diagnostics done ===\n")
