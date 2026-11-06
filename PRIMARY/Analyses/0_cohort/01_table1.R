# 01_table1.R
# -----------------------------------------------------------------------------
# Produces the paper Table 1 (hypertensive §2.2 cohort) plus the `full` and
# `htn_naive` sensitivity characteristic tables — each as a CSV and a
# paper-ready markdown — with between-stratum P-values. Sources the canonical
# definitions in cohort_defs.R.
#
# Retires cohort_characteristics.R (single producer for the cohort chapter).
#
# Run with working directory = project root (bp.Rproj folder).
#
# Outputs (Analyses/0_cohort/derived/results/):
#   table_1.csv / table_1_paper.md            — paper Table 1 (hypertensive)
#   table_1_full.csv / table_1_full_paper.md  — sensitivity: full anchored cohort
#   table_1_htn_naive.csv / ..._paper.md      — sensitivity: HTN-naive cohort
#   01_table1.log                             — sink log
# -----------------------------------------------------------------------------

library(tidyverse)
source("./PRIMARY/Analyses/0_cohort/cohort_defs.R")

project_dir <- "./PRIMARY"
input_path  <- file.path(project_dir, "Build", "derived", "qc", "user_analytic_frame_1X.rds")
results_dir <- file.path(project_dir, "Analyses", "0_cohort", "derived", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "01_table1.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("01_table1 —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Load ---
frame <- read_rds(input_path)
cat("Source:", input_path, "\n")
cat("Loaded:", nrow(frame), "rows x", ncol(frame), "cols\n\n")

# --- Build the three cohort tables -------------------------------------------
hyp   <- apply_cohort(frame, "hypertensive")
full  <- apply_cohort(frame, "full")
naive <- apply_cohort(frame, "htn_naive")

t_hyp   <- make_char_table(hyp,   COHORTS$hypertensive$strata)
t_full  <- make_char_table(full,  COHORTS$full$strata)
t_naive <- make_char_table(naive, COHORTS$htn_naive$strata)

# --- Write CSVs ---
write_csv(t_hyp,   file.path(results_dir, "table_1.csv"))
write_csv(t_full,  file.path(results_dir, "table_1_full.csv"))
write_csv(t_naive, file.path(results_dir, "table_1_htn_naive.csv"))

cat("--- Table 1 (paper, hypertensive §2.2) ---\n")
print(t_hyp, n = Inf, width = Inf)
cat("\n--- Sensitivity: full anchored cohort ---\n")
print(t_full, n = Inf, width = Inf)
cat("\n--- Sensitivity: HTN-naive cohort ---\n")
print(t_naive, n = Inf, width = Inf)

# =============================================================================
# Generic markdown renderer
# -----------------------------------------------------------------------------
# Builds the title + characteristic table for ANY cohort (3 or 4 strata). The
# logical row structure is identical across cohorts; only the number of stratum
# columns varies, so it is derived from the table itself. Footnotes (which carry
# cohort-specific counts) are assembled separately and appended by the caller.
# =============================================================================

get_val <- function(tbl, char, col) {
  tbl %>% filter(Characteristic == char) %>% pull(.data[[col]])
}
fmt_n <- function(x) format(x, big.mark = ",")

# Logical rows of the characteristic table, in display order. `kind`:
#   single = value row carrying its own p; header = group label, blank values,
#   p pulled from `p_from`; sub = indented value row, no p.
TABLE1_ROWS <- list(
  list(kind = "single", label = "Age, mean ± SD, y",              char = "Age, mean ± SD, y"),
  list(kind = "header", label = "Sex, n (%)",                     p_from = "Sex: Male, n (%)"),
  list(kind = "sub",    label = "    Male",                       char = "Sex: Male, n (%)"),
  list(kind = "sub",    label = "    Female",                     char = "Sex: Female, n (%)"),
  list(kind = "single", label = "BMI, mean ± SD, kg/m² ¶",        char = "BMI, mean ± SD, kg/m²"),
  list(kind = "single", label = "Baseline SBP, mean ± SD, mmHg",  char = "Baseline SBP, mean ± SD, mmHg"),
  list(kind = "single", label = "Baseline DBP, mean ± SD, mmHg",  char = "Baseline DBP, mean ± SD, mmHg"),
  list(kind = "header", label = "HTN medication status †, n (%)", p_from = "HTN meds: On at baseline, n (%)"),
  list(kind = "sub",    label = "    On at baseline",             char = "HTN meds: On at baseline, n (%)"),
  list(kind = "sub",    label = "    Initiated during study",     char = "HTN meds: Initiated during study, n (%)"),
  list(kind = "sub",    label = "    Never logged",               char = "HTN meds: Never logged, n (%)")
)

mk_row <- function(cells) paste0("| ", paste(cells, collapse = " | "), " |")

render_char_table_md <- function(tbl, title) {
  # Value columns = everything except the Characteristic label and the P-value
  # (preserves order: Overall, [Normal], Elevated, Stage 1, Stage 2).
  val_cols <- setdiff(names(tbl), c("Characteristic", "P-value"))

  # Header carries per-column N (from the N row).
  hdr_cells <- vapply(val_cols, function(s) {
    paste0(s, "<br>(n = ", get_val(tbl, "N", s), ")")
  }, character(1))
  header <- mk_row(c("Characteristic", hdr_cells, "P-value §"))
  sep    <- paste0("|:---|", strrep(":--:|", length(val_cols) + 1L))

  body <- vapply(TABLE1_ROWS, function(r) {
    if (r$kind == "header") {
      mk_row(c(r$label, rep("", length(val_cols)), get_val(tbl, r$p_from, "P-value")))
    } else {
      vals <- vapply(val_cols, function(s) get_val(tbl, r$char, s), character(1))
      p    <- if (r$kind == "single") get_val(tbl, r$char, "P-value") else ""
      mk_row(c(r$label, vals, p))
    }
  }, character(1))

  c(title, "",
    paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"), "",
    header, sep, body)
}

# --- Cohort-specific footnotes -----------------------------------------------
foot_tests <- function(n_strata, n_cohort) {
  word    <- c("3" = "three", "4" = "four")[as.character(n_strata)]
  n_round <- fmt_n(round(n_cohort / 1000) * 1000)
  paste0("<sub>§ P-values test differences across the ", word,
         " baseline AHA strata: one-way ANOVA (continuous — age, BMI, SBP, DBP) and Pearson chi-square (categorical — sex [Male vs Female], HTN medication status). **Baseline SBP and DBP differences are definitional** (the strata are defined by baseline BP) and are reported for completeness. At n ≈ ", n_round, ", p reflects statistical detectability, not effect magnitude.</sub>")
}

foot_htn <- "<sub>† HTN medication status per §12.8 derivation. **On at baseline** = `stable_on_meds` ∪ `recent_initiator` (qualifying HTN-med journal entry on or before `STUDY_START_DATE`, including the 180-day lookback). **Initiated during study** = `mid_study_initiator` (first qualifying entry on/after `STUDY_START_DATE`). **Never logged** = no qualifying HTN-med journal entry in lookback or study window. Self-reported via app journal; under-ascertainment of users who do not log meds is a known limitation (§12.8).</sub>"

foot_sex <- function(d, with_modeling) {
  n_total <- nrow(d)
  n_mf    <- sum(d$SEX %in% c("male", "female"))
  n_nb    <- sum(d$SEX == "non-binary")
  n_pn    <- sum(d$SEX == "i prefer not to choose")
  s <- paste0("<sub>‡ Sex row sums to ", fmt_n(n_mf), " (Male + Female), which is less than Overall N = ",
              fmt_n(n_total), ". Excluded from the Sex row: ", n_nb, " `non-binary` + ", n_pn,
              " `i prefer not to choose` users. These users are included in the N total and all other rows.")
  if (with_modeling) {
    s <- paste0(s, " Modeling analyses (§5.1 / §5.5 / §6.1 / §6.2 / §7.1 / §7.2 / §8) additionally exclude these users from the analytic cohort (modeled n = ",
                fmt_n(n_mf), "); see the supplementary modeled-n table for per-analysis n's.")
  }
  paste0(s, "</sub>")
}

foot_bmi <- function(d) {
  n_miss <- sum(is.na(d$BASELINE_BMI))
  paste0("<sub>¶ BMI missing for ", n_miss, " (",
         sprintf("%.2f", 100 * n_miss / nrow(d)),
         "%) users; the mean ± SD is computed over non-missing users. All other Table 1 variables are complete.</sub>")
}

foot_abbrev <- "<sub>AHA = American Heart Association; BMI = body mass index; SBP / DBP = systolic / diastolic blood pressure.</sub>"

# --- Assemble + write the three markdowns ------------------------------------
write_md <- function(lines, file) {
  writeLines(lines, file.path(results_dir, file))
  cat("  wrote", file, "\n")
}

cat("\nWriting paper-ready markdown:\n")

write_md(c(
  render_char_table_md(t_hyp, "# Table 1. Cohort characteristics by baseline AHA stratum (§2.2 primary cohort)"),
  "", foot_tests(length(COHORTS$hypertensive$strata), nrow(hyp)),
  "", foot_htn,
  "", foot_sex(hyp, with_modeling = TRUE),
  "", foot_bmi(hyp),
  "", foot_abbrev
), "table_1_paper.md")

write_md(c(
  render_char_table_md(t_full, "# Cohort characteristics by baseline AHA stratum — full anchored cohort (§2.1, incl. Normal) [sensitivity]"),
  "", foot_tests(length(COHORTS$full$strata), nrow(full)),
  "", foot_htn,
  "", foot_sex(full, with_modeling = FALSE),
  "", foot_bmi(full),
  "", foot_abbrev
), "table_1_full_paper.md")

write_md(c(
  render_char_table_md(t_naive, "# Cohort characteristics by baseline AHA stratum — HTN-naive cohort (never-logged) [sensitivity]"),
  "", foot_tests(length(COHORTS$htn_naive$strata), nrow(naive)),
  "", foot_htn,
  "", foot_sex(naive, with_modeling = FALSE),
  "", foot_bmi(naive),
  "", foot_abbrev
), "table_1_htn_naive_paper.md")

cat("\nAll outputs written to:\n  ", results_dir, "\n")
cat("\n=== 01_table1 done ===\n")
