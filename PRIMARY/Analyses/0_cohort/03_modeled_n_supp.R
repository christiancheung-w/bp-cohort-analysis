# 03_modeled_n_supp.R
# -----------------------------------------------------------------------------
# Supplementary "Table A": per-analysis modeled n. Aggregates the standardized
# <chapter>_sizes_std.csv files that each modeled chapter emits via
# emit_sizes() (cohort_defs.R). Numbers are read from the chapters' fitted
# objects (approach A2) — never recomputed here — so Table A cannot disagree
# with the reported model results.
#
# Split by `role` into two tables:
#   Table A    — primary modeled analyses (the headline n's)
#   Table A.2  — sensitivity-analysis n's (companion)
#
# ROLLING deliverable: it contains whichever chapters have been instrumented so
# far. It grows as each chapter / sensitivity is reviewed.
#
# Run with working directory = project root.
#
# Outputs (Analyses/0_cohort/derived/results/):
#   table_A_modeled_n.csv          (primary)
#   table_A2_sensitivity_n.csv     (sensitivity)
#   table_A_modeled_n_paper.md     (both, with descriptions)
#   03_modeled_n_supp.log
# -----------------------------------------------------------------------------

library(tidyverse)

project_dir  <- "./PRIMARY"
analyses_dir <- file.path(project_dir, "Analyses")
results_dir  <- file.path(analyses_dir, "0_cohort", "derived", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(results_dir, "03_modeled_n_supp.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)
cat("03_modeled_n_supp —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Collect every chapter's standardized size file --------------------------
files <- list.files(analyses_dir, pattern = "_sizes_std\\.csv$",
                    recursive = TRUE, full.names = TRUE)
cat("Found", length(files), "chapter size file(s):\n")
for (f in files) cat("  ", f, "\n")
if (length(files) == 0) stop("No *_sizes_std.csv files found — run the instrumented chapter scripts first.")

all_sizes <- map_dfr(files, ~ read_csv(.x, show_col_types = FALSE)) %>%
  mutate(`% of ceiling` = sprintf("%.1f", 100 * n_users / ceiling_n))

primary <- all_sizes %>% filter(role == "primary")     %>% arrange(section, cohort, outcome, predictor_set)
sens    <- all_sizes %>% filter(role == "sensitivity")  %>% arrange(section, cohort, outcome, predictor_set)

write_csv(primary, file.path(results_dir, "table_A_modeled_n.csv"))
write_csv(sens,    file.path(results_dir, "table_A2_sensitivity_n.csv"))
cat("\n--- Table A (primary) ---\n");     print(as.data.frame(primary), row.names = FALSE)
cat("\n--- Table A.2 (sensitivity) ---\n"); print(as.data.frame(sens),  row.names = FALSE)

# --- Markdown rendering -------------------------------------------------------
mk_row <- function(cells) paste0("| ", paste(cells, collapse = " | "), " |")
HDR <- c("Section", "Analysis", "Cohort", "Outcome", "Predictor set",
         "n_users", "n_obs", "n_events", "ceiling_n", "% of ceiling")
render_tbl <- function(df) {
  if (nrow(df) == 0) return("_(none instrumented yet)_")
  fmtnum <- function(x) ifelse(is.na(x), "", format(x, big.mark = ","))
  rows <- map_chr(seq_len(nrow(df)), function(i) {
    r <- df[i, ]
    mk_row(c(r$section, r$analysis, r$cohort, r$outcome, r$predictor_set,
             fmtnum(r$n_users), fmtnum(r$n_obs), fmtnum(r$n_events),
             fmtnum(r$ceiling_n), r$`% of ceiling`))
  })
  c(mk_row(HDR), paste0("|", strrep("---|", length(HDR))), rows)
}

md <- c(
  "# Table A. Per-analysis modeled n (transparency & reproducibility supplement)",
  "",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "_"),
  "",
  "Each row is one fitted model. The n's are **read directly from the fitted model objects** (not recomputed for this table), so they match the reported results exactly. This is a rolling supplement — it contains the analyses instrumented so far, split into the primary analyses (Table A) and the sensitivity analyses (Table A.2).",
  "",
  "**Column definitions**",
  "",
  "- **n_users** — distinct participants contributing to the model.",
  "- **n_obs** — observations (model rows): equals `n_users` for cross-sectional models; for longitudinal LMMs it is participant × timepoint, so `n_obs` > `n_users` when participants contribute multiple timepoints.",
  "- **n_events** — number of outcome events for logistic models (e.g. responders); blank for continuous outcomes.",
  "- **ceiling_n** — the maximum possible N for the analysis: its cohort entry count *before* any outcome/covariate complete-casing. For the sex-restricted modeling cohorts this is the modeled N (e.g. 20,003) — i.e. Table B's descriptive denominator (20,021) minus the sex-excluded users.",
  "- **% of ceiling** — `n_users / ceiling_n`: the share of the cohort the model retained. Attrition below 100% is explained variable-by-variable in **Table B** (e.g. a model requiring ΔHbA1c retains only the ~1.9% of the cohort with paired HbA1c; an LMM needing only baseline covariates retains ~100% bar the few BMI-missing).",
  "",
  "## Table A. Primary modeled analyses",
  "",
  render_tbl(primary),
  "",
  "## Table A.2. Sensitivity analyses",
  "",
  render_tbl(sens)
)
writeLines(md, file.path(results_dir, "table_A_modeled_n_paper.md"))
cat("\nWrote: table_A_modeled_n.csv, table_A2_sensitivity_n.csv, table_A_modeled_n_paper.md\n")
cat("\n=== 03_modeled_n_supp done ===\n")
