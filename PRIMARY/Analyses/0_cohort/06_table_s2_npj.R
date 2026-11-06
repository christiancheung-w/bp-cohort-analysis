# 06_table_s2_npj.R
# Table S2 for npj: the modeled-N table (04_modeled_n_word.R output) restricted
# to analyses retained in the npj manuscript.
#
# Dropped (analyses cut for npj or no longer reported):
#   LE8 lifestyle regressions and their sensitivities; WHOOP-metric trajectory
#   LMMs; Oldham/Blomqvist RTM decomposition (superseded by Barnett); empirical
#   median paired change; joint engagement logistic.
#
# Input:  Analyses/0_cohort/derived/results/table_A_A2_word_pasteable.tsv
# Output: Analyses/0_cohort/derived/results/table_s2_npj_paper.csv

suppressPackageStartupMessages(library(tidyverse))

results_dir <- "./PRIMARY/Analyses/0_cohort/derived/results"
lines <- read_lines(file.path(results_dir, "table_A_A2_word_pasteable.tsv"))

DROP <- c("^LE8 regression", "^WHOOP trajectory LMM", "Oldham/Blomqvist", "^Continuous-baseline-BP LMM",
          "^Empirical median", "^Engagement joint logistic")

section <- NA_character_; header <- NULL; rows <- list()
for (ln in lines) {
  f <- str_split(ln, "\t")[[1]]
  if (str_starts(f[1], "Table A.2")) { section <- "Sensitivity analyses"; next }
  if (str_starts(f[1], "Table A"))   { section <- "Primary analyses";     next }
  if (f[1] == "Analysis") { header <- f; next }
  if (all(f == "")) next
  rows[[length(rows) + 1]] <- c(section, f)
}
tab <- as_tibble(do.call(rbind, rows), .name_repair = "minimal")
names(tab) <- c("section", header)
stopifnot(!anyNA(tab$section), ncol(tab) == 10)

keep <- !Reduce(`|`, lapply(DROP, function(p) str_detect(tab$Analysis, p)))
out  <- tab[keep, ]
cat("Rows kept:", nrow(out), "of", nrow(tab), "| dropped analyses:\n")
print(unique(tab$Analysis[!keep]))
stopifnot(!any(str_detect(out$Cohort, "Normotensive")),
          !any(as.numeric(out$`% of ceiling`) > 100))

# Display labels for npj (numbers unchanged)
ANALYSIS <- c(
  "Primary BP-trajectory LMM"                      = "6-month BP trajectory (LMM)",
  "12-month BP-trajectory LMM (4-timepoint)"       = "12-month BP trajectory (LMM, 4 timepoints)",
  "WHOOP-signal logistic (dual-responder)"         = "Wearable metrics (logistic)",
  "W/B daily LMM (anchored BP − baseline_BP)"      = "Reading-level within/between-person model (LMM)",
  "W/B daily LMM (raw-BP)"                         = "Reading-level model, raw BP outcome (LMM)",
  "Engagement separate logistic"                   = "Engagement (logistic, separate models)",
  "RTM adjustment (Barnett et al. 2005)"           = "RTM correction, 6 months (Barnett)",
  "12-month RTM adjustment (Barnett et al. 2005)"  = "RTM correction, 12 months (Barnett)",
  "Continuous-baseline-BP LMM (cutpoint-free)"     = "6-month BP trajectory, continuous baseline BP (LMM)")
COHORT <- c("Hypertensive"                      = "Primary cohort",
            "Antihypertensive Medication Naive" = "No recorded antihypertensive medication",
            "2×Weekly"                          = "Twice-weekly readings",
            "Hypertensive (12-mo completers)"   = "12-month completers")
PREDSET <- c("Sleep need achieved"             = "Sleep need achieved",
             "Daily app time"                  = "Daily app time",
             "BP readings logged"              = "Blood pressure readings",
             "Days engaged with notifications" = "Notification engagement",
             "Days journaled"                  = "Journal completion",
             "Days met optimal strain"         = "Physical activity recommendation",
             "Days used AI Coach"              = "AI Coach engagement")
relabel <- function(x, map) ifelse(x %in% names(map), map[x], x)
out <- out %>%
  mutate(Analysis = case_when(str_starts(Analysis, "BP-trajectory LMM") ~ "6-month BP trajectory (LMM)",
                              str_starts(Analysis, fixed("W/B daily LMM (")) & !str_detect(Analysis, "raw-BP|anchored") ~
                                "Reading-level within/between-person model (LMM)",
                              TRUE ~ relabel(Analysis, ANALYSIS)),
         Cohort = relabel(Cohort, COHORT),
         `Predictor set` = relabel(`Predictor set`, PREDSET),
         Outcome = recode(Outcome, "BP responder" = "≥5 mmHg reduction"))
stopifnot(!any(str_detect(out$Analysis, "W/B|WHOOP-signal|Hypertensive")))

write_csv(out, file.path(results_dir, "table_s2_npj_paper.csv"))
cat("\nWrote:", file.path(results_dir, "table_s2_npj_paper.csv"), "\n")
