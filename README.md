# Blood Pressure Cohort Study — Analysis Code

R code for a retrospective observational cohort study of blood pressure (BP)
change among users of a consumer wearable device and its companion digital
health platform, and of the association of that change with platform
engagement and wearable-derived physiology.

This repository is the code companion to the manuscript that cites it. It
contains the complete analysis code behind every table, figure and reported
result in that manuscript and its supplement. Analyses that were explored
during the project but are not reported in the manuscript are not included.

---

## No data is included

**This repository contains no participant data of any kind** — no raw inputs,
no derived analytic frames, no intermediate outputs, no figures.

The study operates under an IRB condition requiring the analyst to work only
with de-identified data. The analytic tree replaces the platform user
identifier with a study surrogate (`STUDY_ID`); the identifier crosswalk is
held separately and is not part of this or any other repository.
Participant-level data are not shareable and are not distributed here.

The data-extraction layer (queries against internal source systems) is also not
included. These scripts begin downstream of extraction, reading flat files that
the extraction layer produces. Some cohort-defining steps are therefore applied
upstream of this code and only validated here: the age range, anchoring of the
study period on each user's first cuff-measured entry, exclusion of automated
accounts, the antihypertensive-medication tracker definition and its look-back
window, and the denominators of the pre-aggregated engagement metrics (days with
wearable data, days with an activity recommendation, days of app use). The
scripts re-derive the BP timepoints and the medication groups from the extracted
readings and check them against the extraction layer's values.

`.gitignore` excludes `*.csv`, `*.rds`, `*.RData`, and all image formats, so
data and outputs cannot be committed by accident.

**Consequence: these scripts will not run as-published.** Every script expects
input files that are not in the repository. The code is here to be read and
audited — for the cohort definitions, model specifications, covariate sets,
sensitivity analyses and multiplicity handling — not executed.

---

## Layout

```
run_all.R        — runs every script below in dependency order
PRIMARY/
├── Build/       — construct the analytic frames from extracted flat files
└── Analyses/    — one directory per analysis section
```

Scripts resolve inputs through the relative path `./PRIMARY`. Opening
`bp.Rproj` in RStudio sets the working directory accordingly. `run_all.R`
lists the 58 scripts in the order they were run; it is the authoritative
dependency order.

### Build

Numbered scripts run in order. Each produces one component frame; `99`
assembles them into the participant-level analytic frame and `98` applies the
quality-control bounds (including the BP plausibility range) to it.

| Script | Builds |
|---|---|
| `01_bp_timepoints` | Baseline, 3-, 6- and 12-month BP readings from the reading history |
| `02_bp_outcomes` | Change scores and ≥5 mmHg responder classification |
| `03_demographics_clinical` | Demographics, BMI, antihypertensive-medication classification |
| `04_lab_weight_outcomes` | Laboratory and weight variables (frame columns only; see note below) |
| `05_whoop_window_means` | Wearable metrics averaged over the pre-baseline and study windows |
| `06_le8_tobacco_nicotine` | Tobacco / nicotine use and related journal variables |
| `07_engagement` | Platform engagement measures |
| `08_lmm_covariates` | Season and time-of-day covariates for the mixed models |
| `09_glp1_flag` | GLP-1 receptor agonist flag for the sensitivity cohort |
| `98_qc_user_frame` / `99_assemble_user_frame` | QC bounds and assembly of the participant frame |
| `bp_reading_long_frame` | Reading-level long frame for the within/between-person models |
| `Sensitivity/2xWeekly/build_user_frame` | Frame for the twice-weekly-readings sensitivity cohort |

The Build layer is included whole so that the analytic frame is fully
specified. Two Build scripts (`04`, `06`) create columns for variables that
are summarised descriptively or not reported in the manuscript; no model in
this repository uses the laboratory or weight columns.

### Analyses

Directories keep the section numbering of the study's statistical analysis
plan. The table below maps each manuscript item to the script(s) that produce
it. `cohort_defs.R` holds the cohort definitions and is sourced by most
analysis scripts.

| Manuscript item | Script(s) |
|---|---|
| Table 1 | `0_cohort/07_table1_s1_npj_20003.R` |
| Table 2 (6-month change, RTM-corrected columns, fit rows) | `5_1/5_1_primary_lmm.R` → `5_1/5_1_table_2.R` → `5_1/5_1_table_2_npj.R`; fit rows from `5_1/5_1_diagnostics.R`; RTM from `5_1/Sensitivity/rtm_barnett/` and `5_1/Sensitivity/rtm_barnett_bootstrap/` |
| Figure 1A (trajectories) | `5_1/5_1_figure_1a.R` |
| Figure 1B (baseline vs 6-month distributions) | `5_3/5_3_distributional_comparison.R` |
| Figure 1C (improved / stable / worsened) | `5_2/5_2_clinical_bins.R` |
| Figure 1D (AHA category transitions) | `5_4/5_4_aha_category_shift.R` |
| Figure 2 (engagement odds ratios) | `8_engagement/8_engagement_logistic.R` → `8_engagement/8_figure_forest.R` |
| Figure 3 (wearable-metric odds ratios) | `6_2/6_2_logistic.R` (model and figure) |
| Figure 4 (reading-level between/within coefficients) | `7_2/7_2_wb_lmm.R` → `7_2/7_2_figure_4.R` |
| Twice-weekly-readings sensitivity (Results text) | `Build/Sensitivity/2xWeekly/`, `5_1/Sensitivity/2xWeekly/`, `5_1/5_1_comparison_1X_vs_2xWeekly.R` |
| Table S1 (completeness) | `0_cohort/02_completeness_supp.R`, recounted in `0_cohort/07_table1_s1_npj_20003.R` |
| Table S2 (participants and observations per model) | `0_cohort/03_modeled_n_supp.R` → `04_modeled_n_word.R` → `06_table_s2_npj.R` → `08_table_s2_npj_add_cohorts.R` |
| Table S3 (6-month change: sensitivity cohorts and subgroups) | `5_1/Sensitivity/{full_cohort,htn_naive,age_sex}/`, `glp1_naive_sensitivity/`, `mid_excluded_sensitivity/` → `5_1/5_1_table_sensitivity.R` |
| Table S4 (12-month model) and Figure S1 | `5_5/5_5_primary_lmm.R` → `5_5/5_5_table_s3.R` → `5_5/5_5_table_s3_npj.R`; `5_5/5_5_diagnostics.R`; `5_5/Sensitivity/rtm_barnett/`; `5_5/5_5_figure_2.R` |
| Figure S2 (12-month AHA transitions) | `5_5/5_5_aha_distribution_12mo.R` |
| Table S5 (engagement, sensitivity cohorts) | `8_engagement/8_table_sensitivity.R` |
| Table S6 (wearable metrics, sensitivity cohorts) | `6_2/6_2_table_sensitivity.R` |
| Table S7 (reading-level model, sensitivity analyses) | `7_2/7_2_sensitivities.R` → `7_2/7_2_table_sensitivity.R` |

Figure 1 is assembled from the four panel PNGs outside this repository. The
panel-letter suffixes in the `5_2` and `5_3` output file names predate the
final panel order; the mapping above is the one used in the manuscript.

Notes on scope:

- `glp1_naive_sensitivity.R` and `mid_excluded_sensitivity.R` refit every
  manuscript model in the GLP-1-naive and medication-initiator-excluded
  cohorts. The GLP-1 script also contains sections for two analyses that were
  removed from the manuscript; the source scripts for those analyses are not
  included and those sections are not used by any reported table.
- `03_modeled_n_supp.R` collects model sizes from every analysis directory it
  finds; `06_table_s2_npj.R` then restricts Table S2 to the analyses reported.
- `0_cohort/01_table1.R`, `02_completeness_supp.R` and `05_table1_npj_rows.R`
  are earlier descriptive scripts computed before the final sex restriction
  (20,021 participants); `07_table1_s1_npj_20003.R` recomputes Table 1 and
  Table S1 on the final analytic cohort and is the version the manuscript
  reports.

---

## Requirements

R (the analysis was run under R 4.6.0), with:

```
tidyverse  data.table  lubridate
lme4  lmerTest  emmeans  broom  car  clubSandwich  performance  pROC
ggtext  legendry  scales  knitr
```

Mixed models are fit with `lme4`/`lmerTest`; estimated marginal means and
contrasts with `emmeans`; cluster-robust standard errors with `clubSandwich`;
R² with `performance`.

---

## Provenance

The scripts are exported unchanged from the frozen analysis repository used to
produce the manuscript's results, with two exceptions: an internal project
prefix was removed from the input file names, and two header comments were
generalised. No statistical code was altered.

---

## License

Released under the BSD 3-Clause License (see `LICENSE`). Copyright WHOOP,
Inc. The code is provided as-is, without warranty. The associated manuscript
is the authoritative report of this analysis.
