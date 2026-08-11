# Blood Pressure Cohort Study — Analysis Code

Statistical analysis code for a retrospective observational cohort study of
six-month blood pressure change among WHOOP members with elevated or
hypertensive baseline blood pressure, and its association with wearable-derived
physiology, laboratory biomarkers, and platform engagement.

Published for transparency and methodological review. Manuscript in preparation.

---

## No data is included

**This repository contains no participant data of any kind** — no raw inputs, no
derived analytic frames, no intermediate outputs, no figures.

The study operates under an IRB condition requiring the analyst to work only with
de-identified data. The analytic tree replaces the platform user identifier with a
study surrogate (`STUDY_ID`); the identifier crosswalk is held separately and is
not part of this or any other repository. Participant-level data are not
shareable and are not distributed here.

The data-extraction layer (SQL against internal source systems) is also not
included. These scripts begin downstream of extraction, reading flat files that
the extraction layer produces.

`.gitignore` excludes `*.csv`, `*.rds`, `*.RData`, and all image formats, so data
and outputs cannot be committed by accident.

**Consequence: these scripts will not run as-published.** Every script expects
input files that are not in the repository. The code is here to be read and
audited — for the model specifications, covariate sets, sensitivity analyses, and
multiplicity handling — not executed.

---

## Layout

```
PRIMARY/
├── Build/       — construct the analytic frames from extracted flat files
└── Analyses/    — one directory per analysis section
```

Scripts resolve inputs through the relative path `./PRIMARY`. Opening
`bp.Rproj` in RStudio sets the working directory accordingly.

### Build

Numbered scripts run in order. Each produces one component frame; `99` assembles
them into the participant-level analytic frame and `98` runs quality-control
checks against it.

| Script | Builds |
|---|---|
| `01_bp_timepoints` | Blood pressure at each study timepoint |
| `02_bp_outcomes` | Change scores and responder classification |
| `03_demographics_clinical` | Demographics and clinical covariates |
| `04_lab_weight_outcomes` | Laboratory biomarkers and weight |
| `05_whoop_window_means` | Wearable metrics averaged over study windows |
| `06_le8_tobacco_nicotine` | Life's Essential 8 components |
| `07_engagement` | Platform engagement measures |
| `08_lmm_covariates` | Covariates for the mixed models |
| `98_qc_user_frame` / `99_assemble_user_frame` | QC and assembly |
| `bp_reading_long_frame` | Reading-level long frame for within-person models |

### Analyses

Directories follow the statistical analysis plan's section numbering.

| Section | Analysis |
|---|---|
| `0_cohort` | Cohort definitions, Table 1, completeness |
| `5_1` | Primary blood pressure trajectory (linear mixed model, by baseline category) |
| `5_2` | Proportions crossing clinical change thresholds |
| `5_3` | Cohort-wide distributional shift, baseline vs six months |
| `5_4` | Reclassification across AHA blood pressure categories |
| `5_5` | Twelve-month trajectory (paired subset) |
| `6_1` | Life's Essential 8–anchored multivariable regression |
| `6_2` | Logistic regression of responder status on wearable signals |
| `7_1` | Wearable-metric trajectories by responder status |
| `7_2` | Daily wearable signals predicting same-day blood pressure (within/between decomposition) |
| `8_engagement` | Engagement predictors of clinically meaningful response |
| `Abstract` | Abstract figures |

`Sensitivity/` subdirectories hold pre-specified sensitivity analyses for their
parent section — alternative cohort definitions, regression-to-the-mean checks,
covariate variations, and an alternative measurement cadence.

---

## Requirements

R, with:

```
tidyverse  data.table  lubridate
lme4  lmerTest  emmeans  broom  car  clubSandwich  performance  pROC
ggh4x  ggtext  legendry  patchwork
```

Mixed models are fit with `lme4`/`lmerTest`; estimated marginal means and
contrasts with `emmeans`.

---

## Status and use

Research code accompanying a manuscript in preparation, shared for review and
transparency, provided as-is with no warranty. The associated manuscript is the
authoritative reporting of this analysis.

**No license is granted.** This work was produced in the course of employment at
WHOOP; all rights are reserved by the copyright holder. The code is viewable here
but is not offered for reuse, redistribution, or derivative work. Please get in
touch before using any part of it.
