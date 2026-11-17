# Analysis code: physiological modifiers of plasma Alzheimer disease biomarkers

R code for the data preparation and every analysis, table and figure in "Physiological Variation
in Plasma Alzheimer Disease Biomarkers: Concentration and Assay Effects" (Yang et al.). Pipeline version 1.2.0 (`PAPER1_VERSION` in `00_setup.R`).

## Data

No data are included. Both cohorts are available to qualified investigators on application:

- ADNI, including the FNIH Biomarkers Consortium plasma measures
  (`FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.csv`): https://adni.loni.usc.edu
- A4 and LEARN: https://www.a4studydata.org

The authors are not permitted to redistribute either dataset. The code reads the released files
unchanged and never writes to the raw folders. The ADNI files read are listed in `00_setup.R`
(`ADNI_FILES`): the FNIH plasma file, ADNIMERGE, the Berkeley amyloid PET file, LABDATA, VITALS,
RECCMEDS and RECMHIST. The A4 and LEARN files are listed in `A4_FILES`: the Roche plasma
panel, the p-tau217 file and the clinical laboratory file (External Data), SUBJINFO (Derived Data)
and vitals (Raw Data).

## Requirements

R 4.5 with `dplyr`, `tidyr`, `readr`, `stringr`, `lubridate`, `purrr`, `tibble`, `survival`,
`lme4`, `nlme`, `lavaan` (the MIMIC model only; skipped with a message if absent), `splines` and
`parallel`; `rstudioapi` is optional. The analyses for the paper used R 4.5.1; the session
information of each run is written to `logs/`.

## Paths

`00_setup.R` takes three folders from environment variables:

| Variable | Meaning | Default |
|---|---|---|
| `PAPER1_ADNI_DIR` | folder with the ADNI files (read only); several folders may be given, separated by semicolons | `ADNI/` beside the folder that holds the code |
| `PAPER1_A4_DIR` | folder with the A4 and LEARN release (read only) | `A4&learn/` beside the folder that holds the code |
| `PAPER1_OUT_DIR` | folder for everything the pipeline writes | `paper1_out/` beside the folder that holds the code |

## How to run

The files live in one flat folder and source `00_setup.R` by file name, so keep them together.
`run_all.R` runs everything in order:

```r
Sys.setenv(PAPER1_ADNI_DIR = "<ADNI folder>", PAPER1_A4_DIR = "<A4 and LEARN folder>",
           PAPER1_OUT_DIR = "<output folder>")
STEPS <- "all"     # or, for example, c("09", "10")
QUICK <- FALSE     # TRUE skips the joint growth model bootstraps and stage 08
source("<code folder>/run_all.R")
```

`STEPS` selects stages by their two-digit prefix. Each stage reads what the earlier stages wrote,
so a stage can be rerun on its own once the earlier ones have run. With `QUICK <- FALSE`, the
default, every interval is computed as in the paper: 1,000 resamples of participants, seed 1. That
run takes about a day, almost all of it in the bootstraps of the joint growth model (stages 06, 07
and 08). `QUICK <- TRUE` turns off the bootstraps of the joint growth model and skips stage 08; the
1,000-resample bootstraps of sections 4.2, 5.4 and 7.3, which take minutes, still run, and the
point estimates are unchanged. Write each run to an empty output folder: a rerun into a used one
reuses some earlier results (the joint model bootstrap rows of section 6.1, the comorbidity
interval of section 7.2 when its bootstrap is off, the saved tasks of section 8.2, and the
laboratory extract of section 2.1). The stage 08 tasks run side by side on all but one
core (`PAPER1_MI_CORES` sets the number); on a Mac with Apple's Accelerate BLAS, set
`VECLIB_MAXIMUM_THREADS=1` in `~/.Renviron` so that the parallel fits do not compete for threads.

## Files, in run order

| File | What it does |
|---|---|
| `00_setup.R` | Paths, file names, constants, the definition of the eleven factors and their signs, and the estimation functions used throughout: least squares with standard errors clustered by participant, left-censored regression for C2N p-tau217, the equivalence test, the Deming slope, the reliability model, regression calibration, and the censored eGFR spline. |
| `01_data_adni.R` | The ADNI analysis table: plasma on every platform, one row per participant and draw, with amyloid PET, the laboratory panel, eGFR, body mass index, statin and kidney-acting medication. |
| `02_data_a4learn.R` | The A4 and LEARN tables: the screening cross-section, the repeated p-tau217 draws and every laboratory visit. |
| `03_levels.R` | Research question 1, levels in both cohorts, with the eGFR spline, mutual adjustment, regression calibration, the negative control, sex strata, laboratory recency, and the level models under six specifications. |
| `04_platform_differences.R` | Research question 2, the difference between platforms on the same aliquot for every pair and factor, equivalence, the subsamples, slope correction, the MIMIC model, the C2N detection limit, the creatine kinase interactions, the C2N differences as a difference of coefficients, and the date of the draw. |
| `05_ratios.R` | Research question 3, the %p-tau217 ratio and its components; the Fujirebio reporting floor; the two post hoc Fujirebio ratios; clinical units; the cleared ratio corrected for body mass index against amyloid PET. |
| `06_change.R` | Research question 4, in A4: the hybrid model, the reliability of change, the joint growth model with its bootstrap and span rules, and the joint model on a clinical scale and by trial arm. |
| `07_sensitivity_and_robustness.R` | Classification agreement between p-tau217 platforms, a visit effect shared by both series, creatinine in place of eGFR, E-values, hypertension and diabetes as covariates, and agreement with amyloid PET by physiological factor. |
| `08_intervals.R` | Intervals added in the review, for estimates of stages 05 to 07: the cleared ratio against PET, the joint model with trial arm terms, the corrected within-participant estimates, and two eGFR rate correlations. |
| `09_tables.R` | Tables 1 to 5, the supplementary tables S1 to S18 and S26, the additions to Table 1, S20 with the direction measure, and plain-language labels in every `Table*.csv` and `S*.csv` file (S19 to S25 and S27 are written in stage 07, and S28 to S32 in stage 05). No model is fitted. |
| `10_figures.R` | Figures 1 and 2 and eFigure 2. |
| `run_all.R` | The driver. |
| `PIPELINE.md` | Every stage step by step: inputs, models, outputs, and where each result appears in the paper. |

Each stage file is divided into sections. Each section runs in its own environment and carries
the name and version of the script it was developed as, so the provenance of every table is
visible; a comment such as "as in step 06" inside a section refers to that earlier script, which
the table below maps to its section. Where a section uses only part of an earlier script, the
lines it needs from the top of that script are repeated in it.

| Earlier script | Now |
|---|---|
| `01_build_adni.R` | 1.1 |
| `02_build_a4.R` | 2.1 |
| `03_rq1_levels.R` | 3.1 |
| `04_rq2_platform.R` | 4.1 |
| `05_rq3_ratio.R` | 5.1 |
| `06_rq4_change.R` | 6.1 |
| `07_tables.R` | 9.1 |
| `09_supplementary_tables.R` | 9.2 |
| `10_robustness.R` | 7.1 |
| `11_figures.R` | 10.1 |
| `12_comorbidity_and_labels.R` | 7.2 (sections 1 to 2c), 5.2 (2d), 9.4 (3 and 4) |
| `13_fujirebio_ratios_and_floor.R` | 5.3 |
| `14_review_checks.R` | 3.2 (section 7), 4.2 (1, 1b, 2), 5.4 (5), 6.2 (3), 7.3 (6, 6b), 9.3 (4) |
| `15_review_intervals.R` | 8.1 |
| `16_missing_intervals.R` | 8.2 |

## Outputs

All under the output folder: the analysis tables as `.rds`; `tables/` (every estimate, as written
by the analysis stages, with full precision); `tables/review/` (the analyses added in the review);
`manuscript/` (the formatted tables); `checks/` (counts, matched medication and history strings,
and the raw files read, with their sizes and dates); `figures/`; and `logs/` (the log of each run
and its session information). `PIPELINE.md` names the file behind every table and figure.
