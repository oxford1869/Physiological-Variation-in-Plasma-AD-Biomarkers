# The analysis pipeline, stage by stage

A specification of what the code in this folder does, for someone checking the code against the
paper. Each stage gives its inputs, the models as the code implements them, the files it writes,
and where each result appears in the manuscript or the supplement. Section numbers (3.1, 4.2, ...)
are the numbered sections inside each stage file; each section begins with a banner naming the
script it was developed as. File names such as `T2_platform_difference.csv` are in `tables/` of
the output folder unless another folder is named.

Conventions used throughout (`00_setup.R`):

- Biomarkers are analysed on the natural log scale, so a coefficient is reported as a percentage
  difference, 100 × (exp(b) − 1), per standard deviation of the factor.
- Each factor is standardised within its cohort and signed so that a positive estimate means the
  direction in its label: eGFR, haemoglobin and albumin per SD lower; body mass index, creatine
  kinase (log), urea, AST (log), ALT (log), HbA1c, C-reactive protein (log) and glucose per SD
  higher. The four primary factors are eGFR, body mass index, creatine kinase and haemoglobin.
- Covariates of the level and difference models: age at the draw, sex, APOE4, Centiloid,
  diagnosis and time since the first draw in ADNI (`COVS_ADNI`); age, sex, APOE4 and Centiloid in
  A4 and LEARN (`COVS_A4`).
- Least squares with CR1 standard errors clustered by participant (`fit_cl`). C2N p-tau217, which
  reports 0.65 pg/mL for values below its limit of 1.30 pg/mL, is fitted by left-censored
  regression censored at 1.30 pg/mL with cluster-robust standard errors (`fit_tobit`); the
  %p-tau217 ratio is censored at 100 × 1.30 / np-tau217, and a between-platform difference with
  C2N at log(1.30) minus the other platform's log value.
- Equivalence between two platforms: the 90% interval of the difference lies within ±5% per SD
  (`equivalence`). Verdict: equivalent first; otherwise "differ" when the Benjamini-Hochberg q is
  below 0.05; otherwise no difference detected. False discovery rate families as in eMethods 9.
- Bootstraps resample participants with replacement, 1,000 times, seed 1; for the joint growth
  model each drawn participant is given a new identifier and only converged fits are used.

---

## Stage 1. ADNI data (`01_data_adni.R`, section 1.1)

Reads the FNIH plasma file, ADNIMERGE, the Berkeley amyloid PET file, LABDATA, VITALS and RECCMEDS.
Plasma dates within 7 days are one draw. The sample is every participant with at least one draw
measured by a platform other than Quanterix; draws measured only by Quanterix are dropped.
Amyloid PET is matched to each draw within 180 days and weight within 365 days; height is that
visit's measured height, or the participant's median measured height where none was measured. eGFR (CKD-EPI 2021) is from the creatinine nearest the first
draw; the laboratory panel is from screening, with an unscheduled re-draw used only to fill a test
with no valid screening value. Statin and kidney-acting medication are set at each draw (a record
dated on or before the draw, or undated); every matched string is written out.

Writes `adni_visit.rds`, `adni_long.rds`, `checks/adni_build.csv`, `checks/adni_sources.csv`,
`checks/adni_statin_strings.csv`, `checks/adni_renal_medication_strings.csv`.
Reported in eTables 1, 2 and 13 (statin and kidney-acting classes).

## Stage 2. A4 and LEARN data (`02_data_a4learn.R`, section 2.1)

Reads SUBJINFO, the clinical laboratory file (extract cached as `a4_lab_raw_v3.rds`), the Roche
screening panel and the p-tau217 file. Laboratory tests are selected by code. C-reactive protein
below its limit is set to the limit. The screening cross-section keeps the most recent record of a
re-screened participant; laboratory values are matched within 90 days of the screening draw and
45 days of a later draw. p-tau217 readings outside the quantification range are used, with flags.

Writes `a4_screen.rds`, `a4_long.rds`, `a4_lab.rds`, `checks/a4_build.csv`, `checks/a4_sources.csv`.
Reported in eTables 1 and 2.

## Stage 3. Levels (`03_levels.R`)

**3.1** Every factor against every biomarker on every platform, both cohorts (`T1_levels.csv`;
q within cohort and factor role). Regression calibration with the reliability of the factor given
the covariates, from the A4 repeats (`R1_reliability.csv`, columns `rel_conditional`,
`pct_corrected` in T1). eGFR as a natural cubic spline (knots 45, 60, 75), with a cluster-robust
Wald test against a straight line and the difference from eGFR 90 at 60, 45 and 30; C2N p-tau217
by the censored spline (`T1b_egfr_shape.csv`). The four primary factors in one model
(`T1c_mutual_adjustment.csv`); Centiloid as the outcome (`T1d_negative_control.csv`); creatine
kinase by sex (`T1e_sex_stratified.csv`); ADNI laboratory values within 1 and 3 years of each draw
(`T1f_lab_recency.csv`).

**3.2** The four primary factors in ADNI under six specifications: the main model, without
diagnosis, Centiloid as a natural spline (3 df) with and without diagnosis, the first draw only,
and cognitively unimpaired amyloid-negative participants (`review/review_7_level_sensitivity.csv`).
The main rows must reproduce T1.

Reported in Table 2 (T1), eTables 3 (T1), 4 (T1c), 5 (T1, R1), 6 (T1f), 7 (T1e), 8 (T1b),
9 (T1d), 28 (T1) and 35 (review_7), and the abstract and Results.

## Stage 4. Differences between platforms (`04_platform_differences.R`)

**4.1** For every pair of platforms measuring the same analyte and every factor, the log
difference regressed on the factor and the covariates, with equivalence and verdicts
(`T2_platform_difference.csv`, `T2b_equivalence_summary.csv`). Slope correction by Deming slopes
at error-variance ratios of 0.50, 1.00 and 2.00 (columns `slope_*`, `corrected_*` of T2; pairs with
C2N p-tau217 on the samples above its limit). Creatine kinase in subsamples: statin users and
non-users at each draw, statin as a covariate, kidney-acting medication, men, women, the first
draw with both platforms, laboratory values within 1 and 3 years (`T2c_sensitivity.csv`). The MIMIC
model, omnibus test of any assay-specific effect by the scaled difference test, p-tau217 without
C2N and, as a sensitivity analysis, with C2N values below the limit set to missing
(`T3_latent_direct_effects.csv`). C2N p-tau217 as supplied, with the samples below the limit
excluded, and censored (`T6_floor_sensitivity.csv`); the C2N pairs with the samples below the limit
set aside (`T6b_difference_floor_excluded.csv`).

**4.2** Creatine kinase and the two differences, C2N minus ALZpath p-tau217 and Roche minus
Quanterix p-tau181, by sex, statin use and kidney-acting medication, in fully interacted models
(`review/review_1_ck_interactions.csv`); the C2N differences as the C2N censored level coefficient
minus the other platform's least squares coefficient, with bootstrap intervals
(`review/review_1b_c2n_difference_of_coefficients.csv`); the eGFR and creatine kinase differences
with the calendar date of the draw added (`review/review_2_draw_date.csv`). Each model first
reproduces T2.

Reported in Table 3 (T2), Figure 1 (T2), eTables 10 (T2b), 11 (T2; supplied as a data file),
12 (T2c), 14 (T3), 15 (T2), 16 (T6, with T4b), 17 (T6b), 36 (review_1), 37 (review_1b) and
38 (review_2).

## Stage 5. Ratios (`05_ratios.R`)

**5.1** The %p-tau217 ratio, its numerator (censored) and its denominator, for every factor, with
the verdict that the ratio cancels (equivalent), does not cancel (q below 0.05), or is
inconclusive (`T4_ratio_components.csv`); the three handlings of the C2N limit, each with its
interval and q (`T4b_ratio_floor.csv`).

**5.2** The ratio's body mass index association per SD, per 2 SD and per 5 kg/m² (`T12_ratio_bmi.csv`).

**5.3** The Fujirebio p-tau217 reporting floor (47 values at 0.030 pg/mL) left-censored in the
level models, the three Fujirebio pair differences, the MIMIC model (floor values missing) and the
Deming slopes (`T13a` to `T13d`); the post hoc ratios, the FDA-cleared Fujirebio p-tau217 over
Fujirebio amyloid-beta 42 and Fujirebio p-tau217 over Roche amyloid-beta 42, decomposed as the
%p-tau217 ratio is (`T13e_posthoc_fujirebio_ratios.csv`), with their sensitivity analyses
(`T13f_posthoc_ratio_sensitivity.csv`); the cleared categories against amyloid PET as counts
(`checks/cleared_ratio_calibration.csv`); estimates per doubling of creatine kinase, per
15 mL/min/1.73 m² lower eGFR and per 5 kg/m² higher body mass index (`T13h_clinical_units.csv`).

**5.4** Each sample's cleared ratio corrected for its own body mass index, r × exp(−b × z), and its
categories against amyloid PET at Centiloid 24 and 20, in PET-positive and PET-negative samples,
with bootstrap intervals (`review/review_5_bmi_correction_reclassification.csv`).

Reported in Table 4 (T4), eTables 16 (T4b), 26 (T4b), 30 (T13a, T13b), 31 (T13e), 32 (T13f),
33 (review_5) and 34 (T13h), eMethods 5 and 8.9 (T13c, T13d, the calibration counts), and the
Results (T12).

## Stage 6. Change on change (`06_change.R`)

**6.1** In A4 trial participants: the hybrid model, each factor split into the participant's mean
and the deviation from it (`T5_hybrid_within_between.csv`); the reliability of a within-person
change from the repeated laboratory values (`R2_change_reliability.csv`); the corrected
within-person estimate where that reliability is 0.20 or more. The joint growth model: both series
standardised, time centred at 2.5 years, a random level and rate for each series in each
participant with an unstructured covariance, separate residual variances, and the level and rate
of both series adjusted for age at the first draw, sex, APOE4 and Centiloid; the estimand is the
correlation between the two rates of change, fitted overall and in each trial arm, with a second
optimiser as a check (`T5b_joint_growth.csv`). Eligibility: at least two laboratory values
spanning a year and at least two p-tau217 draws. Bootstrap intervals for eGFR, creatine kinase and
haemoglobin, overall and by arm (`T5c_joint_bootstrap.csv`). The eGFR correlation under no span
rule, a one-year span on both series, out-of-range p-tau217 draws excluded, and a quadratic term
for time, each with a bootstrap interval (`T5d_span_sensitivity.csv`).

**6.2** From the same eGFR model: the change in the p-tau217 rate per 1 mL/min/1.73 m² per year
faster eGFR decline and per SD, the SD of the true eGFR rate, and the difference in correlation
between the solanezumab and placebo arms, with bootstrap intervals
(`review/review_3_joint_model.csv`); the model with trial arm on the level and rate of both series
(`review/review_3b_arm_by_time.csv`). Each fit first reproduces T5b.

Reported in Table 5 (T5, R2, T5b, T5c), eTables 18 (T5, R2), 19 (T5c, T5b, review_3, review_3b)
and 27 (T5d).

## Stage 7. Sensitivity analyses and robustness checks (`07_sensitivity_and_robustness.R`)

**7.1** Each p-tau217 platform given an amyloid threshold by Youden's index against PET at
Centiloid 20 and 24 (`T7b_thresholds.csv`); the proportion of samples two platforms classify
differently, and its change per SD of each factor, with the direction of disagreement
(`T7_classification.csv`). The eGFR joint growth model with the residuals of the two series
correlated within a visit (laboratory values and draws whose follow-up times round to the same
0.10 years), and the within-visit
residual correlation (`T10_occasion_check.csv`). Log creatinine in place of eGFR
(`T9_creatinine.csv`). E-values for the primary associations (`T8_evalues.csv`). Formatted as
`manuscript/S19.csv` to `S25.csv`.

**7.2** Hypertension and diabetes: in ADNI from the medical history, medication and blood pressure
at the nearest vitals visit; in A4 and LEARN from screening blood pressure and HbA1c (definitions
in the section). The level models, the two creatine kinase differences and the eGFR joint growth
model, each with and without the two conditions on the same rows, with a bootstrap interval for
the adjusted correlation (`T11_comorbidity.csv`, `manuscript/S27.csv`,
`checks/comorbidity_build.csv`, and the matched strings in `checks/`).

**7.3** Sensitivity and specificity of each p-tau217 platform at its threshold, and of the cleared
ratio, per SD of each primary factor at a fixed amyloid burden (Centiloid spline, age, sex)
(`review/review_6_accuracy_per_sd.csv`), and by band, descriptively
(`review/review_6b_accuracy_by_band.csv`).

Reported in eTables 13 (hypertension and diabetes strings), 20 (T7b), 21 and 22 (T7), 23 (T8),
24 (T9), 25 (T10), 29 (T11), 39 (review_6) and 40 (review_6b).

## Stage 8. Intervals added in the review (`08_intervals.R`)

**8.1** Of the samples the cleared ratio calls positive, the percentage PET positive, and of those
called negative, the percentage PET negative, with bootstrap intervals; the eGFR rate correlation
from the model with trial arm terms, with its bootstrap interval
(`review/review_8_missing_intervals.csv`).

**8.2** The corrected within-participant estimates with intervals scaled by the same reliability
(part A); the eGFR rate correlation on the participants with both conditions known, without the
two conditions (part C); and with the visit effect of section 7.1, with the within-visit residual
correlation (part D); each resample refits the specification that fitted the full sample
(`review/review_9_missing_intervals.csv`).

Reported in eMethods 8.9 (review_8), eTables 18, 19, 25 and 29.

## Stage 9. Formatted tables (`09_tables.R`)

**9.1** Tables 1 to 5 (`manuscript/Table1.csv` to `Table5.csv`). **9.2** The supplementary tables
`manuscript/S1.csv` to `S18.csv` and `S26.csv`. **9.3** Race and ethnicity as self-reported,
diagnosis and statin use at the first draw, eGFR below 60 and 45, and anaemia by the WHO 2024
thresholds, for Table 1 (`review/review_4_table1_additions.csv`). **9.4** `manuscript/S20.csv` with
the direction measure, then plain-language labels in every formatted table.

Which formatted file belongs to which table:

| Paper | Formatted file | Values from |
|---|---|---|
| Table 1 | Table1.csv, review_4 | stages 1, 2 |
| Tables 2 to 5 | Table2.csv to Table5.csv | T1, T2, T4, T5 with T5b and T5c |
| eTables 1 to 12 | S1 to S12 | as in stages 1 to 4 |
| eTable 13 | checks/*_strings.csv | stages 1, 7 |
| eTables 14 to 17 | S13 to S16 | T3, T2, T6 with T4b, T6b |
| eTable 18 | S17 | T5, R2; intervals of the corrected estimates in review_9 |
| eTable 19 | S18 | T5c, T5b; review_3, review_3b, review_8 |
| eTable 20 | T7b_thresholds.csv | stage 7 |
| eTables 21, 22 | S19, S20 | T7 |
| eTables 23, 24 | S21, S22 | T8, T9 |
| eTable 25 | S23 | T10; intervals in review_9 |
| eTables 26, 27, 28 | S24, S25, S26 | T4b, T5d, T1 |
| eTable 29 | S27 | T11; interval without the conditions in review_9 |
| eTables 30, 31, 32 | S28, S29, S30 | T13a and T13b, T13e, T13f |
| eTable 33 | review_5 | stage 5 |
| eTable 34 | S32 | T13h |
| eTables 35 to 40 | review_7, review_1, review_1b, review_2, review_6, review_6b | stages 3, 4, 7 |

## Stage 10. Figures (`10_figures.R`)

Figure 1, the between-platform differences from T2 with the equivalence margin. Figure 2, the eGFR
spline for p-tau217, p-tau181, NfL and GFAP on every platform and in both cohorts, refitted on the
rows of T1b with its pointwise 95% band (C2N p-tau217 by the censored spline). `eFigure1` is each
curve in a panel of its own, eFigure 2 of the supplement. eFigure 1 of the supplement, the causal
diagram, is drawn by hand and is not produced by the code.
