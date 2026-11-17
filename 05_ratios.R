# =============================================================================
# 05_ratios.R
# Stage 5. Research question 3: the %p-tau217 ratio and its two components; the
# two post hoc Fujirebio ratios and the Fujirebio reporting floor.
#   5.1  the %p-tau217 ratio and the C2N detection limit (Table 4, eTables 16, 26)
#   5.2  the ratio and body mass index in kg/m2 (Results)
#   5.3  the Fujirebio floor, the post hoc ratios and clinical units
#        (eTables 30 to 32, 34)
#   5.4  the cleared ratio corrected for body mass index, against PET (eTable 33)
# Writes paper1_out/tables/T4*.csv, T12, T13*.csv, review/review_5_*.csv.
# =============================================================================

# #############################################################################
# Section 5.1. The %p-tau217 ratio and its components
# Developed as 05_rq3_ratio.R version 1.1.0, whole
# #############################################################################
local({
# =============================================================================
# 05_rq3_ratio.R
# Research question 3: why does the %p-tau217 ratio resist kidney function, and
# does it resist the other factors as well?
#
# C2N reports three things on the same sample: phosphorylated tau217, the
# non-phosphorylated form of the same peptide, and the ratio between them. The
# ratio is 100 x p-tau217 / np-tau217, which was checked against the file and
# holds to 0.005. On the log scale the ratio is therefore the numerator minus
# the denominator, so for least squares on a common sample the effect of a
# factor on the ratio is exactly the difference between its effects on the two.
# Here the numerator and the ratio are censored regressions, and the numerator
# has two more samples than the ratio (the two with no np-tau217 value), so the
# identity holds only approximately in the estimates below (for eGFR -4.48
# against -4.68, for body mass index -9.25 against -8.83).
#
# That gives a direct test of the mechanism. A factor acting on the whole tau
# pool, by clearance or by dilution, moves numerator and denominator together
# and cancels in the ratio. A factor acting on the phosphorylated form alone
# does not cancel. The equivalence test says which of the two has been shown,
# rather than leaving a small estimate to be read as cancellation.
#
# Detection limits are handled measure by measure.
#   p-tau217         reported as half the limit of 1.30 pg/mL in 424 of 1,179
#                    samples, so the censored regression is the main analysis
#   np-tau217        has no censoring to model. The two samples at its limit
#                    carry no value in the file (the laboratory comment is
#                    "LOD"), so they are missing and drop out of the np-tau217
#                    and ratio models; every sample with a value is used
#   the ratio        never flagged itself, but its numerator carries the limit,
#                    so it is censored at 100 x 1.30 / np-tau217, which differs
#                    from sample to sample
#
# False discovery rate: Benjamini-Hochberg within the family of eMethods 9, the
# three measures by all factors of one role. The verdict 'does not cancel'
# needs q < 0.05; the unadjusted p is kept beside it. Each handling in T4b is a
# sensitivity analysis and is adjusted within itself, never pooled with T4.
#
# Reads    paper1_out/adni_visit.rds
# Writes   paper1_out/tables/
#            T4_ratio_components.csv   each factor against numerator, denominator
#                                      and ratio, side by side
#            T4b_ratio_floor.csv       the three ways of handling the limit, each
#                                      with its interval, p, q and verdict
#
# Needs    survival. Takes under a minute.
# =============================================================================
SCRIPT_VERSION <- "1.1.0"
cat("05_rq3_ratio.R version", SCRIPT_VERSION, "\n")

if (!exists("PAPER1_VERSION")) source("00_setup.R")

d <- readRDS(file.path(OUT_DIR, "adni_visit.rds"))
FAC <- FACTORS %>% filter(paste0("z_", var) %in% names(d))

NUM   <- "C2N_plasma_ptau217"
DEN   <- "C2N_plasma_nptau217"
RATIO <- "C2N_plasma_ptau217_ratio"
if (!all(c(NUM, DEN, RATIO) %in% names(d)))
  stop("The C2N p-tau217 family is not in adni_visit.rds.")

# The ratio is the numerator over the denominator, so its limit moves with the
# denominator. Below the limit the reported numerator is a substituted value,
# which makes the reported ratio a substituted value too.
d$ratio_limit <- 100 * C2N_LOD / exp(d[[DEN]])
d$ratio_censored <- coalesce(d$c2n_at_floor, FALSE)

say("Samples at the p-tau217 detection limit: %d of %d (%.2f%%)",
    sum(d$ratio_censored), nrow(d), 100 * mean(d$ratio_censored))

# ---- The censored regression on the rows it actually uses ---------------------
# fit_tobit is given only the complete cases of the model (outcome, covariates,
# cluster, and for the ratio its limit), so that the numbers of samples, of
# censored samples and of people it reports are those of the fitted sample.
# The censoring indicator is the laboratory's flag for both measures; for the
# numerator it agrees exactly with the values (every flagged value is 0.65 and
# no unflagged value is at or below 1.30).
fit_censored <- function(formula, data, term, ratio = FALSE, cluster = "RID") {
  vars <- unique(c(all.vars(formula), cluster, "ratio_censored", if (ratio) "ratio_limit"))
  dd <- as.data.frame(data)
  dd <- dd[stats::complete.cases(dd[, vars, drop = FALSE]), , drop = FALSE]
  rownames(dd) <- NULL
  r <- fit_tobit(formula, dd, term,
                 censor_at = if (ratio) dd$ratio_limit else C2N_LOD,
                 cluster = cluster, censored = dd$ratio_censored)
  tibble(n_obs = nrow(dd), n_censored = sum(dd$ratio_censored),
         n_people = dplyr::n_distinct(dd[[cluster]]),
         b = r$b, se = r$se, lo = r$lo, hi = r$hi, p = r$p)
}

# The verdict on a ratio row, from its estimate, standard error and adjusted q.
# Cancels: the 90% interval lies inside the equivalence margin. Does not cancel:
# the association survives the false discovery rate adjustment. Otherwise
# neither has been shown.
ratio_verdict <- function(is_ratio, equivalent, q) dplyr::case_when(
  !is_ratio            ~ NA_character_,
  is.na(q)             ~ NA_character_,
  equivalent           ~ "cancels in the ratio",
  q < 0.05             ~ "does not cancel",
  TRUE                 ~ "not settled either way")

# ---- Each factor against the three measures ---------------------------------
T4 <- purrr::map_dfr(seq_len(nrow(FAC)), function(j) {
  zv <- paste0("z_", FAC$var[j])
  f  <- function(y) as.formula(sprintf("`%s` ~ %s + %s", y, zv, COVS_ADNI))

  num <- fit_censored(f(NUM), d, zv)
  den <- fit_cl(f(DEN), d, zv)
  rat <- fit_censored(f(RATIO), d, zv, ratio = TRUE)

  bind_rows(
    tibble(measure = "p-tau217 (numerator)", method = "censored regression",
           n_obs = num$n_obs, n_censored = num$n_censored, n_people = num$n_people,
           b = num$b, se = num$se, lo = num$lo, hi = num$hi, p = num$p),
    tibble(measure = "np-tau217 (denominator)", method = "least squares, every quantified sample",
           n_obs = den$n_obs, n_censored = NA_integer_, n_people = den$n_people,
           b = den$b, se = den$se, lo = den$lo, hi = den$hi, p = den$p),
    tibble(measure = "%p-tau217 (ratio)", method = "censored regression",
           n_obs = rat$n_obs, n_censored = rat$n_censored, n_people = rat$n_people,
           b = rat$b, se = rat$se, lo = rat$lo, hi = rat$hi, p = rat$p)
  ) %>%
    mutate(factor = FAC$label[j], var = FAC$var[j], role = FAC$role[j], .before = 1)
}) %>%
  bind_cols(purrr::map_dfr(seq_len(nrow(.)), function(i) equivalence(.$b[i], .$se[i]))) %>%
  group_by(role) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(pct_diff = pct(b), pct_lo = pct(lo), pct_hi = pct(hi),
         verdict = ratio_verdict(measure == "%p-tau217 (ratio)", equivalent, q))
write_csv(T4, file.path(TAB_DIR, "T4_ratio_components.csv"))

say("Each factor on the numerator, the denominator and the ratio, %% per SD:")
print(as.data.frame(
  T4 %>% transmute(factor = sub(" \\(.*", "", factor), measure,
                   estimate = sprintf("%.2f (%.2f to %.2f)", pct_diff, pct_lo, pct_hi)) %>%
    pivot_wider(names_from = measure, values_from = estimate)), row.names = FALSE)

say("What the ratio does with each factor:")
print(as.data.frame(
  T4 %>% filter(measure == "%p-tau217 (ratio)") %>%
    transmute(factor = sub(" \\(.*", "", factor), role,
              ratio = sprintf("%.2f (%.2f to %.2f)", pct_diff, pct_lo, pct_hi),
              p = signif(p, 2), q = signif(q, 2), verdict)),
  row.names = FALSE)

# ---- How much the detection limit matters ------------------------------------
# The same three measures under each way of handling the limit, so that the
# choice of method can be seen rather than assumed, each with its interval, p,
# q and, for the ratio, the verdict under the same rule as T4.
#   as_delivered    least squares on the values as supplied, 0.65 included
#   floor_excluded  least squares with the samples below the limit dropped.
#                   Excluding means dropping the sample, not just its p-tau217
#                   value, because that is what the earlier analyses did and
#                   what the comparison is against.
#   censored        the main analysis, exactly as in T4: censored regression
#                   for the numerator and the ratio, and least squares for
#                   np-tau217, which has no censoring to model.
# Each handling is its own false discovery rate family per factor role.
d_excl <- d %>% filter(!ratio_censored)

T4b_long <- purrr::map_dfr(seq_len(nrow(FAC)), function(j) {
  zv <- paste0("z_", FAC$var[j])
  f  <- function(y) as.formula(sprintf("`%s` ~ %s + %s", y, zv, COVS_ADNI))
  purrr::map_dfr(c(NUM, DEN, RATIO), function(y) {
    as_is <- fit_cl(f(y), d, zv)
    excl  <- fit_cl(f(y), d_excl, zv)
    cens  <- if (y == DEN) fit_cl(f(y), d, zv) else fit_censored(f(y), d, zv, ratio = (y == RATIO))
    bind_rows(
      as_is %>% mutate(handling = "as_delivered"),
      excl  %>% mutate(handling = "floor_excluded"),
      cens  %>% mutate(handling = "censored")) %>%
      transmute(factor = FAC$label[j], var = FAC$var[j], role = FAC$role[j],
                measure = gsub("_plasma", "", y), handling, n_obs, b, se, lo, hi, p)
  })
}) %>%
  bind_cols(purrr::map_dfr(seq_len(nrow(.)), function(i) equivalence(.$b[i], .$se[i]) %>%
                             select(equivalent))) %>%
  group_by(handling, role) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(pct = pct(b), lo = pct(lo), hi = pct(hi),
         verdict = ratio_verdict(measure == gsub("_plasma", "", RATIO), equivalent, q))

T4b <- T4b_long %>%
  select(factor, var, role, measure, handling, n_obs, pct, lo, hi, p, q, verdict) %>%
  pivot_wider(names_from = handling, values_from = c(n_obs, pct, lo, hi, p, q, verdict),
              names_glue = "{.value}_{handling}") %>%
  transmute(factor, var, role, measure,
            n_all = n_obs_as_delivered,
            # The samples remaining after the exclusion, not the number dropped.
            n_after_exclusion = n_obs_floor_excluded,
            n_censored_fit = n_obs_censored,
            method_censored = if_else(measure == gsub("_plasma", "", DEN),
                                      "least squares (no censoring to model)",
                                      "censored regression"),
            pct_as_delivered, lo_as_delivered, hi_as_delivered,
            p_as_delivered, q_as_delivered, verdict_as_delivered,
            pct_floor_excluded, lo_floor_excluded, hi_floor_excluded,
            p_floor_excluded, q_floor_excluded, verdict_floor_excluded,
            pct_censored, lo_censored, hi_censored,
            p_censored, q_censored, verdict_censored)
write_csv(T4b, file.path(TAB_DIR, "T4b_ratio_floor.csv"))

# The censored handling is the main analysis, so it has to reproduce T4.
chk <- T4b %>% select(var, measure, pct_censored, q_censored, verdict_censored) %>%
  inner_join(T4 %>% mutate(measure = rep(gsub("_plasma", "", c(NUM, DEN, RATIO)), nrow(FAC))) %>%
               select(var, measure, pct_diff, q, verdict), by = c("var", "measure"))
if (max(abs(chk$pct_censored - chk$pct_diff)) > 1e-8 || max(abs(chk$q_censored - chk$q)) > 1e-8 ||
    !identical(chk$verdict_censored, chk$verdict))
  stop("The censored handling of T4b does not reproduce T4.")

say("The ratio under each handling of the limit, primary factors:")
print(as.data.frame(
  T4b %>% filter(measure == gsub("_plasma", "", RATIO), role == "primary") %>%
    transmute(factor = sub(" \\(.*", "", factor),
              as_delivered = sprintf("%.2f (%.2f to %.2f) q=%.3f %s", pct_as_delivered,
                                     lo_as_delivered, hi_as_delivered, q_as_delivered,
                                     verdict_as_delivered),
              floor_excluded = sprintf("%.2f (%.2f to %.2f) q=%.3f %s", pct_floor_excluded,
                                       lo_floor_excluded, hi_floor_excluded, q_floor_excluded,
                                       verdict_floor_excluded),
              censored = sprintf("%.2f (%.2f to %.2f) q=%.3f %s", pct_censored,
                                 lo_censored, hi_censored, q_censored, verdict_censored))),
  row.names = FALSE)

say("np-tau217 has no detection limit of its own. All rows against the samples")
say("kept when p-tau217 was quantified:")
print(as.data.frame(
  T4b %>% filter(measure == "C2N_nptau217", var %in% c("egfr", "bmi", "l_ck", "hb")) %>%
    transmute(factor = sub(" \\(.*", "", factor), n_all, n_after_exclusion,
              all_rows = round(pct_as_delivered, 2),
              floor_excluded = round(pct_floor_excluded, 2))), row.names = FALSE)

say("Written to %s", TAB_DIR)
say("T4 %d rows, T4b %d rows", nrow(T4), nrow(T4b))
})


# #############################################################################
# Section 5.2. The ratio and body mass index per 5 kg/m2
# Developed as 12_comorbidity_and_labels.R version 1.0.4, section 2d
# The header lines of 12_comorbidity_and_labels.R that this section needs are repeated here.
# #############################################################################
local({
SCRIPT_VERSION <- "1.0.4"
cat("12_comorbidity_and_labels.R version", SCRIPT_VERSION, "\n")

if (!exists("PAPER1_VERSION")) source("00_setup.R")
for (pkg in c("survival", "nlme")) if (!requireNamespace(pkg, quietly = TRUE))
  stop("The ", pkg, " package is needed for step 12.")

short <- function(x) gsub("_plasma", "", x)
short_factor <- function(x) sub(" \\(.*", "", x)
fmt <- function(e, lo, hi) ifelse(is.na(e), "", sprintf("%.2f (%.2f to %.2f)", e, lo, hi))

adni <- readRDS(file.path(OUT_DIR, "adni_visit.rds"))

# =============================================================================
# 2d. The ratio's body mass index association in kg/m2
# =============================================================================
# Table 4 gives the %p-tau217 association per standard deviation of body mass
# index. For a clinical reader it is put here on the scale of body mass index
# itself, with no new model: the log-scale coefficient and its interval from
# T4_ratio_components.csv are rescaled by the standard deviation that step 01
# used to z-score body mass index, which is that of every ADNI draw with a
# value (zf in 00_setup.R).
say("2d. The ratio and body mass index, per kg/m2")
T4r <- read_csv(file.path(TAB_DIR, "T4_ratio_components.csv"), show_col_types = FALSE)
rat <- T4r %>% filter(var == "bmi", measure == "%p-tau217 (ratio)")
bmi_sd <- stats::sd(adni$bmi, na.rm = TRUE)
if (nrow(rat) == 1 && is.finite(bmi_sd)) {
  steps <- tibble(contrast = c("1 standard deviation higher", "2 standard deviations higher",
                               "5 kg/m2 higher"),
                  sds = c(1, 2, 5 / bmi_sd))
  T12 <- steps %>% mutate(bmi_sd_kg_m2 = bmi_sd,
                          pct = pct(rat$b * sds), lo = pct(rat$lo * sds), hi = pct(rat$hi * sds))
  write_csv(T12, file.path(TAB_DIR, "T12_ratio_bmi.csv"))
  say("  standard deviation of body mass index %.2f kg/m2", bmi_sd)
  print(as.data.frame(T12 %>% transmute(contrast, ratio = fmt(pct, lo, hi))), row.names = FALSE)
} else say("  the ratio row for body mass index was not found in T4_ratio_components.csv; skipped")
})


# #############################################################################
# Section 5.3. The Fujirebio floor, the post hoc ratios and clinical units
# Developed as 13_fujirebio_ratios_and_floor.R version 1.0.1
# The reclassification under a fixed shift (T13g, S31), no longer reported, is removed.
# #############################################################################
local({
# =============================================================================
# 13_fujirebio_ratios_and_floor.R
# Analyses added after the main results had been seen (eMethods 14). Nothing
# here changes a primary or secondary result, and every output file says in its
# name, and in an analysis_type column, which of three kinds it is.
#
#   SENSITIVITY  (re-estimates published results under another handling)
#     1  Fujirebio p-tau217 has a reporting floor: 47 of 1,179 values are
#        exactly 0.030 pg/mL, against 1 to 7 at each neighbouring value, and 42
#        of the 47 carry the laboratory comment "Below Calibrator Range" in the
#        raw plasma file. Step 01 records comments but carries them forward only
#        for C2N, so adni_visit.rds has no Fujirebio flag; the flags are counted
#        here from the raw file, which is read only. The main analysis (steps 03
#        and 04) treats the 47 values as measurements. Here all 47 values
#        are left-censored at 0.030, as C2N values below its limit are:
#          1a  the Fujirebio p-tau217 level models
#          1b  the three p-tau217 pair differences that include Fujirebio
#          1c  the MIMIC model with the floor values set to missing
#          1d  the Deming slopes with the floor values set aside
#
#   POST HOC     (new questions, reported as such)
#     2  Two ratios decomposed as the %p-tau217 ratio is in step 05:
#          the FDA-cleared Lumipulse ratio, Fujirebio p-tau217 over Fujirebio
#            amyloid-beta 42 (510(k) K242706)
#          the ratio of Mondesert 2026 (Neurology), Fujirebio p-tau217 over
#            Roche amyloid-beta 42
#        On the log scale a ratio is the numerator minus the denominator, so a
#        factor that acts on both alike cancels. The numerator is left-censored
#        at the floor, and the ratio at 0.030 over that sample's denominator.
#        Components and ratio are fitted on the same rows. With sensitivity
#        analyses of their own: values as supplied; samples below the cleared
#        measuring interval (0.047 pg/mL) excluded; participants with MCI or
#        dementia at the first draw, the cleared intended use.
#     3  The cleared categories (negative 0.00370 or below, indeterminate
#        0.00371 to 0.00737, positive 0.00738 or above) against amyloid PET,
#        as counts, for eMethods 8.9 and section 8.1. (The reclassification
#        under a 1 or 2 SD shift, T13g and S31, is no longer reported and was
#        removed on reorganisation.)
#
#   RE-EXPRESSION (existing estimates in clinical units; no new model)
#     4  creatine kinase per doubling, eGFR per 15 mL/min/1.73 m2 lower, and
#        body mass index per 5 kg/m2 for the two new ratios
#
# False discovery rate (eMethods 9). Each sensitivity analysis is adjusted
# within itself and never pooled with the analysis it checks. Each post hoc
# ratio is its own family, the three measures by all factors of one role, as
# the %p-tau217 ratio is in step 05; each of its sensitivity analyses is
# adjusted within itself.
#
# Checks, each of which stops the step if it fails:
#   the least squares Fujirebio level models reproduce T1_levels.csv
#   the pair differences with the Fujirebio floor not censored reproduce
#     T2_platform_difference.csv
#   the MIMIC model with the floor untouched reproduces T3 (when lavaan is here)
#   for least squares on common rows, ratio = numerator minus denominator
#   the standardised factors are exactly (x - mean) / SD, so that the clinical
#     units use the SD the models used
#
# Reads    paper1_out/adni_visit.rds, a4_screen.rds
#          ADNI FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.csv, read only, for the
#            laboratory comments on Fujirebio p-tau217
#          paper1_out/tables/T1_levels.csv, T2_platform_difference.csv,
#            T3_latent_direct_effects.csv
# Writes   paper1_out/tables/
#            T13a_sensitivity_fujirebio_floor_levels.csv
#            T13b_sensitivity_fujirebio_floor_differences.csv
#            T13c_sensitivity_fujirebio_floor_mimic.csv
#            T13d_sensitivity_fujirebio_floor_slopes.csv
#            T13e_posthoc_fujirebio_ratios.csv
#            T13f_posthoc_ratio_sensitivity.csv
#            T13h_clinical_units.csv
#          paper1_out/checks/fujirebio_floor.csv, cleared_ratio_calibration.csv
#          paper1_out/manuscript/
#            S28_sensitivity_fujirebio_floor.csv
#            S29_posthoc_fujirebio_ratios.csv
#            S30_posthoc_ratio_sensitivity.csv
#            S32_clinical_units.csv
#
# Needs    survival; lavaan for 1c (skipped with a message if absent).
# Takes a few minutes.
# =============================================================================
# Changes in 1.0.1: the laboratory's "Below Calibrator Range" comments on
# Fujirebio p-tau217 are counted from the raw plasma file. Version 1.0.0 looked
# for a flag column that does not exist and reported 0 flagged values. No
# estimate changes; the censoring rule (a value at 0.030) is the same.
SCRIPT_VERSION <- "1.0.1"
cat("13_fujirebio_ratios_and_floor.R version", SCRIPT_VERSION, "\n")

if (!exists("PAPER1_VERSION")) source("00_setup.R")
if (!requireNamespace("survival", quietly = TRUE)) stop("The survival package is needed for step 13.")
HAVE_LAVAAN <- requireNamespace("lavaan", quietly = TRUE)

MS_DIR <- file.path(OUT_DIR, "manuscript")
dir.create(MS_DIR, showWarnings = FALSE, recursive = TRUE)

d  <- readRDS(file.path(OUT_DIR, "adni_visit.rds"))
a4 <- readRDS(file.path(OUT_DIR, "a4_screen.rds"))
T1 <- read_csv(file.path(TAB_DIR, "T1_levels.csv"), show_col_types = FALSE)
T2 <- read_csv(file.path(TAB_DIR, "T2_platform_difference.csv"), show_col_types = FALSE)
FAC <- FACTORS %>% filter(paste0("z_", var) %in% names(d))

C2N_P217   <- "C2N_plasma_ptau217"
FUJI_P217  <- "Fuji_plasma_ptau217"
ALZ_P217   <- "AlzPath_plasma_ptau217"
JAN_P217   <- "Janssen_plasma_ptau217"
FUJI_AB42  <- "Fuji_plasma_Ab42"
ROCHE_AB42 <- "Roche_plasma_Ab42"
stopifnot(all(c(C2N_P217, FUJI_P217, ALZ_P217, JAN_P217, FUJI_AB42, ROCHE_AB42) %in% names(d)))

# The floor, and the cleared test's lower limit for p-tau217 (K242706).
FUJI_FLOOR     <- 0.030
CLEARED_P217_LOWER <- 0.047
CUT_NEGATIVE   <- 0.00370   # ratio at or below: negative
CUT_POSITIVE   <- 0.00738   # ratio at or above: positive

PLATFORM <- c(C2N = "C2N", Fuji = "Fujirebio", AlzPath = "ALZpath", Janssen = "Janssen",
              Roche = "Roche", QX = "Quanterix")
platform_of <- function(m) unname(PLATFORM[sub("_.*", "", m)])
short_factor <- function(x) sub(" \\(.*", "", x)
fmt <- function(e, lo, hi) ifelse(is.na(e), "", sprintf("%.2f (%.2f to %.2f)", e, lo, hi))
# Two significant figures below 0.1, keeping a trailing zero (0.050, 0.0020);
# two decimals from 0.1 up; <0.001 below that.
fmt_p <- function(x) vapply(x, function(v) {
  if (is.na(v)) return("")
  if (v < 0.001) return("<0.001")
  if (v >= 0.1) return(sprintf("%.2f", v))
  formatC(signif(v, 2), digits = 2, format = "fg", flag = "#")
}, character(1))
put <- function(x, name, what) {
  write_csv(x, file.path(MS_DIR, paste0(name, ".csv")))
  say("  %-36s %-44s %3d rows", paste0(name, ".csv"), what, nrow(x))
}

fuji_at_floor <- function(data) !is.na(data[[FUJI_P217]]) & data[[FUJI_P217]] <= log(FUJI_FLOOR) + 1e-8
c2n_at_floor  <- function(data) coalesce(as.logical(data$c2n_at_floor), FALSE)

# ---- The floor itself ---------------------------------------------------------
fv <- round(exp(d[[FUJI_P217]]), 3)
# The laboratory's own comments, from the raw plasma file (read only).
pl_raw <- read_source(adni_file("plasma"), "plasma")
fr <- pl_raw %>% filter(PLASMA_BIOMARKER == "Fuji_plasma_ptau217")
fr_value <- suppressWarnings(readr::parse_number(fr$TESTVALUE))
fr_flag  <- grepl("below calibrator range", tolower(coalesce(fr$COMMENTS, "")), fixed = TRUE)
floor_chk <- tibble(
  item = c("Fujirebio p-tau217 values", "at 0.030 pg/mL", "people with a value at 0.030",
           "at 0.031", "at 0.032", "at 0.033",
           "below 0.047 pg/mL, the cleared lower limit", "people with a value below 0.047",
           "raw file rows for Fujirebio p-tau217",
           "flagged Below Calibrator Range in the raw file",
           "flagged and at 0.030 pg/mL"),
  value = c(sum(!is.na(fv)), sum(fv == 0.030, na.rm = TRUE),
            n_distinct(d$RID[fuji_at_floor(d)]),
            sum(fv == 0.031, na.rm = TRUE), sum(fv == 0.032, na.rm = TRUE),
            sum(fv == 0.033, na.rm = TRUE),
            sum(fv < CLEARED_P217_LOWER, na.rm = TRUE),
            n_distinct(d$RID[!is.na(fv) & fv < CLEARED_P217_LOWER]),
            nrow(fr), sum(fr_flag), sum(fr_flag & !is.na(fr_value) & abs(fr_value - FUJI_FLOOR) < 1e-9)))
write_csv(floor_chk, file.path(CHECK_DIR, "fujirebio_floor.csv"))
say("The Fujirebio p-tau217 floor:")
print(as.data.frame(floor_chk), row.names = FALSE)
if (min(fv, na.rm = TRUE) != FUJI_FLOOR)
  stop("The lowest Fujirebio p-tau217 value is not 0.030; the floor has changed. Look before going on.")

# ---- Estimation helpers --------------------------------------------------------
# The complete cases of a model, so that the counts are those of the fitted rows.
complete_rows <- function(data, formula, extra = character()) {
  vars <- unique(c(all.vars(formula), "RID", extra))
  dd <- as.data.frame(data)
  dd <- dd[stats::complete.cases(dd[, vars, drop = FALSE]), , drop = FALSE]
  rownames(dd) <- NULL
  dd
}
# A regression on an outcome known only as an interval: exact where lower equals
# upper, left-censored where lower is NA, right-censored where upper is NA.
# Standard errors are clustered by person, as in fit_tobit.
fit_interval <- function(rhs, data, lower, upper, term, cluster = "RID") {
  dd <- as.data.frame(data)
  dd$.l <- lower; dd$.u <- upper
  keep <- !(is.na(dd$.l) & is.na(dd$.u))
  dd <- dd[keep, , drop = FALSE]
  f <- stats::as.formula(sprintf("survival::Surv(.l, .u, type = 'interval2') ~ %s + cluster(%s)",
                                 rhs, cluster))
  fit <- survival::survreg(f, data = dd, dist = "gaussian")
  s <- summary(fit)$table
  n_cens <- sum(is.na(dd$.l) | is.na(dd$.u) | dd$.l != dd$.u)
  if (!term %in% rownames(s))
    return(tibble(n_obs = nrow(dd), n_censored = n_cens, n_people = n_distinct(dd[[cluster]]),
                  b = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
  b <- s[term, 1]; se <- s[term, 2]
  tibble(n_obs = nrow(dd), n_censored = n_cens, n_people = n_distinct(dd[[cluster]]),
         b = b, se = se, lo = b - 1.96 * se, hi = b + 1.96 * se,
         p = 2 * stats::pnorm(-abs(b / se)))
}
lvl_formula <- function(y, zv, covs = COVS_ADNI) as.formula(sprintf("`%s` ~ %s + %s", y, zv, covs))

# =============================================================================
# 1a. Fujirebio p-tau217 levels, floor censored (sensitivity)
# =============================================================================
say("1a. Fujirebio p-tau217 levels with the floor left-censored")
T13a <- purrr::map_dfr(seq_len(nrow(FAC)), function(j) {
  zv <- paste0("z_", FAC$var[j])
  f  <- lvl_formula(FUJI_P217, zv)
  dd <- complete_rows(d, f)
  ls <- fit_cl(f, dd, zv)
  cz <- fit_tobit(f, dd, zv, censor_at = FUJI_FLOOR, censored = fuji_at_floor(dd))
  tibble(factor = FAC$label[j], var = FAC$var[j], role = FAC$role[j],
         n_obs = nrow(dd), n_people = n_distinct(dd$RID), n_censored = sum(fuji_at_floor(dd)),
         b_ls = ls$b, b = cz$b, se = cz$se, lo = cz$lo, hi = cz$hi, p = cz$p)
}) %>%
  group_by(role) %>% mutate(q = p.adjust(p, "BH")) %>% ungroup() %>%
  mutate(pct_diff = pct(b), pct_lo = pct(lo), pct_hi = pct(hi),
         analysis_type = "sensitivity analysis",
         fdr_family = "within this sensitivity analysis, factors of one role")
main_lvl <- T1 %>% filter(cohort == "ADNI", measure == "Fuji_ptau217") %>%
  select(var, main_b = b, main_pct = pct_diff, main_lo = pct_lo, main_hi = pct_hi, main_q = q)
T13a <- T13a %>% left_join(main_lvl, by = "var")
if (anyNA(T13a$main_b) || max(abs(T13a$b_ls - T13a$main_b)) > 1e-8)
  stop("The least squares Fujirebio level models do not reproduce T1_levels.csv.")
say("  check: least squares reproduces T1 for all %d factors", nrow(T13a))
write_csv(T13a, file.path(TAB_DIR, "T13a_sensitivity_fujirebio_floor_levels.csv"))

# =============================================================================
# 1b. The pair differences that include Fujirebio, floor censored (sensitivity)
# =============================================================================
# The difference is platform 1 minus platform 2, in the order of step 04.
#   platform 1 below its limit: the difference is at most (limit 1 - value 2),
#     left-censored. For C2N the limit is 1.30, not the reported 0.65; for
#     Fujirebio the reported 0.030 is the floor itself.
#   platform 2 (Fujirebio) at its floor: the difference is at least
#     (value 1 - log 0.030), right-censored.
#   both: the sample says nothing about the difference and is left out.
say("1b. p-tau217 pair differences that include Fujirebio")
FUJI_PAIRS <- list(c(C2N_P217, FUJI_P217), c(FUJI_P217, ALZ_P217), c(FUJI_P217, JAN_P217))
pair_label <- function(pr) paste(platform_of(pr[1]), "minus", platform_of(pr[2]))
pair_raw   <- function(pr) paste(gsub("_plasma", "", pr[1]), "minus", gsub("_plasma", "", pr[2]))

fit_pair <- function(data, pr, zv, censor_fuji) {
  f  <- as.formula(sprintf("`%s` ~ %s + %s", pr[1], zv, COVS_ADNI))
  dd <- complete_rows(data, f, extra = pr[2])
  y1 <- dd[[pr[1]]]; y2 <- dd[[pr[2]]]; diff <- y1 - y2
  c1 <- if (pr[1] == C2N_P217) c2n_at_floor(dd) else if (censor_fuji) fuji_at_floor(dd) else rep(FALSE, nrow(dd))
  c2 <- if (pr[2] == FUJI_P217 && censor_fuji) fuji_at_floor(dd) else rep(FALSE, nrow(dd))
  lim1 <- if (pr[1] == C2N_P217) log(C2N_LOD) else log(FUJI_FLOOR)
  lower <- ifelse(c1, NA_real_, ifelse(c2, y1 - log(FUJI_FLOOR), diff))
  upper <- ifelse(c2, NA_real_, ifelse(c1, lim1 - y2, diff))
  both  <- c1 & c2
  rhs <- paste(zv, "+", COVS_ADNI)
  r <- fit_interval(rhs, dd, lower, upper, zv)
  r %>% mutate(n_both_censored_dropped = sum(both),
               method = if (any(c1 | c2)) "censored regression" else "maximum likelihood, no censoring")
}

T13b <- purrr::map_dfr(FUJI_PAIRS, function(pr) purrr::map_dfr(seq_len(nrow(FAC)), function(j) {
  zv <- paste0("z_", FAC$var[j])
  chk <- fit_pair(d, pr, zv, censor_fuji = FALSE)
  sen <- fit_pair(d, pr, zv, censor_fuji = TRUE)
  # The check fit uses the same estimator with the floor not censored, so a
  # verdict that changes between it and the sensitivity fit is due to the
  # censoring, and one that changes between it and T2 is due to the estimator.
  bind_cols(tibble(pair = pair_label(pr), pair_raw = pair_raw(pr),
                   factor = FAC$label[j], var = FAC$var[j], role = FAC$role[j],
                   b_check = chk$b, se_check = chk$se,
                   equivalent_check = equivalence(chk$b, chk$se)$equivalent),
            sen, equivalence(sen$b, sen$se))
})) %>%
  group_by(role) %>% mutate(q = p.adjust(p, "BH")) %>% ungroup() %>%
  mutate(verdict = case_when(equivalent ~ "equivalent",
                             !is.na(q) & q < 0.05 ~ "differ",
                             TRUE ~ "no difference detected"),
         pct_diff = pct(b), pct_lo = pct(lo), pct_hi = pct(hi),
         analysis_type = "sensitivity analysis",
         fdr_family = "within this sensitivity analysis, the three pairs by factors of one role")
main_diff <- T2 %>% filter(sample == "main") %>%
  transmute(pair_raw = pair, var, main_b = b, main_pct = pct_diff, main_lo = pct_lo,
            main_hi = pct_hi, main_q = p_fdr, main_verdict = verdict)
T13b <- T13b %>% left_join(main_diff, by = c("pair_raw", "var"))
if (anyNA(T13b$main_b) || max(abs(T13b$b_check - T13b$main_b)) > 1e-4)
  stop("The pair differences without the Fujirebio floor censored do not reproduce T2_platform_difference.csv.")
say("  check: without the floor censored, every difference reproduces T2 (largest gap %.1e)",
    max(abs(T13b$b_check - T13b$main_b)))
write_csv(T13b, file.path(TAB_DIR, "T13b_sensitivity_fujirebio_floor_differences.csv"))
say("  equivalence, primary factors: T2 (least squares) / same estimator, floor not censored / floor censored")
print(as.data.frame(T13b %>% filter(role == "primary") %>%
  transmute(pair, factor = short_factor(factor),
            T2 = main_verdict == "platforms agree to within the margin",
            estimator_only = equivalent_check, censored = equivalent,
            hi90 = round(hi90, 2), lo90 = round(lo90, 2))), row.names = FALSE)

# =============================================================================
# 1c. The MIMIC model, floor values set to missing (sensitivity)
# =============================================================================
# mimic_omnibus is step 04's, version 1.4.1, unchanged.
mimic_omnibus <- function(data, ms, zv) {
  if (!HAVE_LAVAAN || length(ms) < 3) return(NULL)
  dd <- data %>%
    select(RID, all_of(ms), all_of(zv), age_at_visit, female, e4, centiloid, dx, plasma_time) %>%
    filter(!is.na(.data[[zv]]), !is.na(centiloid), !is.na(age_at_visit), !is.na(e4),
           !is.na(female), !is.na(dx), !is.na(plasma_time))
  if (nrow(dd) < 50) return(NULL)
  ind <- paste0("y", seq_along(ms))
  names(dd)[match(ms, names(dd))] <- ind
  dx_d <- model.matrix(~ dx, dd)[, -1, drop = FALSE]
  colnames(dx_d) <- paste0("dx", seq_len(ncol(dx_d)))
  dd <- bind_cols(dd, as.data.frame(dx_d))
  dd$x <- dd[[zv]]
  covs <- paste(c("age_at_visit", "female", "e4", "centiloid", colnames(dx_d), "plasma_time"),
                collapse = " + ")
  base <- sprintf("eta =~ %s\n eta ~ x + %s", paste(ind, collapse = " + "), covs)
  fitm <- function(model) {
    w <- character(0)
    f <- tryCatch(withCallingHandlers(
      lavaan::sem(model, data = dd, estimator = "MLR", missing = "fiml",
                  cluster = "RID", fixed.x = TRUE),
      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }),
      error = function(e) { w <<- c(w, paste("error:", conditionMessage(e))); NULL })
    if (!is.null(f) && !isTRUE(lavaan::lavInspect(f, "converged"))) w <- c(w, "did not converge")
    list(fit = f, warnings = w)
  }
  f0    <- fitm(base)
  f_all <- fitm(paste0(base, "\n", paste(sprintf("%s ~ x", ind[-1]), collapse = "\n")))
  if (is.null(f0$fit)) return(NULL)
  fm <- function(f, what) if (is.null(f)) NA_real_ else as.numeric(lavaan::fitMeasures(f, what))
  lrt <- if (!is.null(f_all$fit))
    tryCatch(as.data.frame(lavaan::lavTestLRT(f0$fit, f_all$fit)), error = function(e) NULL) else NULL
  tibble(n_obs = as.integer(lavaan::nobs(f0$fit)), n_people = n_distinct(dd$RID),
         omnibus_chisq = if (is.null(lrt)) NA_real_ else lrt$`Chisq diff`[2],
         omnibus_df = if (is.null(lrt)) NA_real_ else lrt$`Df diff`[2],
         omnibus_p = if (is.null(lrt)) NA_real_ else lrt$`Pr(>Chisq)`[2],
         fit_warnings = paste(unique(c(f0$warnings, f_all$warnings)), collapse = " | "))
}
MIMIC_P217 <- c(FUJI_P217, ALZ_P217, JAN_P217)
T13c <- NULL
if (HAVE_LAVAAN) {
  say("1c. MIMIC model for p-tau217 with the Fujirebio floor set to missing")
  T3 <- read_csv(file.path(TAB_DIR, "T3_latent_direct_effects.csv"), show_col_types = FALSE) %>%
    filter(specification == "main", analyte == "pTau217")
  d_m <- d; d_m[[FUJI_P217]][fuji_at_floor(d_m)] <- NA_real_
  T13c <- purrr::map_dfr(seq_len(nrow(FAC)), function(j) {
    zv <- paste0("z_", FAC$var[j])
    chk <- mimic_omnibus(d, MIMIC_P217, zv)
    sen <- mimic_omnibus(d_m, MIMIC_P217, zv)
    if (is.null(sen)) return(NULL)
    bind_cols(tibble(factor = FAC$label[j], var = FAC$var[j], role = FAC$role[j],
                     check_chisq = if (is.null(chk)) NA_real_ else chk$omnibus_chisq), sen)
  }) %>%
    group_by(role) %>% mutate(omnibus_q = p.adjust(omnibus_p, "BH")) %>% ungroup() %>%
    left_join(T3 %>% select(var, main_chisq = omnibus_chisq, main_p = omnibus_p, main_q = omnibus_q),
              by = "var") %>%
    mutate(analysis_type = "sensitivity analysis",
           fdr_family = "within this sensitivity analysis, factors of one role")
  if (anyNA(T13c$main_chisq) || max(abs(T13c$check_chisq - T13c$main_chisq)) > 1e-6)
    stop("The MIMIC model with the floor untouched does not reproduce T3.")
  say("  check: with the floor untouched the omnibus test reproduces T3")
  write_csv(T13c, file.path(TAB_DIR, "T13c_sensitivity_fujirebio_floor_mimic.csv"))
} else say("1c skipped: lavaan is not installed.")

# =============================================================================
# 1d. Deming slopes with both floors set aside (sensitivity)
# =============================================================================
say("1d. Deming slopes for the Fujirebio pairs with the floor values set aside")
d_s <- d
d_s[[FUJI_P217]][fuji_at_floor(d_s)] <- NA_real_
d_s[[C2N_P217]][c2n_at_floor(d_s)]   <- NA_real_
T13d <- purrr::map_dfr(FUJI_PAIRS, function(pr) {
  main <- T2 %>% filter(sample == "main", pair == pair_raw(pr)) %>% slice(1)
  tibble(pair = pair_label(pr),
         samples = sum(stats::complete.cases(d_s[[pr[1]]], d_s[[pr[2]]])),
         slope_050 = deming_slope(d_s[[pr[2]]], d_s[[pr[1]]], 0.5),
         slope_100 = deming_slope(d_s[[pr[2]]], d_s[[pr[1]]], 1),
         slope_200 = deming_slope(d_s[[pr[2]]], d_s[[pr[1]]], 2),
         main_slope_050 = main$slope_050, main_slope_100 = main$slope_100,
         main_slope_200 = main$slope_200, analysis_type = "sensitivity analysis")
})
write_csv(T13d, file.path(TAB_DIR, "T13d_sensitivity_fujirebio_floor_slopes.csv"))
print(as.data.frame(T13d %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

# =============================================================================
# 2. The two ratios (post hoc)
# =============================================================================
say("2. The FDA-cleared ratio and the Mondesert ratio")
RATIOS <- list(
  cleared   = list(num = FUJI_P217, den = FUJI_AB42,
                   label = "Fujirebio p-tau217 / Fujirebio amyloid-beta 42 (FDA-cleared)"),
  mondesert = list(num = FUJI_P217, den = ROCHE_AB42,
                   label = "Fujirebio p-tau217 / Roche amyloid-beta 42 (Mondesert 2026)"))

ratio_verdict <- function(equivalent, q) dplyr::case_when(
  is.na(q)   ~ NA_character_,
  equivalent ~ "cancels in the ratio",
  q < 0.05   ~ "does not cancel",
  TRUE       ~ "inconclusive")

# One factor, one ratio: numerator, denominator and ratio on the same rows.
# censor = TRUE is the main handling; FALSE treats the floor as a measurement.
decompose <- function(data, r, zv, censor = TRUE) {
  f0 <- lvl_formula(r$num, zv)
  dd <- complete_rows(data, f0, extra = r$den)
  dd$.ratio <- dd[[r$num]] - dd[[r$den]]
  cens <- fuji_at_floor(dd)
  use_c <- censor && any(cens)
  num <- if (use_c) fit_tobit(f0, dd, zv, censor_at = FUJI_FLOOR, censored = cens) else fit_cl(f0, dd, zv)
  den <- fit_cl(lvl_formula(r$den, zv), dd, zv)
  rat <- if (use_c) fit_tobit(lvl_formula(".ratio", zv), dd, zv,
                              censor_at = FUJI_FLOOR / exp(dd[[r$den]]), censored = cens)
         else fit_cl(lvl_formula(".ratio", zv), dd, zv)
  meth <- if (use_c) "censored regression" else "least squares"
  bind_rows(
    tibble(measure = "numerator", platform_measure = "Fujirebio p-tau217", method = meth,
           b = num$b, se = num$se, lo = num$lo, hi = num$hi, p = num$p),
    tibble(measure = "denominator",
           platform_measure = paste(platform_of(r$den), "amyloid-beta 42"), method = "least squares",
           b = den$b, se = den$se, lo = den$lo, hi = den$hi, p = den$p),
    tibble(measure = "ratio", platform_measure = "ratio", method = meth,
           b = rat$b, se = rat$se, lo = rat$lo, hi = rat$hi, p = rat$p)) %>%
    mutate(n_obs = nrow(dd), n_people = n_distinct(dd$RID), n_censored = if (use_c) sum(cens) else 0L)
}

run_ratios <- function(data, handling, censor) {
  purrr::map_dfr(names(RATIOS), function(k) purrr::map_dfr(seq_len(nrow(FAC)), function(j) {
    zv <- paste0("z_", FAC$var[j])
    decompose(data, RATIOS[[k]], zv, censor = censor) %>%
      mutate(ratio_key = k, ratio_label = RATIOS[[k]]$label, handling = handling,
             factor = FAC$label[j], var = FAC$var[j], role = FAC$role[j], .before = 1)
  })) %>%
    bind_cols(purrr::map_dfr(seq_len(nrow(.)), function(i) equivalence(.$b[i], .$se[i]))) %>%
    group_by(ratio_key, role) %>% mutate(q = p.adjust(p, "BH")) %>% ungroup() %>%
    mutate(pct_diff = pct(b), pct_lo = pct(lo), pct_hi = pct(hi),
           verdict = if_else(measure == "ratio", ratio_verdict(equivalent, q), NA_character_))
}

# How the ratio behaves, from the point estimates (descriptive), given only
# where the ratio does not cancel, since a pattern read off a null estimate
# would mean nothing: magnified when
# the components move in opposite directions and the ratio keeps the sign of the
# numerator; reversed (over-correction) when the ratio takes the other sign;
# attenuated when it keeps the sign and is smaller than the numerator.
pattern_of <- function(bn, bd, br) dplyr::case_when(
  is.na(br) ~ NA_character_,
  sign(br) != sign(bn) ~ "reversed",
  abs(br) > abs(bn) ~ "magnified",
  TRUE ~ "attenuated")

T13e <- run_ratios(d, "main: numerator floor left-censored", censor = TRUE)
patt <- T13e %>% select(ratio_key, var, measure, b) %>%
  pivot_wider(names_from = measure, values_from = b) %>%
  transmute(ratio_key, var, pattern = pattern_of(numerator, denominator, ratio))
T13e <- T13e %>% left_join(patt, by = c("ratio_key", "var")) %>%
  mutate(pattern = if_else(measure == "ratio" & verdict == "does not cancel", pattern, NA_character_),
         analysis_type = "post hoc analysis",
         fdr_family = "this ratio: numerator, denominator and ratio by all factors of one role")
write_csv(T13e, file.path(TAB_DIR, "T13e_posthoc_fujirebio_ratios.csv"))

# Sensitivity analyses of the post hoc ratios, each adjusted within itself.
d_in  <- d %>% filter(!is.na(.data[[FUJI_P217]]), exp(.data[[FUJI_P217]]) >= CLEARED_P217_LOWER - 1e-9)
d_sym <- d %>% filter(dx %in% c("MCI", "Dementia"))
T13f <- bind_rows(
  run_ratios(d,     "values as supplied, floor treated as measured", censor = FALSE),
  run_ratios(d_in,  "samples below 0.047 pg/mL excluded", censor = FALSE),
  run_ratios(d_sym, "MCI or dementia at the first draw", censor = TRUE)) %>%
  mutate(analysis_type = "sensitivity analysis of a post hoc analysis",
         fdr_family = "within this handling and ratio, measures by factors of one role")
write_csv(T13f, file.path(TAB_DIR, "T13f_posthoc_ratio_sensitivity.csv"))

# The least squares decomposition on common rows is an identity.
ident <- T13f %>% filter(handling == "values as supplied, floor treated as measured") %>%
  select(ratio_key, var, measure, b) %>% pivot_wider(names_from = measure, values_from = b) %>%
  mutate(gap = ratio - (numerator - denominator))
if (max(abs(ident$gap)) > 1e-8) stop("Least squares ratio is not numerator minus denominator.")
say("  check: least squares ratio = numerator - denominator (largest gap %.1e)", max(abs(ident$gap)))

say("The two ratios, primary factors, %% per SD (numerator / denominator / ratio, verdict, pattern):")
print(as.data.frame(
  T13e %>% filter(role == "primary") %>%
    mutate(est = sprintf("%.2f (%.2f to %.2f)", pct_diff, pct_lo, pct_hi)) %>%
    group_by(ratio_key, factor = short_factor(factor)) %>%
    summarise(numerator = est[measure == "numerator"], denominator = est[measure == "denominator"],
              ratio_est = est[measure == "ratio"], q = signif(q[measure == "ratio"], 2),
              verdict = verdict[measure == "ratio"], pattern = pattern[measure == "ratio"],
              n = n_obs[1], censored = n_censored[1], .groups = "drop")), row.names = FALSE)

# =============================================================================
# 3. Reclassification against the cleared cutoffs (post hoc, descriptive)
# =============================================================================
# The ratio is reported to five decimals, so it is rounded before it is
# categorised, as the device's report would be.
cleared_category <- function(r) {
  r5 <- round(r, 5)
  dplyr::case_when(is.na(r5) ~ NA_character_,
                   r5 <= CUT_NEGATIVE ~ "negative",
                   r5 >= CUT_POSITIVE ~ "positive",
                   TRUE ~ "indeterminate")
}
cr <- d %>% filter(!is.na(.data[[FUJI_P217]]), !is.na(.data[[FUJI_AB42]])) %>%
  mutate(r = exp(.data[[FUJI_P217]] - .data[[FUJI_AB42]]), cat0 = cleared_category(r))

# Agreement of the cleared categories with amyloid PET: a calibration check that
# the cutoffs can be read on these values, written to checks/ only.
calib <- purrr::map_dfr(c(20, 24), function(cl) cr %>% filter(!is.na(centiloid)) %>%
  mutate(pet = if_else(centiloid >= cl, "PET positive", "PET negative")) %>%
  count(centiloid_threshold = cl, cat0, pet) %>%
  pivot_wider(names_from = pet, values_from = n, values_fill = 0))
write_csv(calib, file.path(CHECK_DIR, "cleared_ratio_calibration.csv"))
say("3. The cleared categories against amyloid PET (samples):")
print(as.data.frame(calib), row.names = FALSE)

say("  baseline categories, all samples: %s",
    paste(names(table(cr$cat0)), table(cr$cat0), collapse = ", "))
# =============================================================================
# 4. Clinical units (re-expression)
# =============================================================================
say("4. Clinical units")
# The factors were standardised by zf, (x - mean) / SD with the sign of the
# label, within each cohort. That SD is checked here so that the conversion
# uses the SD the models used.
sd_used <- function(data, v) {
  x <- data[[v]]; z <- data[[paste0("z_", v)]]
  ok <- !is.na(x) & !is.na(z)
  s <- stats::sd(x, na.rm = TRUE)
  slope <- unname(stats::coef(stats::lm(z[ok] ~ x[ok]))[2])
  if (abs(abs(slope) - 1 / s) > 1e-8) stop("z_", v, " is not (x - mean) / SD; the conversion would be wrong.")
  s
}
SD <- list(ADNI = c(l_ck = sd_used(d, "l_ck"), egfr = sd_used(d, "egfr"), bmi = sd_used(d, "bmi")),
           A4   = c(l_ck = sd_used(a4, "l_ck"), egfr = sd_used(a4, "egfr")))
say("  SD used: ADNI log creatine kinase %.4f (%.2f doublings), eGFR %.2f, BMI %.2f;",
    SD$ADNI[["l_ck"]], SD$ADNI[["l_ck"]] / log(2), SD$ADNI[["egfr"]], SD$ADNI[["bmi"]])
say("           A4 and LEARN log creatine kinase %.4f, eGFR %.2f", SD$A4[["l_ck"]], SD$A4[["egfr"]])
unit_k <- function(cohort, v) {
  s <- SD[[if (cohort == "ADNI") "ADNI" else "A4"]]
  switch(v, l_ck = log(2) / s[["l_ck"]], egfr = 15 / s[["egfr"]], bmi = 5 / s[["bmi"]])
}
unit_name <- c(l_ck = "per doubling of creatine kinase", egfr = "per 15 mL/min/1.73 m2 lower eGFR",
               bmi = "per 5 kg/m2 higher body mass index")
conv <- function(tab) tab %>% rowwise() %>% mutate(k = unit_k(cohort, var)) %>% ungroup() %>%
  mutate(unit = unit_name[var], pct_per_sd = pct(b), pct_unit = pct(k * b),
         lo_unit = pct(k * lo), hi_unit = pct(k * hi))
T13h <- bind_rows(
  T1 %>% filter(var %in% c("l_ck", "egfr")) %>%
    transmute(source = "level (T1)", cohort, measure, var, b, lo, hi),
  T2 %>% filter(sample == "main", var == "l_ck", analyte %in% c("pTau217", "pTau181")) %>%
    transmute(source = "difference (T2)", cohort = "ADNI", measure = pair, var, b, lo, hi),
  T13e %>% filter(measure == "ratio", var %in% c("egfr", "bmi", "l_ck")) %>%
    transmute(source = "post hoc ratio (T13e)", cohort = "ADNI", measure = ratio_label, var, b, lo, hi)) %>%
  conv() %>%
  mutate(analysis_type = "re-expression of an existing estimate")
write_csv(T13h, file.path(TAB_DIR, "T13h_clinical_units.csv"))

# =============================================================================
# 5. Formatted tables for the supplement
# =============================================================================
say("5. Formatted tables")
label_measure <- function(m) {
  m <- gsub("^C2N_ptau217_ratio$", "C2N %p-tau217", m)
  m <- gsub("^C2N_nptau217$", "C2N non-phosphorylated tau217", m)
  m <- gsub("^pTau181$", "p-tau181", gsub("^Ab42_40$", "amyloid-beta 42/40", m))
  m <- gsub("ptau217", "p-tau217", gsub("ptau181", "p-tau181", m))
  m <- gsub("Abeta42_Abeta40|Ab42_Ab40", "amyloid-beta 42/40", m)
  m <- gsub("Abeta42|Ab42", "amyloid-beta 42", gsub("Abeta40|Ab40", "amyloid-beta 40", m))
  for (k in names(PLATFORM)) m <- gsub(paste0("\\b", k, "_"), paste0(PLATFORM[[k]], " "), m)
  m <- gsub("_", " ", m)
  # In a pair the analyte is named once, in its own column, as in eTable 11.
  gsub("^(\\w+) [^ ]+ minus (\\w+) [^ ]+$", "\\1 minus \\2", m)
}

S28 <- bind_rows(
  T13a %>% transmute(analysis = "level", measure = "Fujirebio p-tau217",
                     factor = short_factor(factor), role, samples = n_obs, censored = n_censored,
                     main_estimate = fmt(main_pct, main_lo, main_hi), main_q = fmt_p(main_q),
                     main_verdict = "",
                     censored_estimate = fmt(pct_diff, pct_lo, pct_hi), censored_q = fmt_p(q),
                     censored_verdict = ""),
  T13b %>% transmute(analysis = "difference", measure = pair,
                     factor = short_factor(factor), role, samples = n_obs, censored = n_censored,
                     main_estimate = fmt(main_pct, main_lo, main_hi), main_q = fmt_p(main_q),
                     main_verdict = recode(main_verdict,
                                           "platforms agree to within the margin" = "equivalent",
                                           "platforms differ" = "differ"),
                     censored_estimate = fmt(pct_diff, pct_lo, pct_hi), censored_q = fmt_p(q),
                     censored_verdict = verdict))
put(S28, "S28_sensitivity_fujirebio_floor", "Fujirebio floor, levels and differences")

S29 <- T13e %>%
  mutate(est = fmt(pct_diff, pct_lo, pct_hi)) %>%
  group_by(ratio = ratio_label, factor = short_factor(factor), role) %>%
  summarise(samples = n_obs[1], censored = n_censored[1],
            numerator = est[measure == "numerator"],
            denominator = est[measure == "denominator"],
            ratio_estimate = est[measure == "ratio"],
            ratio_q = fmt_p(q[measure == "ratio"]),
            verdict = verdict[measure == "ratio"],
            pattern = pattern[measure == "ratio"], .groups = "drop") %>%
  mutate(analysis_type = "post hoc") %>%
  arrange(ratio, desc(role == "primary"))
put(S29, "S29_posthoc_fujirebio_ratios", "the two ratios and their components")

S30 <- T13f %>% filter(measure == "ratio") %>%
  transmute(ratio = ratio_label, handling, factor = short_factor(factor), role,
            samples = n_obs, censored = n_censored,
            ratio_estimate = fmt(pct_diff, pct_lo, pct_hi), ratio_q = fmt_p(q), verdict,
            analysis_type = "sensitivity analysis of a post hoc analysis")
put(S30, "S30_posthoc_ratio_sensitivity", "the two ratios under three handlings")

S32 <- T13h %>% transmute(source, cohort = if_else(cohort == "ADNI", "ADNI", "A4 and LEARN"),
                          measure = if_else(grepl("^post hoc", source), measure, label_measure(measure)),
                          unit, per_sd = sprintf("%.2f", pct_per_sd),
                          per_unit = fmt(pct_unit, lo_unit, hi_unit),
                          analysis_type = "re-expression")
put(S32, "S32_clinical_units", "creatine kinase, eGFR and BMI in clinical units")

say("Fujirebio floor, primary factors: least squares (main) against censored (sensitivity):")
print(as.data.frame(S28 %>% filter(role == "primary") %>%
  select(analysis, measure, factor, main_estimate, main_verdict, censored_estimate, censored_q,
         censored_verdict)), row.names = FALSE)
if (!is.null(T13c)) {
  say("MIMIC omnibus test for p-tau217, floor untouched (main) and set to missing:")
  print(as.data.frame(T13c %>% filter(role == "primary") %>%
    transmute(factor = short_factor(factor), main_chisq = round(main_chisq, 2),
              main_q = signif(main_q, 2), chisq = round(omnibus_chisq, 2), df = omnibus_df,
              q = signif(omnibus_q, 2))), row.names = FALSE)
}
say("Clinical units for creatine kinase and eGFR, selected rows:")
print(as.data.frame(S32 %>% filter(grepl("p-tau", measure) | grepl("ratio", source))), row.names = FALSE)
say("Step 13 finished. Written to %s and %s", TAB_DIR, MS_DIR)
})


# #############################################################################
# Section 5.4. The cleared ratio corrected for body mass index, against amyloid PET (eTable 33)
# Developed as 14_review_checks.R version 1.1.0, section 5
# The header lines of 14_review_checks.R that this section needs are repeated here.
# #############################################################################
local({
SCRIPT_VERSION <- "1.1.0"
cat("14_review_checks.R version", SCRIPT_VERSION, "\n")

if (!exists("PAPER1_VERSION")) source("00_setup.R")
for (pkg in c("survival", "nlme", "splines"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("The ", pkg, " package is needed for step 14.")

REV_DIR <- file.path(TAB_DIR, "review")
dir.create(REV_DIR, showWarnings = FALSE, recursive = TRUE)
put_rev <- function(x, name) {
  write_csv(x, file.path(REV_DIR, paste0(name, ".csv")))
  say("  written review/%s.csv, %d rows", name, nrow(x))
}
env_int <- function(name, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (is.na(v) || v < 0) default else v
}
REVIEW_BOOT     <- env_int("PAPER1_REVIEW_BOOT", 1000)
REVIEW_ARM_BOOT <- env_int("PAPER1_REVIEW_ARM_BOOT", 1000)
RESAMPLES <- 1000


adni <- readRDS(file.path(OUT_DIR, "adni_visit.rds"))
read_tab <- function(name) read_csv(file.path(TAB_DIR, paste0(name, ".csv")), show_col_types = FALSE)
T13e <- read_tab("T13e_posthoc_fujirebio_ratios")

C2N_P217  <- "C2N_plasma_ptau217"
P217      <- c("C2N_plasma_ptau217", "Fuji_plasma_ptau217", "AlzPath_plasma_ptau217",
               "Janssen_plasma_ptau217")
P181      <- c("Roche_plasma_ptau181", "QX_plasma_ptau181")
short     <- function(x) gsub("_plasma", "", x)
fmt       <- function(e, lo, hi) ifelse(is.na(e), "", sprintf("%.2f (%.2f to %.2f)", e, lo, hi))
# The covariates with terms removed, for a subsample in which a term is
# constant or for a model in which it is the modifier.
covs_without <- function(drop, covs = COVS_ADNI) {
  terms <- trimws(strsplit(covs, "+", fixed = TRUE)[[1]])
  paste(setdiff(terms, drop), collapse = " + ")
}
# Centiloid as a natural spline with 3 degrees of freedom in place of the
# linear term.
covs_spline_cl <- function(covs = COVS_ADNI) {
  terms <- trimws(strsplit(covs, "+", fixed = TRUE)[[1]])
  terms[terms == "centiloid"] <- "splines::ns(centiloid, df = 3)"
  paste(terms, collapse = " + ")
}

# Participant resampling for any statistic. fun takes a data frame and returns
# a named numeric vector; the interval is the 2.5th to 97.5th percentile over
# B resamples of participants (RID).
boot_pct <- function(d, fun, B = RESAMPLES, seed = 1) {
  est <- fun(d)
  idx <- split(seq_len(nrow(d)), d$RID)
  ids <- names(idx)
  set.seed(seed)
  res <- vapply(seq_len(B), function(b) {
    take <- unlist(idx[sample(ids, length(ids), replace = TRUE)], use.names = FALSE)
    out <- tryCatch(fun(d[take, , drop = FALSE]), error = function(e) rep(NA_real_, length(est)))
    if (length(out) != length(est)) rep(NA_real_, length(est)) else out
  }, numeric(length(est)))
  draws <- if (is.matrix(res)) t(res) else matrix(res, ncol = 1)
  tibble(stat = names(est), est = unname(est),
         lo = unname(apply(draws, 2, stats::quantile, 0.025, na.rm = TRUE)),
         hi = unname(apply(draws, 2, stats::quantile, 0.975, na.rm = TRUE)),
         usable = unname(colSums(!is.na(draws))))
}

# =============================================================================
# 5. The cleared ratio corrected for each sample's own body mass index
# =============================================================================
# Each sample's cleared ratio is corrected to the mean body mass index of the
# ADNI draws, by the fitted body mass index coefficient of the ratio (T13e):
# r_adjusted = r x exp(-b x z_bmi). The three categories (negative,
# indeterminate, positive) before and after are compared with amyloid PET at a
# Centiloid of 24, and 20 as a sensitivity cut, separately in PET-positive and
# PET-negative samples:
#   PET positive: % moving up a category, % moving down, and their difference;
#   PET negative: % moving down, % moving up, and their difference;
#   % of PET-positive samples called positive and of PET-negative samples
#   called negative, and % indeterminate, before and after.
# The coefficient is estimated in the same samples and held fixed in the
# resamples, so this is descriptive and somewhat optimistic.
say("5. The cleared ratio corrected for body mass index, against amyloid PET")
FUJI_P217 <- "Fuji_plasma_ptau217"; FUJI_AB42 <- "Fuji_plasma_Ab42"
CUT_NEGATIVE <- 0.00370; CUT_POSITIVE <- 0.00738
cleared_category <- function(r) {
  r5 <- round(r, 5)
  dplyr::case_when(is.na(r5) ~ NA_character_, r5 <= CUT_NEGATIVE ~ "negative",
                   r5 >= CUT_POSITIVE ~ "positive", TRUE ~ "indeterminate")
}
cat_rank <- function(x) match(x, c("negative", "indeterminate", "positive"))
bmi_b <- T13e %>% filter(ratio_key == "cleared", measure == "ratio", var == "bmi") %>% pull(b)
if (length(bmi_b) != 1) stop("The cleared ratio estimate for body mass index was not found in T13e.")
cr <- adni %>% filter(!is.na(.data[[FUJI_P217]]), !is.na(.data[[FUJI_AB42]])) %>%
  mutate(r = exp(.data[[FUJI_P217]] - .data[[FUJI_AB42]]), cat0 = cleared_category(r))
nri_stats <- function(d) {
  up <- cat_rank(d$cat1) > cat_rank(d$cat0)
  dn <- cat_rank(d$cat1) < cat_rank(d$cat0)
  pos <- d$pet; neg <- !d$pet
  c(pet_pos_up = 100 * mean(up[pos]), pet_pos_down = 100 * mean(dn[pos]),
    pet_pos_net = 100 * (mean(up[pos]) - mean(dn[pos])),
    pet_neg_down = 100 * mean(dn[neg]), pet_neg_up = 100 * mean(up[neg]),
    pet_neg_net = 100 * (mean(dn[neg]) - mean(up[neg])),
    pos_called_pos_before = 100 * mean(d$cat0[pos] == "positive"),
    pos_called_pos_after  = 100 * mean(d$cat1[pos] == "positive"),
    neg_called_neg_before = 100 * mean(d$cat0[neg] == "negative"),
    neg_called_neg_after  = 100 * mean(d$cat1[neg] == "negative"),
    indeterminate_before  = 100 * mean(d$cat0 == "indeterminate"),
    indeterminate_after   = 100 * mean(d$cat1 == "indeterminate"))
}
R5 <- purrr::map_dfr(c("all samples", "MCI or dementia at the first draw"), function(sub) {
  base <- if (sub == "all samples") cr else cr %>% filter(dx %in% c("MCI", "Dementia"))
  purrr::map_dfr(c(24, 20), function(cut) {
    d <- base %>% filter(!is.na(centiloid), !is.na(z_bmi)) %>%
      mutate(pet = centiloid >= cut, cat1 = cleared_category(r * exp(-bmi_b * z_bmi)))
    boot_pct(d, nri_stats) %>%
      mutate(subsample = sub, centiloid_cut = cut, samples = nrow(d), people = n_distinct(d$RID),
             pet_positive = sum(d$pet), changed = sum(d$cat1 != d$cat0), .before = 1)
  })
})
put_rev(R5, "review_5_bmi_correction_reclassification")
print(as.data.frame(R5 %>% filter(centiloid_cut == 24) %>%
  transmute(subsample, samples, pet_positive, changed, stat, value = fmt(est, lo, hi))), row.names = FALSE)
})
