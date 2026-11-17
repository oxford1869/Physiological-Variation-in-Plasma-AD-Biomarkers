# =============================================================================
# 09_tables.R
# Stage 9. The formatted tables. Nothing is estimated here.
#   9.1  the main tables, Tables 1 to 5
#   9.2  the supplementary tables built from stages 1 to 6 (S1 to S18, S26)
#   9.3  the additions to Table 1: race, ethnicity, diagnosis, statin use, eGFR
#        bands and anaemia
#   9.4  eTable 22 with the direction measure, and plain-language labels in
#        every formatted table
# Writes paper1_out/manuscript/*.csv and review/review_4_table1_additions.csv.
# =============================================================================

# #############################################################################
# Section 9.1. Tables 1 to 5
# Developed as 07_tables.R version 1.0.0, whole
# #############################################################################
local({
# =============================================================================
# 07_tables.R
# The tables for the manuscript, built from what steps 03 to 06 wrote.
#
# Five main tables, one for the sample and one for each research question, each
# carrying the primary factors only. Everything else goes to supplementary.
# Nothing is estimated here: this script formats results that already exist, so
# a table can never disagree with the analysis that produced it.
#
#   Table 1  the people in each sample
#   Table 2  question 1, levels, both cohorts
#   Table 3  question 2, differences between platforms
#   Table 4  question 3, the p-tau217 ratio and its two parts
#   Table 5  question 4, change on change
#
# The supplementary tables are built separately, by 09_supplementary_tables.R.
#
# Numbers are given to two decimal places, as percentage difference per
# standard deviation of the factor with a 95% interval, unless a table says
# otherwise.
#
# Reads    paper1_out/*.rds and paper1_out/tables/*.csv
# Writes   paper1_out/manuscript/Table1.csv to Table5.csv
#
# Takes a few seconds.
# =============================================================================
SCRIPT_VERSION <- "1.0.0"
cat("07_tables.R version", SCRIPT_VERSION, "\n")

if (!exists("PAPER1_VERSION")) source("00_setup.R")

MS_DIR <- file.path(OUT_DIR, "manuscript")
dir.create(MS_DIR, showWarnings = FALSE, recursive = TRUE)

# A result table from an earlier step, or NULL with a warning if that step has
# not been run. Missing input leaves a table out rather than stopping.
grab <- function(name) {
  f <- file.path(TAB_DIR, paste0(name, ".csv"))
  if (!file.exists(f)) {
    say("  %s.csv is not there, so the tables that need it are skipped.", name)
    return(NULL)
  }
  read_csv(f, show_col_types = FALSE)
}
put <- function(x, name, what) {
  if (is.null(x) || nrow(x) == 0) return(invisible(NULL))
  write_csv(x, file.path(MS_DIR, paste0(name, ".csv")))
  say("  %-8s %-58s %3d rows", paste0(name, ".csv"), what, nrow(x))
}

est <- function(e, lo, hi) sprintf("%.2f (%.2f to %.2f)", e, lo, hi)
short_factor <- function(x) sub(" \\(.*", "", x)
PRIMARY <- c("eGFR", "Body mass index", "Creatine kinase", "Haemoglobin")

adni <- readRDS(file.path(OUT_DIR, "adni_visit.rds"))
a4s  <- if (file.exists(file.path(OUT_DIR, "a4_screen.rds")))
  readRDS(file.path(OUT_DIR, "a4_screen.rds")) else NULL
a4l  <- if (file.exists(file.path(OUT_DIR, "a4_long.rds")))
  readRDS(file.path(OUT_DIR, "a4_long.rds")) else NULL

say("Building the tables")

# ---- Table 1. The people -----------------------------------------------------
# One row per characteristic, one column per sample. Continuous measures are
# median and interquartile range, counts are n and percentage.
med_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("")
  sprintf("%.2f (%.2f to %.2f)", median(x), quantile(x, 0.25), quantile(x, 0.75))
}
n_pct <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("")
  sprintf("%d (%.2f%%)", sum(x == 1), 100 * mean(x == 1))
}

adni_person <- adni %>% arrange(RID, visit_no) %>% group_by(RID) %>% slice(1) %>% ungroup()
a4l_person  <- if (!is.null(a4l)) a4l %>% arrange(BID, draw_day) %>%
  group_by(BID) %>% slice(1) %>% ungroup() else NULL

column_for <- function(d, label, id, draws) {
  tibble(characteristic = c(
    "People", "Blood draws", "Age at the first draw, years", "Women",
    "APOE4 carriers", "Centiloid", "eGFR, mL/min/1.73m2", "Body mass index, kg/m2",
    "Creatine kinase, U/L", "Haemoglobin, g/L", "Follow-up, years"),
    value = c(
      as.character(n_distinct(d[[id]])),
      as.character(draws),
      med_iqr(d$age_at_visit),
      n_pct(d$female),
      n_pct(d$e4),
      med_iqr(d$centiloid),
      med_iqr(d$egfr),
      med_iqr(d$bmi),
      if ("l_ck" %in% names(d)) med_iqr(exp(d$l_ck)) else med_iqr(d$ck),
      med_iqr(d$hb),
      if ("plasma_time" %in% names(d)) "" else "")) %>%
    rename(!!label := value)
}

T1_people <- column_for(adni_person, "ADNI", "RID", nrow(adni))
if (!is.null(a4s))
  T1_people <- left_join(T1_people, column_for(a4s, "A4 screening", "BID", nrow(a4s)),
                         by = "characteristic")
if (!is.null(a4l_person))
  T1_people <- left_join(T1_people, column_for(a4l_person, "A4 repeated draws", "BID", nrow(a4l)),
                         by = "characteristic")

# Follow-up is the span of the draws, which only the repeated samples have.
span <- function(d, id, t) {
  s <- d %>% group_by(.data[[id]]) %>% summarise(s = max(.data[[t]]) - min(.data[[t]]), .groups = "drop")
  med_iqr(s$s)
}
T1_people$ADNI[T1_people$characteristic == "Follow-up, years"] <- span(adni, "RID", "plasma_time")
if (!is.null(a4l))
  T1_people$`A4 repeated draws`[T1_people$characteristic == "Follow-up, years"] <-
    span(a4l, "BID", "plasma_time")
put(T1_people, "Table1", "the people in each sample")

# ---- Table 2. Question 1, levels ---------------------------------------------
T1r <- grab("T1_levels")
MAIN_MEASURES <- c("Roche_ptau181", "QX_ptau181", "C2N_ptau217", "AlzPath_ptau217",
                   "Roche_NfL", "QX_NfL", "Roche_GFAP", "QX_GFAP",
                   "C2N_Abeta42_Abeta40", "pTau181", "NfL", "GFAP", "Ab42_40")
if (!is.null(T1r)) {
  # The sample differs a little between factors, because not everyone has every
  # laboratory value, so the column gives the range rather than one number.
  people_range <- function(d) d %>%
    group_by(cohort, measure) %>%
    summarise(people = if (min(n_people) == max(n_people)) as.character(min(n_people))
              else sprintf("%d to %d", min(n_people), max(n_people)), .groups = "drop")
  t2_src <- T1r %>%
    mutate(factor = short_factor(factor)) %>%
    filter(role == "primary", measure %in% MAIN_MEASURES)
  t2 <- t2_src %>%
    transmute(cohort, measure, factor, cell = est(pct_diff, pct_lo, pct_hi)) %>%
    pivot_wider(names_from = factor, values_from = cell) %>%
    left_join(people_range(t2_src), by = c("cohort", "measure")) %>%
    select(cohort, measure, people, any_of(PRIMARY)) %>%
    arrange(cohort, match(measure, MAIN_MEASURES))
  put(t2, "Table2", "question 1, levels, % per SD")


}

# ---- Table 3. Question 2, differences between platforms ----------------------
T2 <- grab("T2_platform_difference")
T3l <- grab("T3_latent_direct_effects")
if (!is.null(T2)) {
  t3_src <- T2 %>% mutate(factor = short_factor(factor)) %>% filter(role == "primary")
  wide <- t3_src %>%
    transmute(analyte, pair, factor, cell = est(pct_diff, pct_lo, pct_hi)) %>%
    pivot_wider(names_from = factor, values_from = cell) %>%
    left_join(t3_src %>% group_by(analyte, pair) %>%
                summarise(people = if (min(n_people) == max(n_people))
                  as.character(min(n_people)) else
                    sprintf("%d to %d", min(n_people), max(n_people)), .groups = "drop"),
              by = c("analyte", "pair"))
  # Which factors this pair agrees on, to within the margin
  agrees <- T2 %>%
    mutate(factor = short_factor(factor)) %>%
    filter(role == "primary", equivalent) %>%
    group_by(analyte, pair) %>%
    summarise(equivalent_for = paste(factor, collapse = ", "), .groups = "drop")
  t3 <- wide %>%
    left_join(agrees, by = c("analyte", "pair")) %>%
    mutate(equivalent_for = coalesce(equivalent_for, "none")) %>%
    select(analyte, pair, people, any_of(PRIMARY), equivalent_for)
  put(t3, "Table3", "question 2, difference between platforms, % per SD")


}



# ---- Table 4. Question 3, the ratio and its two parts ------------------------
T4 <- grab("T4_ratio_components")
if (!is.null(T4)) {
  t4 <- T4 %>%
    mutate(factor = short_factor(factor)) %>%
    filter(role == "primary") %>%
    transmute(factor, measure, n_obs, cell = est(pct_diff, pct_lo, pct_hi),
              verdict = coalesce(verdict, "")) %>%
    pivot_wider(names_from = measure, values_from = c(cell, n_obs, verdict)) %>%
    select(factor,
           `p-tau217 (numerator)` = `cell_p-tau217 (numerator)`,
           `np-tau217 (denominator)` = `cell_np-tau217 (denominator)`,
           `%p-tau217 (ratio)` = `cell_%p-tau217 (ratio)`,
           verdict = `verdict_%p-tau217 (ratio)`) %>%
    arrange(match(factor, PRIMARY))
  put(t4, "Table4", "question 3, the ratio and its two parts, % per SD")
}

# ---- Table 5. Question 4, change on change -----------------------------------
T5  <- grab("T5_hybrid_within_between")
R2  <- grab("R2_change_reliability")
T5b <- grab("T5b_joint_growth")
T5c <- grab("T5c_joint_bootstrap")
if (!is.null(T5) && !is.null(T5b)) {
  joint <- T5b %>% filter(group == "everyone") %>% select(var, slope_correlation)
  if (!is.null(T5c))
    joint <- joint %>%
      left_join(T5c %>% filter(group == "everyone") %>% select(var, lo, hi), by = "var")
  # T5 already carries rel_within from step 06, so R2 is used only when it does not.
  if (!"rel_within" %in% names(T5) && !is.null(R2))
    T5 <- T5 %>% left_join(R2 %>% select(var, rel_within), by = "var")
  joined <- T5 %>%
    left_join(FACTORS %>% select(var, label), by = "var") %>%
    left_join(joint, by = "var")
  joined$rate_cell <- if (all(c("lo", "hi") %in% names(joined)))
    ifelse(is.na(joined$lo), sprintf("%.2f", joined$slope_correlation),
           est(joined$slope_correlation, joined$lo, joined$hi)) else
    sprintf("%.2f", joined$slope_correlation)
  t5 <- joined %>%
    transmute(factor = short_factor(label),
              role = FACTORS$role[match(var, FACTORS$var)],
              people = n_people, draws = n_draws,
              between_people = est(between_pct, between_lo, between_hi),
              within_person = est(within_pct, within_lo, within_hi),
              reliability_of_change = sprintf("%.2f", rel_within),
              rate_correlation = rate_cell) %>%
    arrange(role != "primary", match(factor, PRIMARY))
  put(t5 %>% filter(role == "primary") %>% select(-role), "Table5",
      "question 4, change on change")

}


say("Written to %s", MS_DIR)
say("Every table is formatted from the results of the earlier steps, so rerun")
say("those steps and then this one after any change to the analysis.")
})


# #############################################################################
# Section 9.2. Supplementary tables
# Developed as 09_supplementary_tables.R version 1.1.5, whole
# #############################################################################
local({
# =============================================================================
# 09_supplementary_tables.R
# The supplementary tables, built from what steps 01 to 06 wrote.
#
# The supplement carries every analysis the paper ran, not a selection. The
# argument of the paper is that method choices change answers, so a reader has
# to be able to see each choice and what it gave. Like step 7, this script only
# formats: nothing is estimated here, so a supplementary table cannot disagree
# with the analysis behind it.
#
#   Reproducibility
#     S1   participants and samples, how each was arrived at
#     S2   the data files each cohort was built from, with versions
#   Question 1, levels
#     S3   the exploratory factors
#     S4   the primary factors in one mutually adjusted model
#     S5   corrected for measurement error, with the reliabilities
#     S6   restricted to recent laboratory values
#     S7   creatine kinase in men and women
#     S8   eGFR modelled flexibly
#     S9   amyloid PET as the outcome
#   Question 2, platforms
#     S10  whether platforms are equivalent, counted by analyte and factor
#     S11  every platform pair and every factor
#     S12  creatine kinase in each subsample
#     S13  the latent concentration model
#     S14  the slope correction across three error-variance ratios
#   Question 3, the ratio
#     S15  the detection limit, three ways of handling it
#     S16  platform differences with the floored C2N values set aside
#   Question 4, change
#     S17  the hybrid model and the reliability of a within-person change
#     S18  by trial arm, with bootstrap intervals
#
# Reads    paper1_out/tables/*.csv and paper1_out/checks/*.csv
# Writes   paper1_out/manuscript/S1.csv to S18.csv
#
# Takes a few seconds.
# =============================================================================
SCRIPT_VERSION <- "1.1.5"
cat("09_supplementary_tables.R version", SCRIPT_VERSION, "\n")

if (!exists("PAPER1_VERSION")) source("00_setup.R")

MS_DIR <- file.path(OUT_DIR, "manuscript")
dir.create(MS_DIR, showWarnings = FALSE, recursive = TRUE)

grab <- function(name, dir = TAB_DIR) {
  f <- file.path(dir, paste0(name, ".csv"))
  if (!file.exists(f)) {
    say("  %-34s missing, so the table that needs it is skipped", paste0(name, ".csv"))
    return(NULL)
  }
  read_csv(f, show_col_types = FALSE)
}
put <- function(x, name, what) {
  if (is.null(x) || nrow(x) == 0) return(invisible(NULL))
  write_csv(x, file.path(MS_DIR, paste0(name, ".csv")))
  say("  %-7s %-58s %3d rows", paste0(name, ".csv"), what, nrow(x))
}
est <- function(e, lo, hi) sprintf("%.2f (%.2f to %.2f)", e, lo, hi)
short_factor <- function(x) sub(" \\(.*", "", x)
MAIN_MEASURES <- c("Roche_ptau181", "QX_ptau181", "C2N_ptau217", "AlzPath_ptau217",
                   "Roche_NfL", "QX_NfL", "Roche_GFAP", "QX_GFAP",
                   "C2N_Abeta42_Abeta40", "pTau181", "NfL", "GFAP", "Ab42_40")

say("Building the supplementary tables")

# ---- Reproducibility ---------------------------------------------------------
adni_build <- grab("adni_build", CHECK_DIR)
a4_build   <- grab("a4_build",   CHECK_DIR)
put(bind_rows(
      if (!is.null(adni_build)) adni_build %>% mutate(cohort = "ADNI", .before = 1) else NULL,
      if (!is.null(a4_build))   a4_build   %>% mutate(cohort = "A4 and LEARN", .before = 1) else NULL),
    "S1", "participants and samples, how each was arrived at")

adni_src <- grab("adni_sources", CHECK_DIR)
a4_src   <- grab("a4_sources",   CHECK_DIR)
put(bind_rows(
      if (!is.null(adni_src)) adni_src %>% mutate(cohort = "ADNI", .before = 1) else NULL,
      if (!is.null(a4_src))   a4_src   %>% mutate(cohort = "A4 and LEARN", .before = 1) else NULL) %>%
      select(any_of(c("cohort", "label", "file", "rows", "columns", "modified", "size_mb", "folder"))),
    "S2", "the data files each cohort was built from")

# ---- Question 1 --------------------------------------------------------------
T1r <- grab("T1_levels")
if (!is.null(T1r)) {
  put(T1r %>%
        mutate(factor = short_factor(factor)) %>%
        filter(role == "exploratory", measure %in% MAIN_MEASURES) %>%
        transmute(cohort, measure, factor, n_people, cell = est(pct_diff, pct_lo, pct_hi)) %>%
        pivot_wider(names_from = factor, values_from = cell),
      "S3", "levels, the exploratory factors")

  R1 <- grab("R1_reliability")
  if (!is.null(R1))
    put(T1r %>%
          mutate(factor = short_factor(factor)) %>%
          filter(measure %in% MAIN_MEASURES, !is.na(pct_corrected)) %>%
          arrange(role != "primary", cohort, match(measure, MAIN_MEASURES)) %>%
          # The correction divides by the reliability conditional on the model
          # covariates, which step 03 computes for each row as rel_conditional.
          transmute(cohort, measure, factor, role,
                    reliability_given_covariates = sprintf("%.2f", rel_conditional),
                    uncorrected = est(pct_diff, pct_lo, pct_hi),
                    corrected = est(pct_corrected, pct_corrected_lo, pct_corrected_hi)),
        "S5", "levels corrected for measurement error")
}

# Every level estimate: all biomarkers on every platform, by every factor, in
# both cohorts. Table 2 and eTable 3 show a selection; this table reports
# everything that was estimated, so nothing analysed goes unreported. q is the
# false discovery rate adjusted value within the family of eMethods 9.
if (!is.null(T1r))
  put(T1r %>%
        transmute(cohort, measure, factor = short_factor(factor), role,
                  method = if ("method" %in% names(T1r)) method else NA_character_,
                  n_people, n_obs,
                  estimate = est(pct_diff, pct_lo, pct_hi),
                  p = signif(p, 2),
                  p_false_discovery = if ("q" %in% names(T1r)) signif(q, 2) else NA_real_) %>%
        arrange(cohort, role != "primary", measure, factor),
      "S26", "levels, every biomarker, platform and factor")

T1c <- grab("T1c_mutual_adjustment")
if (!is.null(T1c))
  put(T1c %>%
        left_join(FACTORS %>% select(var, label), by = "var") %>%
        transmute(cohort, measure, factor = short_factor(label), n_people,
                  mutually_adjusted = est(pct_diff, pct_lo, pct_hi)) %>%
        pivot_wider(names_from = factor, values_from = mutually_adjusted),
      "S4", "levels, the primary factors in one model")

T1f <- grab("T1f_lab_recency")
if (!is.null(T1f))
  put(T1f %>%
        mutate(within = ifelse(is.na(within_years), "any",
                               sprintf("%.0f year", within_years))) %>%
        transmute(measure, factor = short_factor(factor), laboratory_value_within = within,
                  n_people, estimate = est(pct_diff, pct_lo, pct_hi)) %>%
        arrange(measure, factor),
      "S6", "levels, restricted to recent laboratory values")

T1e <- grab("T1e_sex_stratified")
if (!is.null(T1e))
  put(T1e %>% transmute(cohort, measure, sex, n_people,
                        creatine_kinase = est(pct_diff, pct_lo, pct_hi)),
      "S7", "creatine kinase in men and women")

T1b <- grab("T1b_egfr_shape")
if (!is.null(T1b))
  put(T1b %>% transmute(cohort, measure, n_obs, samples_used = floor_handling,
                        p_against_a_straight_line = signif(p_curved, 2),
                        # 1.1.5: each difference carries its 95% interval (step 03, 1.3.0)
                        at_eGFR_60 = sprintf("%.2f (%.2f to %.2f)", pct_at_60, pct_at_60_lo, pct_at_60_hi),
                        at_eGFR_45 = sprintf("%.2f (%.2f to %.2f)", pct_at_45, pct_at_45_lo, pct_at_45_hi),
                        at_eGFR_30 = sprintf("%.2f (%.2f to %.2f)", pct_at_30, pct_at_30_lo, pct_at_30_hi)),
      "S8", "eGFR modelled flexibly, difference from eGFR 90")

T1d <- grab("T1d_negative_control")
if (!is.null(T1d))
  put(T1d %>% transmute(cohort, factor = short_factor(factor), role, n_people,
                        centiloid_per_sd = est(centiloid_per_sd, lo, hi),
                        p = signif(p, 2)),
      "S9", "amyloid PET as the outcome")

# ---- Question 2 --------------------------------------------------------------
put(grab("T2b_equivalence_summary"), "S10",
    "whether platforms are equivalent, by analyte and factor")

T2 <- grab("T2_platform_difference")
if (!is.null(T2)) {
  put(T2 %>%
        transmute(analyte, pair, factor = short_factor(factor), role, n_people,
                  method = if ("method" %in% names(T2)) method else "least squares",
                  difference = est(pct_diff, pct_lo, pct_hi),
                  p = signif(p, 2), p_false_discovery = signif(p_fdr, 2),
                  interval_90 = sprintf("%.2f to %.2f", lo90, hi90),
                  margin_pct, verdict),
      "S11", "every platform pair and every factor")

  if (all(c("corrected_050", "corrected_100", "corrected_200") %in% names(T2)))
    put(T2 %>%
          filter(role == "primary") %>%
          transmute(analyte, pair, factor = short_factor(factor),
                    uncorrected = est(pct_diff, pct_lo, pct_hi),
                    ratio_0.50 = est(corrected_050, corrected_lo_050, corrected_hi_050),
                    ratio_1.00 = est(corrected_100, corrected_lo_100, corrected_hi_100),
                    ratio_2.00 = est(corrected_200, corrected_lo_200, corrected_hi_200),
                    slope_at_0.50 = sprintf("%.2f", slope_050),
                    slope_at_1.00 = sprintf("%.2f", slope_100),
                    slope_at_2.00 = sprintf("%.2f", slope_200),
                    slope_estimated_on = slope_sample),
        "S14", "the slope correction across three error-variance ratios")
}

T2c <- grab("T2c_sensitivity")
if (!is.null(T2c))
  put(T2c %>% transmute(sample, analyte, pair, n_people,
                        creatine_kinase = est(pct_diff, pct_lo, pct_hi),
                        p_false_discovery = signif(p_fdr, 2), verdict),
      "S12", "creatine kinase in each subsample")

T3l <- grab("T3_latent_direct_effects")
if (!is.null(T3l))
  # Only the omnibus test is reported. The common effect and the effect on any
  # one platform are not identified without assuming which assay is unbiased;
  # the omnibus test does not depend on that choice (see step 04).
  put(T3l %>% transmute(specification = if ("specification" %in% names(T3l)) specification else "main",
                        analyte,
                        platforms = if ("platforms" %in% names(T3l)) platforms else NA_character_,
                        factor = short_factor(factor), role,
                        n_obs, n_people,
                        chi_square = sprintf("%.2f", omnibus_chisq), df = omnibus_df,
                        omnibus_p = signif(omnibus_p, 2),
                        p_false_discovery = signif(omnibus_q, 2)),
      "S13", "the latent concentration model")

# ---- Question 3 --------------------------------------------------------------
T6  <- grab("T6_floor_sensitivity")
T4b <- grab("T4b_ratio_floor")
s15 <- bind_rows(
  if (!is.null(T6)) T6 %>% transmute(
    measure = "C2N_ptau217", factor = short_factor(factor), role,
    n_samples = n_obs, n_at_the_limit = n_censored,
    as_delivered = sprintf("%.2f", pct_as_delivered),
    samples_excluded = sprintf("%.2f", pct_floor_excluded),
    censored = est(pct_censored, lo_censored, hi_censored)) else NULL,
  # np-tau217 has no detection limit of its own, so its third handling is the
  # least-squares fit on every quantified sample and it has no samples at the
  # limit; the ratio is censored wherever its numerator is at the limit.
  if (!is.null(T4b)) T4b %>% filter(measure != "C2N_ptau217") %>% transmute(
    measure, factor = short_factor(factor), role,
    n_samples = n_all,
    n_at_the_limit = ifelse(measure == "C2N_nptau217", NA_integer_,
                            as.integer(n_all - n_after_exclusion)),
    as_delivered = sprintf("%.2f", pct_as_delivered),
    samples_excluded = sprintf("%.2f", pct_floor_excluded),
    censored = est(pct_censored, lo_censored, hi_censored)) else NULL)
put(s15, "S15", "the detection limit, three ways of handling it")

put(grab("T6b_difference_floor_excluded") %>%
      (function(x) if (is.null(x)) NULL else x %>%
         transmute(analyte, pair, factor = short_factor(factor), role, n_people,
                   difference = est(pct_diff, pct_lo, pct_hi),
                   p_false_discovery = signif(p_fdr, 2), verdict)),
    "S16", "platform differences with the floored C2N values set aside")

# ---- Question 4 --------------------------------------------------------------
T5 <- grab("T5_hybrid_within_between")
R2 <- grab("R2_change_reliability")
if (!is.null(T5))
  put(T5 %>%
        left_join(FACTORS %>% select(var, label, role), by = "var") %>%
        left_join(if (!is.null(R2)) R2 %>% select(var, draws_per_person, span_years) else
                    tibble(var = character()), by = "var") %>%
        transmute(factor = short_factor(label), role, people = n_people, draws = n_draws,
                  values_per_person = sprintf("%.2f", draws_per_person),
                  span_years = sprintf("%.2f", span_years),
                  between_people = est(between_pct, between_lo, between_hi),
                  within_person = est(within_pct, within_lo, within_hi),
                  p_for_the_difference = signif(p_difference, 2),
                  reliability_of_change = sprintf("%.2f", rel_within),
                  within_corrected = ifelse(is.na(within_corrected_pct), "",
                                            sprintf("%.2f", within_corrected_pct))) %>%
        arrange(role != "primary"),
      "S17", "the hybrid model and the reliability of a within-person change")

T5c <- grab("T5c_joint_bootstrap")
if (!is.null(T5c))
  put(T5c %>%
        left_join(FACTORS %>% select(var, label), by = "var") %>%
        transmute(factor = short_factor(label), group, draws, usable,
                  correlation = sprintf("%.2f", correlation),
                  bootstrap_median = sprintf("%.2f", bootstrap_median),
                  interval = sprintf("%.2f to %.2f", lo, hi),
                  minutes = sprintf("%.1f", minutes)),
      "S18", "change on change by trial arm, with bootstrap intervals")

say("Written to %s", MS_DIR)
})


# #############################################################################
# Section 9.3. Additions to Table 1
# Developed as 14_review_checks.R version 1.1.0, section 4
# The check of haemoglobin at the first repeated draw (section 4b), not reported, is removed.
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
a4s  <- readRDS(file.path(OUT_DIR, "a4_screen.rds"))
L    <- readRDS(file.path(OUT_DIR, "a4_long.rds"))
read_tab <- function(name) read_csv(file.path(TAB_DIR, paste0(name, ".csv")), show_col_types = FALSE)

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
# 4. Additions to Table 1
# =============================================================================
# Race and ethnicity as self-reported: ADNIMERGE PTRACCAT and PTETHCAT in ADNI,
# and SUBJINFO RACE and ETHNIC in A4/LEARN, coded White 1, Black or African
# American 2, Asian 58, Native Hawaiian or Other Pacific Islander 79, American
# Indian or Alaskan Native 84, unknown 97, more than one race 100; Hispanic or
# Latino 50, not Hispanic or Latino 56, unknown 97 (A4 derived data dictionary).
# Anaemia is haemoglobin below 130 g/L in men and 120 g/L in women (WHO 2024).
# eGFR is from a single creatinine value (CKD-EPI 2021). Diagnosis and statin
# use are at the first draw, ADNI only: A4/LEARN enrolled cognitively
# unimpaired people and released no medication file.
say("4. Additions to Table 1")
mg <- read_chr(adni_file("merge"))
race_adni <- mg %>%
  transmute(RID = as.integer(RID), race = na_if(str_trim(PTRACCAT), ""),
            eth = na_if(str_trim(PTETHCAT), "")) %>%
  group_by(RID) %>%
  summarise(race = dplyr::first(race[!is.na(race)]), eth = dplyr::first(eth[!is.na(eth)]),
            .groups = "drop") %>%
  mutate(race = case_when(race == "White" ~ "White",
                          race == "Black" ~ "Black or African American",
                          race == "Asian" ~ "Asian",
                          race %in% c("Am Indian/Alaskan", "Hawaiian/Other PI", "More than one") ~
                            "Other or more than one race",
                          TRUE ~ "Unknown or not reported"),
         hispanic = case_when(eth == "Hisp/Latino" ~ 1L, eth == "Not Hisp/Latino" ~ 0L,
                              TRUE ~ NA_integer_))
si <- read_chr(a4_file("subjinfo"))
race_a4 <- si %>%
  transmute(BID, race = suppressWarnings(as.integer(RACE)), eth = suppressWarnings(as.integer(ETHNIC))) %>%
  mutate(race = case_when(race == 1 ~ "White", race == 2 ~ "Black or African American",
                          race == 58 ~ "Asian", race %in% c(79, 84, 100) ~ "Other or more than one race",
                          TRUE ~ "Unknown or not reported"),
         hispanic = case_when(eth == 50 ~ 1L, eth == 56 ~ 0L, TRUE ~ NA_integer_)) %>%
  select(BID, race, hispanic)

RACES <- c("White", "Black or African American", "Asian", "Other or more than one race",
           "Unknown or not reported")
anaemia <- function(hb, female) as.integer(ifelse(female == 1, hb < 120, hb < 130))
n_of <- function(x, known = !is.na(x)) {
  if (!any(known)) return("")
  sprintf("%d (%.2f%%)", sum(x[known] == 1), 100 * mean(x[known] == 1))
}
describe <- function(p, sample_label, adni_only = FALSE) {
  rows <- list(
    tibble(characteristic = paste("Race:", RACES),
           value = vapply(RACES, function(k) n_of(as.integer(p$race == k)), character(1))),
    tibble(characteristic = c("Hispanic or Latino, of those reported", "Ethnicity unknown or not reported"),
           value = c(n_of(p$hispanic), as.character(sum(is.na(p$hispanic))))),
    tibble(characteristic = c("eGFR below 60, of those with eGFR", "eGFR below 45, of those with eGFR",
                              "People with eGFR"),
           value = c(n_of(as.integer(p$egfr < 60)), n_of(as.integer(p$egfr < 45)),
                     as.character(sum(!is.na(p$egfr))))),
    tibble(characteristic = c("Anaemia, of those with haemoglobin", "People with haemoglobin"),
           value = c(n_of(anaemia(p$hb, p$female)), as.character(sum(!is.na(p$hb))))))
  if (adni_only) rows <- c(rows, list(
    tibble(characteristic = paste("Diagnosis at the first draw:", c("CN", "MCI", "Dementia")),
           value = vapply(c("CN", "MCI", "Dementia"), function(k)
             n_of(as.integer(as.character(p$dx) == k), !is.na(p$dx)), character(1))),
    tibble(characteristic = c("Statin at the first draw, of those with medication records",
                              "No medication record"),
           value = c(n_of(p$statin_at_draw), as.character(sum(is.na(p$statin_at_draw)))))))
  bind_rows(rows) %>% mutate(sample = sample_label, people = nrow(p), .before = 1)
}
p_adni <- adni %>% arrange(RID, visit_no) %>% group_by(RID) %>% slice(1) %>% ungroup() %>%
  left_join(race_adni, by = "RID")
p_a4s  <- a4s %>% left_join(race_a4, by = "BID")
p_a4l  <- L %>% arrange(BID, draw_day) %>% group_by(BID) %>% slice(1) %>% ungroup() %>%
  left_join(race_a4, by = "BID")
R4 <- bind_rows(describe(p_adni, "ADNI", adni_only = TRUE),
                describe(p_a4s, "A4 and LEARN screening"),
                describe(p_a4l, "A4 and LEARN repeated draws"))
put_rev(R4, "review_4_table1_additions")
print(as.data.frame(R4 %>% select(sample, characteristic, value)), row.names = FALSE)
})


# #############################################################################
# Section 9.4. eTable 22 with the direction measure, and plain-language labels
# Developed as 12_comorbidity_and_labels.R version 1.0.4, sections 3 and 4
# The header lines of 12_comorbidity_and_labels.R that this section needs are repeated here.
# The labels must be applied after every other formatted table is written, so this section runs last.
# #############################################################################
local({
SCRIPT_VERSION <- "1.0.4"
cat("12_comorbidity_and_labels.R version", SCRIPT_VERSION, "\n")

if (!exists("PAPER1_VERSION")) source("00_setup.R")
for (pkg in c("survival", "nlme")) if (!requireNamespace(pkg, quietly = TRUE))
  stop("The ", pkg, " package is needed for step 12.")

MS_DIR <- file.path(OUT_DIR, "manuscript")
dir.create(MS_DIR, showWarnings = FALSE, recursive = TRUE)
put <- function(x, name, what) {
  if (is.null(x) || nrow(x) == 0) { say("  %s skipped, nothing to write", name); return(invisible(NULL)) }
  write_csv(x, file.path(MS_DIR, paste0(name, ".csv")))
  say("  %-8s %-58s %3d rows", paste0(name, ".csv"), what, nrow(x))
}
short <- function(x) gsub("_plasma", "", x)
short_factor <- function(x) sub(" \\(.*", "", x)
fmt <- function(e, lo, hi) ifelse(is.na(e), "", sprintf("%.2f (%.2f to %.2f)", e, lo, hi))

# =============================================================================
# 3. eTable 22 (S20) with the direction measure
# =============================================================================
say("3. S20 with the direction of disagreement")
T7 <- read_csv(file.path(TAB_DIR, "T7_classification.csv"), show_col_types = FALSE)
put(T7 %>% filter(!(sample == "all samples" & centiloid_cut == 20)) %>%
      transmute(sample, centiloid_cut, pair, pair_samples = n_samples,
                disagree_pct = round(disagree_pct, 2), factor = short_factor(factor),
                model_samples = n_samples_model, model_people = n_people,
                disagreement_per_sd = fmt(disagree_pp, disagree_lo, disagree_hi),
                disagreement_q = signif(disagree_q, 2),
                direction_per_sd = fmt(net_pp, net_lo, net_hi),
                direction_q = signif(net_q, 2)),
    "S20", "classification disagreement and its direction, other thresholds and samples")

# =============================================================================
# 4. Plain-language labels in every formatted table
# =============================================================================
say("4. Relabelling the formatted tables")
PLATFORM <- c(C2N = "C2N", Fuji = "Fujirebio", AlzPath = "ALZpath", Janssen = "Janssen",
              Roche = "Roche", QX = "Quanterix")
MEASURE_LABEL <- c(
  C2N_ptau217 = "C2N p-tau217", Fuji_ptau217 = "Fujirebio p-tau217",
  AlzPath_ptau217 = "ALZpath p-tau217", Janssen_ptau217 = "Janssen p-tau217",
  C2N_nptau217 = "C2N non-phosphorylated tau217", C2N_ptau217_ratio = "C2N %p-tau217",
  Roche_ptau181 = "Roche p-tau181", QX_ptau181 = "Quanterix p-tau181",
  Roche_NfL = "Roche NfL", QX_NfL = "Quanterix NfL",
  Roche_GFAP = "Roche GFAP", QX_GFAP = "Quanterix GFAP",
  C2N_Abeta40 = "C2N amyloid-beta 40", Fuji_Ab40 = "Fujirebio amyloid-beta 40",
  Roche_Ab40 = "Roche amyloid-beta 40", QX_Ab40 = "Quanterix amyloid-beta 40",
  C2N_Abeta42 = "C2N amyloid-beta 42", Fuji_Ab42 = "Fujirebio amyloid-beta 42",
  Roche_Ab42 = "Roche amyloid-beta 42", QX_Ab42 = "Quanterix amyloid-beta 42",
  C2N_Abeta42_Abeta40 = "C2N amyloid-beta 42/40", Fuji_Ab42_Ab40 = "Fujirebio amyloid-beta 42/40",
  Roche_Ab42_Ab40 = "Roche amyloid-beta 42/40", QX_Ab42_Ab40 = "Quanterix amyloid-beta 42/40",
  # A4/LEARN measures and the analyte column
  pTau217 = "p-tau217", pTau181 = "p-tau181", npTau217 = "non-phosphorylated tau217",
  Abeta40 = "amyloid-beta 40", Abeta42 = "amyloid-beta 42", Ab42_40 = "amyloid-beta 42/40")
WHOLE_CELL <- c("A4 screening" = "A4 and LEARN",
                "pTau217 ratio" = "%p-tau217",
                "not settled either way" = "inconclusive",
                # version 1.0.4: the between-platform verdicts in the words the
                # manuscript and eMethods 8.6 use
                "platforms agree to within the margin" = "equivalent",
                "platforms differ" = "differ")
HEADER <- c("A4 screening" = "A4 and LEARN screening",
            "A4 repeated draws" = "A4 and LEARN repeated draws")

codes <- names(MEASURE_LABEL)[order(-nchar(names(MEASURE_LABEL)))]
code_rx <- paste0("(?<![A-Za-z0-9_])(", paste(codes, collapse = "|"), ")(?![A-Za-z0-9_])")
platform_of <- function(code) PLATFORM[[sub("_.*$", "", code)]]

relabel_cell <- function(x) {
  out <- x
  hit <- !is.na(x) & x %in% names(WHOLE_CELL)
  out[hit] <- WHOLE_CELL[x[hit]]
  # "A minus B" of two measure codes becomes the two platform names; the
  # analyte is in its own column, or in the table's title for the p-tau217
  # classification tables.
  pr <- !is.na(out) & str_detect(out, "^[A-Za-z0-9_]+ minus [A-Za-z0-9_]+$")
  for (i in which(pr)) {
    parts <- strsplit(out[i], " minus ", fixed = TRUE)[[1]]
    if (all(parts %in% names(MEASURE_LABEL)) && all(sub("_.*$", "", parts) %in% names(PLATFORM)))
      out[i] <- paste(platform_of(parts[1]), "minus", platform_of(parts[2]))
  }
  tok <- !is.na(out) & grepl(code_rx, out, perl = TRUE)
  for (i in which(tok)) {
    m <- gregexpr(code_rx, out[i], perl = TRUE)
    regmatches(out[i], m) <- list(unname(MEASURE_LABEL[regmatches(out[i], m)[[1]]]))
  }
  out
}

files <- list.files(MS_DIR, pattern = "^(Table[0-9]+|S[0-9]+)\\.csv$", full.names = TRUE)
changed <- 0L
for (f in files) {
  # Every cell is read and written as text exactly as it is, so that no number
  # is reformatted and no blank cell becomes NA.
  x <- readr::read_csv(f, col_types = readr::cols(.default = readr::col_character()),
                       na = character(), show_col_types = FALSE, progress = FALSE)
  y <- x %>% mutate(across(everything(), relabel_cell))
  nm <- names(y); hit <- nm %in% names(HEADER); nm[hit] <- HEADER[nm[hit]]; names(y) <- nm
  n <- sum(as.matrix(y) != as.matrix(x)) + sum(names(y) != names(x))
  if (n > 0) { readr::write_csv(y, f, na = ""); changed <- changed + 1L }
  say("  %-11s %4d labels changed", basename(f), n)
}
left <- unlist(lapply(files, function(f) {
  x <- readr::read_csv(f, col_types = readr::cols(.default = readr::col_character()),
                       na = character(), show_col_types = FALSE, progress = FALSE)
  v <- unlist(x, use.names = FALSE)
  v[grepl(code_rx, v, perl = TRUE)]
}))
if (length(left)) {
  say("  WARNING: variable names still present:")
  print(unique(left))
} else say("  no variable names remain in the formatted tables")

say("Written to %s, %s and %s", TAB_DIR, CHECK_DIR, MS_DIR)
say("Steps 07, 09 and 10 write the tables with variable names, so rerun this step after them.")
})
