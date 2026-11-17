# =============================================================================
# 08_intervals.R
# Stage 8. Confidence intervals for estimates of stages 5 to 7 that were first
# reported without one. Each section first reproduces the estimate, then
# resamples participants (1,000 resamples, seed 1).
#   8.1  the cleared ratio against PET (eMethods 8.9) and the eGFR rate
#        correlation with trial arm terms (eTable 19)
#   8.2  the corrected within-participant estimates (eTable 18), the eGFR rate
#        correlation without hypertension and diabetes (eTable 29) and with a
#        visit effect (eTable 25)
# Writes paper1_out/tables/review/review_8_*.csv, review_9_*.csv.
# =============================================================================

# #############################################################################
# Section 8.1. Intervals for the cleared ratio against PET and for the arm-term model
# Developed as 15_review_intervals.R version 1.0.0, whole
# #############################################################################
local({
SCRIPT_VERSION <- "1.0.0"
cat("15_review_intervals.R version", SCRIPT_VERSION, "\n")
# =============================================================================
# 15_review_intervals.R
# Intervals that the revised supplement still lacks: two agreement
# percentages (at two Centiloid thresholds) and one slope correlation.
#
#   A  Agreement of the FDA-cleared ratio categories with amyloid PET, with the
#      samples called positive or negative as the denominator: of the samples
#      called positive, the percentage PET positive (85.23% at a Centiloid of
#      24), and of those called negative, the percentage PET negative (93.81%).
#      Step 13 wrote the counts to checks/cleared_ratio_calibration.csv
#      without intervals. Here the counts are reproduced and the percentages
#      get intervals from 1,000 resamples of participants (seed 1). Both
#      Centiloid thresholds are given.
#   B  The eGFR slope correlation from the joint growth model with trial arm
#      on the level and the rate of each series (0.23 in
#      review_3b_arm_by_time.csv, against 0.22 without arm terms). Interval
#      from PAPER1_ARM_R_BOOT resamples of participants (default 1000), seed 1,
#      so the resamples are the same as those of the main model in T5c and
#      review_3; only converged fits are used, and the number usable is given.
#
# The joint model functions are copied unchanged from section 3 of
# 14_review_checks.R 1.1.0. Before anything is reported the script reproduces
# the main correlation (T5b), the correlation with arm terms (review_3b) and
# the calibration counts; a mismatch stops it.
#
# Reads    paper1_out/adni_visit.rds, a4_long.rds, a4_lab.rds
#          paper1_out/tables/T5b_joint_growth.csv,
#            review/review_3b_arm_by_time.csv
#          paper1_out/checks/cleared_ratio_calibration.csv
# Writes   paper1_out/tables/review/review_8_missing_intervals.csv
#
# Run time: part A takes under a minute; part B about two hours at 1,000
# draws. Set Sys.setenv(PAPER1_ARM_R_BOOT = "0") before running to skip B.
# =============================================================================

if (!exists("PAPER1_VERSION")) source("00_setup.R")
if (!requireNamespace("nlme", quietly = TRUE)) stop("The nlme package is needed for step 15.")

REV_DIR <- file.path(TAB_DIR, "review")
B_ARM_R <- suppressWarnings(as.integer(Sys.getenv("PAPER1_ARM_R_BOOT", "1000")))
if (is.na(B_ARM_R) || B_ARM_R < 0) B_ARM_R <- 1000
RESAMPLES <- 1000

adni <- readRDS(file.path(OUT_DIR, "adni_visit.rds"))
L    <- readRDS(file.path(OUT_DIR, "a4_long.rds"))
lab  <- readRDS(file.path(OUT_DIR, "a4_lab.rds"))
T5b  <- read_csv(file.path(TAB_DIR, "T5b_joint_growth.csv"), show_col_types = FALSE)
R3b  <- read_csv(file.path(REV_DIR, "review_3b_arm_by_time.csv"), show_col_types = FALSE)
calib_ref <- read_csv(file.path(CHECK_DIR, "cleared_ratio_calibration.csv"), show_col_types = FALSE)

# Participant resampling, as in step 14: the interval is the 2.5th to 97.5th
# percentile over B resamples of participants (RID).
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
# A. The cleared categories against amyloid PET, with intervals
# =============================================================================
say("A. Agreement of the cleared categories with amyloid PET, by category called")
FUJI_P217 <- "Fuji_plasma_ptau217"; FUJI_AB42 <- "Fuji_plasma_Ab42"
CUT_NEGATIVE <- 0.00370; CUT_POSITIVE <- 0.00738
cleared_category <- function(r) {
  r5 <- round(r, 5)
  dplyr::case_when(is.na(r5) ~ NA_character_, r5 <= CUT_NEGATIVE ~ "negative",
                   r5 >= CUT_POSITIVE ~ "positive", TRUE ~ "indeterminate")
}
cr <- adni %>% filter(!is.na(.data[[FUJI_P217]]), !is.na(.data[[FUJI_AB42]]), !is.na(centiloid)) %>%
  mutate(r = exp(.data[[FUJI_P217]] - .data[[FUJI_AB42]]), cat0 = cleared_category(r))
RA <- purrr::map_dfr(c(24, 20), function(cut) {
  d <- cr %>% mutate(pet = centiloid >= cut)
  # the counts must equal those step 13 wrote
  ref <- calib_ref %>% filter(centiloid_threshold == cut)
  here <- d %>% count(cat0, pet)
  if (nrow(ref) != 3 || sum(here$n) != sum(ref$`PET negative`, ref$`PET positive`))
    stop("The calibration totals differ at Centiloid ", cut, ". Nothing is reported.")
  for (k in c("positive", "negative", "indeterminate")) {
    want_pos <- ref$`PET positive`[ref$cat0 == k]; want_neg <- ref$`PET negative`[ref$cat0 == k]
    got_pos <- sum(here$n[here$cat0 == k & here$pet]); got_neg <- sum(here$n[here$cat0 == k & !here$pet])
    if (length(want_pos) != 1 || got_pos != want_pos || got_neg != want_neg)
      stop("The cleared categories do not reproduce cleared_ratio_calibration.csv at Centiloid ", cut,
           " (", k, "). Nothing is reported.")
  }
  fun <- function(x) c(pet_positive_of_called_positive = 100 * mean(x$pet[x$cat0 == "positive"]),
                       pet_negative_of_called_negative = 100 * mean(!x$pet[x$cat0 == "negative"]))
  boot_pct(d, fun) %>%
    mutate(part = "A", centiloid_cut = cut, samples = nrow(d), people = n_distinct(d$RID),
           n_called = c(sum(d$cat0 == "positive"), sum(d$cat0 == "negative")),
           draws = RESAMPLES, .before = 1)
})
say("  check: counts reproduce step 13 at Centiloid 20 and 24")
print(as.data.frame(RA %>% transmute(centiloid_cut, stat, samples, n_called,
                                     value = sprintf("%.2f (%.2f to %.2f)", est, lo, hi), usable)),
      row.names = FALSE)

# =============================================================================
# B. The slope correlation with trial arm on the level and rate of each series
# =============================================================================
# Copied unchanged from section 3 of 14_review_checks.R 1.1.0.
say("B. The eGFR slope correlation with trial arm in the model")
TIME_CENTRE <- 2.5
MIN_VALUES <- 2
MIN_SPAN_YEARS <- 1
JOINT_COVS <- c("age_z", "female", "e4", "cl_z")
NONCONVERGENCE <- "convergence error code|without convergence"
eligible_people <- function(v) {
  f <- lab %>% filter(!is.na(.data[[v]])) %>% group_by(BID) %>%
    summarise(n = n(), span = max(lab_time) - min(lab_time), .groups = "drop") %>%
    filter(n >= MIN_VALUES, span >= MIN_SPAN_YEARS)
  b <- L %>% filter(!is.na(pTau217)) %>% group_by(BID) %>%
    summarise(n = n(), .groups = "drop") %>% filter(n >= MIN_VALUES)
  intersect(f$BID, b$BID)
}
person_covs <- L %>% distinct(BID, age0, female, e4, centiloid, arm) %>%
  group_by(BID) %>% slice(1) %>% ungroup()
joint_fixed <- function(extra = character()) {
  base <- c("ix", "ip", "ix:t", "ip:t")
  adj  <- as.vector(outer(base, c(JOINT_COVS, extra), paste, sep = ":"))
  stats::as.formula(paste("y ~ 0 +", paste(c(base, adj), collapse = " + ")))
}
joint_data <- function(v, people = NULL) {
  keep <- eligible_people(v)
  if (!is.null(people)) keep <- intersect(keep, people)
  cv <- person_covs %>% filter(BID %in% keep) %>%
    filter(!is.na(age0), !is.na(female), !is.na(e4), !is.na(centiloid)) %>%
    mutate(age_z = as.numeric(scale(age0)), cl_z = as.numeric(scale(centiloid)),
           sola = as.integer(arm == "Solanezumab"), learn = as.integer(arm == "LEARN")) %>%
    select(BID, all_of(JOINT_COVS), sola, learn, arm)
  keep <- intersect(keep, cv$BID)
  bio <- L %>% filter(!is.na(pTau217), BID %in% keep) %>%
    transmute(BID, t = plasma_time, y = pTau217, o = "bio")
  fac <- lab %>% filter(!is.na(.data[[v]]), BID %in% keep) %>%
    transmute(BID, t = lab_time, y = FACTOR_SIGN[[v]] * .data[[v]], o = "fac")
  raw <- bind_rows(bio, fac)
  sds <- raw %>% group_by(o) %>% summarise(s = stats::sd(y), .groups = "drop")
  out <- raw %>%
    group_by(o) %>% mutate(y = as.numeric(scale(y))) %>% ungroup() %>%
    left_join(cv, by = "BID") %>%
    mutate(t = t - TIME_CENTRE, ip = as.integer(o == "bio"), ix = as.integer(o == "fac")) %>%
    arrange(BID, o, t)
  attr(out, "sd_bio") <- sds$s[sds$o == "bio"]
  attr(out, "sd_fac") <- sds$s[sds$o == "fac"]
  out
}
fit_joint <- function(d, extra = character(), optimiser = "optim") {
  ctrl <- nlme::lmeControl(opt = optimiser, maxIter = 200, msMaxIter = 200,
                           returnObject = TRUE, apVar = FALSE)
  warned <- character()
  m <- withCallingHandlers(tryCatch(
    nlme::lme(joint_fixed(extra),
              random = list(BID = nlme::pdSymm(~ 0 + ix + ip + ix:t + ip:t)),
              weights = nlme::varIdent(form = ~ 1 | o),
              data = d, control = ctrl, method = "REML"),
    error = function(e) NULL),
    warning = function(w) { warned <<- c(warned, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (!is.null(m)) attr(m, "converged") <- !any(grepl(NONCONVERGENCE, warned))
  m
}
converged <- function(m) !is.null(m) && isTRUE(attr(m, "converged"))
slope_correlation <- function(m) {
  if (is.null(m)) return(NA_real_)
  vc <- as.matrix(nlme::getVarCov(m)); nm <- colnames(vc)
  i <- grep("ix:t|t:ix", nm)[1]; j <- grep("ip:t|t:ip", nm)[1]
  if (is.na(i) || is.na(j)) return(NA_real_)
  vc[i, j] / sqrt(vc[i, i] * vc[j, j])
}
resample_people <- function(d, people) {
  ids <- sample(people, length(people), replace = TRUE)
  tibble(BID = ids, .new = seq_along(ids)) %>%
    inner_join(d, by = "BID", relationship = "many-to-many") %>%
    mutate(BID = .new) %>% select(-.new)
}

d3 <- joint_data("egfr")
arm_tab <- d3 %>% distinct(BID, arm) %>% count(arm)
ARM_TERMS <- c("sola", if (any(arm_tab$arm == "LEARN")) "learn")
r_main <- slope_correlation(fit_joint(d3))
ref_main <- T5b$slope_correlation[T5b$var == "egfr" & T5b$group == "everyone"]
m_arm <- fit_joint(d3, extra = ARM_TERMS)
r_arm <- slope_correlation(m_arm)
ref_arm <- unique(R3b$r_with_arm)
say("  check: correlation %.4f here, %.4f in T5b; with arm terms %.4f here, %.4f in review_3b",
    r_main, ref_main, r_arm, ref_arm)
if (length(ref_main) != 1 || length(ref_arm) != 1 || is.na(r_main) || is.na(r_arm) ||
    abs(r_main - ref_main) > 1e-4 || abs(r_arm - ref_arm) > 1e-4 || !converged(m_arm))
  stop("The joint models do not reproduce steps 06 and 14. Nothing further is reported.")

RB <- tibble(part = "B", centiloid_cut = NA_real_, samples = nrow(d3),
             people = n_distinct(d3$BID), n_called = NA_integer_, draws = B_ARM_R,
             stat = "slope correlation with arm terms", est = r_arm,
             lo = NA_real_, hi = NA_real_, usable = NA_integer_)
if (B_ARM_R > 0) {
  people <- unique(d3$BID)
  set.seed(1)
  t0 <- Sys.time()
  res <- vapply(seq_len(B_ARM_R), function(b) {
    if (b %% 50 == 0) say("  arm-term model: %d of %d, %.0f minutes so far", b, B_ARM_R,
                          as.numeric(difftime(Sys.time(), t0, units = "mins")))
    dd <- resample_people(d3, people)
    ma <- tryCatch(fit_joint(dd, extra = ARM_TERMS), error = function(e) NULL)
    if (converged(ma)) slope_correlation(ma) else NA_real_
  }, numeric(1))
  RB$lo <- unname(stats::quantile(res, 0.025, na.rm = TRUE))
  RB$hi <- unname(stats::quantile(res, 0.975, na.rm = TRUE))
  RB$usable <- sum(!is.na(res))
}
out <- bind_rows(RA, RB)
write_csv(out, file.path(REV_DIR, "review_8_missing_intervals.csv"))
say("  written review/review_8_missing_intervals.csv, %d rows", nrow(out))
print(as.data.frame(out %>% transmute(part, centiloid_cut, stat,
                                      value = sprintf("%.2f (%.2f to %.2f)", est, lo, hi), usable)),
      row.names = FALSE)
say("Step 15 finished.")
})


# #############################################################################
# Section 8.2. Intervals for eTables 18, 25 and 29
# Developed as 16_missing_intervals.R version 1.1.0
# The tasks no longer reported (part B, and creatine kinase and haemoglobin in C and D)
# are removed; the two eGFR tasks are unchanged.
# #############################################################################
local({
SCRIPT_VERSION <- "1.1.0"
cat("16_missing_intervals.R version", SCRIPT_VERSION, "\n")
# =============================================================================
# 16_missing_intervals.R
# Confidence intervals for the estimates in the supplement that were reported
# without one. Nothing here changes an estimate: every model is first refitted
# on the full sample and must reproduce the value already reported, or the
# step stops.
#
# Changes in 1.1.0 (decision of 16 Nov 2026): only the two eGFR tasks are run,
#   C_egfr_without and D_egfr_visit. The creatine kinase and haemoglobin rate
#   correlations are removed from eTables 25, 27 and 29, because no statement
#   in the paper rests on them and their main estimates already have intervals
#   (Table 5, eTable 19). On reorganisation the code for the other tasks (part
#   B and the creatine kinase and haemoglobin tasks of C and D) was removed.
#
#   A  eTable 18: the within-participant estimate corrected for the
#      reliability of change. The interval is the within-participant interval
#      divided by the same reliability on the log scale, as eMethods 8.4 does
#      for the corrected level estimates, so it does not carry the uncertainty
#      in the reliability. From T5_hybrid_within_between.csv; no model is fitted.
#   C  eTable 29: the eGFR slope correlation of step 12 on the people with
#      hypertension and diabetes known, without the two conditions (the
#      interval with them comes from step 12). Reproduces T11_comorbidity.csv.
#   D  eTable 25: the eGFR slope correlation of step 10 with a disturbance
#      shared by both series at one visit, and the within-visit residual
#      correlation. Each resample uses the specification that fitted the full
#      sample (correlated residuals within visit). Reproduces
#      T10_occasion_check.csv.
#
# Every interval is the 2.5th to 97.5th percentile over 1,000 resamples of
# people (PAPER1_MI_BOOT), each drawn person a new identifier, seed 1 for every
# task, as in steps 06, 12 and 14; only converged fits are used, and the number
# usable is reported. The tasks run side by side on PAPER1_MI_CORES cores
# (default: all but one). Each finished task is saved at once in
# tables/review/missing_intervals_parts/, and a rerun skips the tasks already
# saved, so an interrupted run can be resumed by sourcing the step again.
#
# Reads    paper1_out/a4_long.rds, a4_lab.rds; A4 Raw Data/vitals.csv (read only)
#          paper1_out/tables/T5_hybrid_within_between.csv,
#            T10_occasion_check.csv, T11_comorbidity.csv
# Writes   paper1_out/tables/review/review_9_missing_intervals.csv
#          paper1_out/tables/review/missing_intervals_parts/*.rds and *.log
#
# Run time: the two eGFR tasks run side by side. The eGFR task without
# hypertension and diabetes takes about two hours; the visit-effect task takes
# several hours (15 at the speed of the 20-resample trial, less with
# VECLIB_MAXIMUM_THREADS=1 in ~/.Renviron), so run it overnight with the Mac on
# mains power; the step keeps the Mac from sleeping while it runs. Before a
# rerun on changed data, move missing_intervals_parts to _to_delete, because
# saved tasks are reused.
# Set Sys.setenv(PAPER1_MI_BOOT = "20") first for a quick trial of the code.
# =============================================================================

if (!exists("PAPER1_VERSION")) source("00_setup.R")
for (pkg in c("nlme", "parallel"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("The ", pkg, " package is needed for step 16.")

REV_DIR  <- file.path(TAB_DIR, "review")
PART_DIR <- file.path(REV_DIR, "missing_intervals_parts")
dir.create(PART_DIR, showWarnings = FALSE, recursive = TRUE)
env_int <- function(name, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (is.na(v) || v < 1) default else v
}
B     <- env_int("PAPER1_MI_BOOT", 1000)
DC    <- suppressWarnings(parallel::detectCores())
if (is.na(DC)) DC <- 2L
CORES <- env_int("PAPER1_MI_CORES", max(1L, DC - 1L))
if (CORES < 1) CORES <- 1
say("Resamples per task: %d; cores: %d", B, CORES)
say("BLAS: %s", sessionInfo()$BLAS)
say("  If the BLAS line names vecLib or Accelerate and the step stalls, rerun with Sys.setenv(PAPER1_MI_CORES = \"1\").")
# Keep the Mac awake while this R session runs (macOS only; nothing is changed
# when R exits).
if (Sys.info()[["sysname"]] == "Darwin")
  try(system2("caffeinate", c("-i", "-w", Sys.getpid()), wait = FALSE), silent = TRUE)

L   <- readRDS(file.path(OUT_DIR, "a4_long.rds"))
lab <- readRDS(file.path(OUT_DIR, "a4_lab.rds"))
rd  <- function(f) read_csv(file.path(TAB_DIR, f), show_col_types = FALSE)
T5  <- rd("T5_hybrid_within_between.csv")
T10 <- rd("T10_occasion_check.csv")
T11 <- rd("T11_comorbidity.csv")

# =============================================================================
# A. eTable 18, corrected within-participant estimates
# =============================================================================
say("A. Corrected within-participant estimates (eTable 18)")
RA <- T5 %>% filter(!is.na(within_corrected_pct)) %>%
  transmute(part = "A", table = "eTable 18", var, analysis = "within participant, corrected",
            estimate = within_corrected_pct,
            lo = pct(log1p(within_lo / 100) / rel_within),
            hi = pct(log1p(within_hi / 100) / rel_within),
            draws = NA_integer_, usable = NA_integer_, bootstrap_median = NA_real_)
chk <- abs(RA$estimate - pct(T5$within_b[match(RA$var, T5$var)] / T5$rel_within[match(RA$var, T5$var)]))
if (any(chk > 1e-6)) stop("The corrected estimates do not reproduce T5. Nothing is reported.")
print(as.data.frame(RA %>% transmute(var, value = sprintf("%.2f (%.2f to %.2f)", estimate, lo, hi))),
      row.names = FALSE)

# =============================================================================
# The joint growth model, as in steps 06, 10 and 12
# =============================================================================
TIME_CENTRE <- 2.5
MIN_VALUES <- 2
MIN_SPAN_YEARS <- 1
OCC_TOLERANCE <- 0.10
JOINT_COVS <- c("age_z", "female", "e4", "cl_z")
NONCONVERGENCE <- "convergence error code|without convergence"

# Hypertension and diabetes for the people of the joint model, as step 12
# builds them: screening blood pressure and the HbA1c at the earliest
# laboratory visit that has one.
a4_vt <- read_source(file.path(A4_DIR, "Raw Data", "vitals.csv"), "A4 vitals") %>%
  filter(str_trim(VISCODE) == "001") %>%
  transmute(BID,
            sbp = suppressWarnings(as.numeric(VSBPSYS)),
            dbp = suppressWarnings(as.numeric(VSBPDIA))) %>%
  mutate(sbp = if_else(sbp >= 60 & sbp <= 260, sbp, NA_real_),
         dbp = if_else(dbp >= 30 & dbp <= 150, dbp, NA_real_)) %>%
  group_by(BID) %>%
  summarise(sbp = mean(sbp, na.rm = TRUE), dbp = mean(dbp, na.rm = TRUE), .groups = "drop") %>%
  mutate(htn = case_when(is.na(sbp) & is.na(dbp) ~ NA_integer_,
                         coalesce(sbp >= 140, FALSE) | coalesce(dbp >= 90, FALSE) ~ 1L,
                         TRUE ~ 0L)) %>%
  select(BID, htn)
a4_person_cm <- lab %>%
  filter(!is.na(hba1c)) %>%
  group_by(BID) %>% slice_min(day, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(BID, dm = as.integer(hba1c >= 6.5)) %>%
  full_join(a4_vt, by = "BID")

person_covs <- L %>% distinct(BID, age0, female, e4, centiloid) %>%
  group_by(BID) %>% slice(1) %>% ungroup() %>%
  left_join(a4_person_cm, by = "BID")

# Eligibility and data, step 06 (span rules and the out-of-range check) and
# step 12 (people with both conditions known).
bio_draws <- function(quantified_only = FALSE) if (quantified_only) L %>% filter(quantified) else L
eligible_people <- function(v, lab_span = MIN_SPAN_YEARS, bio_span = 0, quantified_only = FALSE) {
  f <- lab %>% filter(!is.na(.data[[v]])) %>% group_by(BID) %>%
    summarise(n = n(), span = max(lab_time) - min(lab_time), .groups = "drop") %>%
    filter(n >= MIN_VALUES, span >= lab_span)
  b <- bio_draws(quantified_only) %>% filter(!is.na(pTau217)) %>% group_by(BID) %>%
    summarise(n = n(), span = max(plasma_time) - min(plasma_time), .groups = "drop") %>%
    filter(n >= MIN_VALUES, span >= bio_span)
  intersect(f$BID, b$BID)
}
joint_data <- function(v, lab_span = MIN_SPAN_YEARS, bio_span = 0, quantified_only = FALSE,
                       need_cm = FALSE) {
  keep <- eligible_people(v, lab_span, bio_span, quantified_only)
  cv <- person_covs %>% filter(BID %in% keep) %>%
    filter(!is.na(age0), !is.na(female), !is.na(e4), !is.na(centiloid))
  if (need_cm) cv <- cv %>% filter(!is.na(htn), !is.na(dm))
  cv <- cv %>% mutate(age_z = as.numeric(scale(age0)), cl_z = as.numeric(scale(centiloid))) %>%
    select(BID, all_of(JOINT_COVS), any_of(if (need_cm) c("htn", "dm")))
  keep <- intersect(keep, cv$BID)
  bio <- bio_draws(quantified_only) %>% filter(!is.na(pTau217), BID %in% keep) %>%
    transmute(BID, t = plasma_time, y = pTau217, o = "bio")
  fac <- lab %>% filter(!is.na(.data[[v]]), BID %in% keep) %>%
    transmute(BID, t = lab_time, y = FACTOR_SIGN[[v]] * .data[[v]], o = "fac")
  bind_rows(bio, fac) %>%
    group_by(o) %>% mutate(y = as.numeric(scale(y))) %>% ungroup() %>%
    left_join(cv, by = "BID") %>%
    mutate(t = t - TIME_CENTRE, ip = as.integer(o == "bio"), ix = as.integer(o == "fac")) %>%
    arrange(BID, o, t)
}
joint_fixed <- function(extra = character(), quadratic = FALSE) {
  base <- c("ix", "ip", "ix:t", "ip:t", if (quadratic) c("ix:I(t^2)", "ip:I(t^2)"))
  adj  <- as.vector(outer(c("ix", "ip", "ix:t", "ip:t"), c(JOINT_COVS, extra), paste, sep = ":"))
  stats::as.formula(paste("y ~ 0 +", paste(c(base, adj), collapse = " + ")))
}
lme_checked <- function(...) {
  warned <- character()
  m <- withCallingHandlers(tryCatch(nlme::lme(...), error = function(e) NULL),
    warning = function(w) { warned <<- c(warned, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (!is.null(m)) attr(m, "converged") <- !any(grepl(NONCONVERGENCE, warned))
  m
}
converged <- function(m) !is.null(m) && isTRUE(attr(m, "converged"))
ctrl <- function(opt = "optim") nlme::lmeControl(opt = opt, maxIter = 200, msMaxIter = 200,
                                                 returnObject = TRUE, apVar = FALSE)
# visit = "none": the model of steps 06 and 12; "corr": residuals correlated
# within a visit; "intercept": a random intercept at the visit (step 10).
fit_model <- function(d, extra = character(), quadratic = FALSE, visit = "none", opt = "optim") {
  fx <- joint_fixed(extra, quadratic)
  if (visit == "none")
    return(lme_checked(fx, random = list(BID = nlme::pdSymm(~ 0 + ix + ip + ix:t + ip:t)),
                       weights = nlme::varIdent(form = ~ 1 | o),
                       data = d, control = ctrl(opt), method = "REML"))
  d$occ <- factor(paste(d$BID, round(d$t / OCC_TOLERANCE)))
  d <- d %>% arrange(BID, occ)
  if (visit == "corr")
    return(lme_checked(fx, random = list(BID = nlme::pdSymm(~ 0 + ix + ip + ix:t + ip:t)),
                       correlation = nlme::corCompSymm(form = ~ 1 | BID/occ),
                       weights = nlme::varIdent(form = ~ 1 | o),
                       data = d, control = ctrl(opt), method = "REML"))
  lme_checked(fx, random = list(BID = nlme::pdSymm(~ 0 + ix + ip + ix:t + ip:t),
                                occ = nlme::pdIdent(~ 1)),
              weights = nlme::varIdent(form = ~ 1 | o),
              data = d, control = ctrl(opt), method = "REML")
}
# The covariance of the person-level random effects; getVarCov cannot return it
# for a model with two levels of nesting, so it is taken from the fitted
# structure in that case, as step 10 does.
person_cov <- function(m) {
  rs <- m$modelStruct$reStruct
  pd <- if ("BID" %in% names(rs)) rs[["BID"]] else rs[[length(rs)]]
  tryCatch(as.matrix(nlme::pdMatrix(pd)) * m$sigma^2, error = function(e) NULL)
}
slope_correlation <- function(m) {
  if (is.null(m)) return(NA_real_)
  vc <- tryCatch(as.matrix(nlme::getVarCov(m)), error = function(e) NULL)
  if (is.null(vc) || !is.matrix(vc) || is.null(colnames(vc))) vc <- person_cov(m)
  if (is.null(vc) || is.null(colnames(vc))) return(NA_real_)
  nm <- colnames(vc)
  i <- grep("ix:t|t:ix", nm)[1]; j <- grep("ip:t|t:ip", nm)[1]
  if (is.na(i) || is.na(j)) return(NA_real_)
  vc[i, j] / sqrt(vc[i, i] * vc[j, j])
}
within_visit_rho <- function(m) {
  if (is.null(m) || is.null(m$modelStruct$corStruct)) return(NA_real_)
  tryCatch(unname(coef(m$modelStruct$corStruct, unconstrained = FALSE))[1], error = function(e) NA_real_)
}
stats_of <- function(m, rho) {
  if (!converged(m)) return(c(r = NA_real_, rho = NA_real_))
  c(r = slope_correlation(m), rho = if (rho) within_visit_rho(m) else NA_real_)
}

# =============================================================================
# The tasks
# =============================================================================
tasks <- list()
jm <- T11 %>% filter(analysis == "joint growth model")
for (v in "egfr") {
  tasks[[length(tasks) + 1]] <- list(
    id = paste("C", v, "without", sep = "_"), part = "C", table = "eTable 29", var = v,
    analysis = "people with both conditions known, without them",
    data = list(need_cm = TRUE), extra = character(), quadratic = FALSE, visit = "none", rho = FALSE,
    ref = jm$r0[jm$var == v], ref_rho = NA_real_)
}
for (v in "egfr") {
  how <- T10$fitted_by[T10$var == v]
  vis <- if (identical(how, "correlated residuals within visit")) "corr" else "intercept"
  tasks[[length(tasks) + 1]] <- list(
    id = paste("D", v, "visit", sep = "_"), part = "D", table = "eTable 25", var = v,
    analysis = paste("with a visit effect:", how), data = list(need_cm = FALSE),
    extra = character(), quadratic = FALSE, visit = vis, rho = vis == "corr",
    ref = T10$correlation_with_visit_effect[T10$var == v],
    ref_rho = T10$within_visit_residual_correlation[T10$var == v])
}
say("%d tasks", length(tasks))

# Full-sample fits first, in this session: each must reproduce its table.
say("Checking that every full-sample fit reproduces the reported estimate")
for (k in seq_along(tasks)) {
  tk <- tasks[[k]]
  d <- do.call(joint_data, c(list(v = tk$var), tk$data))
  s <- stats_of(fit_model(d, tk$extra, tk$quadratic, tk$visit), tk$rho)
  ok <- length(tk$ref) == 1 && !is.na(s[["r"]]) && abs(s[["r"]] - tk$ref) < 1e-4 &&
    (!tk$rho || (!is.na(s[["rho"]]) && abs(s[["rho"]] - tk$ref_rho) < 1e-4))
  say("  %-26s %.4f here, %.4f reported%s", tk$id, s[["r"]], tk$ref,
      if (tk$rho) sprintf("; within-visit %.4f here, %.4f reported", s[["rho"]], tk$ref_rho) else "")
  tasks[[k]]$opt <- "optim"
  # Step 10 refits with nlminb when optim does not give an answer; do the same.
  if (!ok && tk$part == "D") {
    s2 <- stats_of(fit_model(d, tk$extra, tk$quadratic, tk$visit, opt = "nlminb"), tk$rho)
    ok2 <- !is.na(s2[["r"]]) && abs(s2[["r"]] - tk$ref) < 1e-4 &&
      (!tk$rho || (!is.na(s2[["rho"]]) && abs(s2[["rho"]] - tk$ref_rho) < 1e-4))
    if (ok2) { s <- s2; ok <- TRUE; tasks[[k]]$opt <- "nlminb"; say("    reproduced with nlminb, which the resamples will use") }
  }
  if (!ok) stop("Task ", tk$id, " does not reproduce the reported estimate. Nothing further is run.")
  tasks[[k]]$estimate <- s
  tasks[[k]]$n_people <- n_distinct(d$BID)
}

# =============================================================================
# The resamples, one task per core
# =============================================================================
run_task <- function(tk) {
  f <- file.path(PART_DIR, paste0(tk$id, "_B", B, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  d <- do.call(joint_data, c(list(v = tk$var), tk$data))
  people <- unique(d$BID)
  set.seed(1)
  t0 <- Sys.time()
  res <- vapply(seq_len(B), function(b) {
    ids <- sample(people, length(people), replace = TRUE)
    dd <- tibble(BID = ids, .new = seq_along(ids)) %>%
      inner_join(d, by = "BID", relationship = "many-to-many") %>%
      mutate(BID = .new) %>% select(-.new)
    m <- tryCatch(fit_model(dd, tk$extra, tk$quadratic, tk$visit, opt = tk$opt), error = function(e) NULL)
    if (b %% 50 == 0) cat(sprintf("%s  %d of %d, %.0f minutes\n", format(Sys.time(), "%H:%M"), b, B,
                                  as.numeric(difftime(Sys.time(), t0, units = "mins"))),
                          file = file.path(PART_DIR, paste0(tk$id, ".log")), append = TRUE)
    stats_of(m, tk$rho)
  }, numeric(2))
  out <- list(id = tk$id, r = res["r", ], rho = res["rho", ],
              minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(out, f)
  out
}
todo <- vapply(tasks, function(tk) !file.exists(file.path(PART_DIR, paste0(tk$id, "_B", B, ".rds"))), logical(1))
say("Resampling: %d tasks to run, %d already saved. Progress is not printed while the", sum(todo), sum(!todo))
say("  tasks run; each writes its progress every 50 resamples to a .log file, and its")
say("  result to an .rds file when it finishes, in %s", PART_DIR)
t_all <- Sys.time()
# Anything not yet written to the run log is written now, so that the forked
# workers do not carry a copy of it.
if (exists("log_con")) try(flush(log_con), silent = TRUE)
# The visit-effect tasks are the slowest, so they start first.
ord <- order(vapply(tasks, function(tk) tk$part != "D", logical(1)))
res <- if (.Platform$OS.type == "windows" || CORES == 1) lapply(tasks[ord], run_task) else
  parallel::mclapply(tasks[ord], run_task, mc.cores = CORES, mc.preschedule = FALSE)
done <- vector("list", length(tasks)); done[ord] <- res
bad <- vapply(done, function(x) inherits(x, "try-error") || is.null(x$r), logical(1))
if (any(bad)) stop("Some tasks failed: ", paste(vapply(tasks[bad], `[[`, "", "id"), collapse = ", "),
                   ". Source the step again to rerun only those.")
say("Resampling finished in %.0f minutes", as.numeric(difftime(Sys.time(), t_all, units = "mins")))

pq <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE))
RB <- purrr::map2_dfr(tasks, done, function(tk, dn) {
  rows <- tibble(part = tk$part, table = tk$table, var = tk$var, analysis = tk$analysis,
                 estimate = tk$estimate[["r"]], lo = pq(dn$r, 0.025), hi = pq(dn$r, 0.975),
                 draws = B, usable = sum(!is.na(dn$r)), bootstrap_median = stats::median(dn$r, na.rm = TRUE),
                 n_people = tk$n_people, minutes = round(dn$minutes, 1))
  if (tk$rho) rows <- bind_rows(rows, tibble(
    part = tk$part, table = tk$table, var = tk$var, analysis = "within-visit residual correlation",
    estimate = tk$estimate[["rho"]], lo = pq(dn$rho, 0.025), hi = pq(dn$rho, 0.975),
    draws = B, usable = sum(!is.na(dn$rho)), bootstrap_median = stats::median(dn$rho, na.rm = TRUE),
    n_people = tk$n_people, minutes = round(dn$minutes, 1)))
  rows
})
out <- bind_rows(RA, RB)
write_csv(out, file.path(REV_DIR, "review_9_missing_intervals.csv"))
say("  written review/review_9_missing_intervals.csv, %d rows", nrow(out))
print(as.data.frame(out %>% transmute(table, var, analysis = substr(analysis, 1, 40),
                                      value = sprintf("%.2f (%.2f to %.2f)", estimate, lo, hi), usable)),
      row.names = FALSE)
say("Step 16 finished.")
})
