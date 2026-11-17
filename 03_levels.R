# =============================================================================
# 03_levels.R
# Stage 3. Research question 1: do routine physiological factors change plasma
# biomarker levels, in both cohorts?
#   3.1  every factor against every biomarker (Table 2, eTables 3 to 9, 28)
#   3.2  the four primary factors under six specifications (eTable 35)
# Writes paper1_out/tables/T1*.csv, R1_reliability.csv, review/review_7_*.csv.
# =============================================================================

# #############################################################################
# Section 3.1. Levels in both cohorts
# Developed as 03_rq1_levels.R version 1.3.0, whole
# #############################################################################
local({
# =============================================================================
# 03_rq1_levels.R
# Research question 1: do routine physiological factors change plasma biomarker
# levels, and does the same pattern appear in both cohorts?
#
# Both cohorts go through the same function, so the z-scoring, the signs, the
# covariates and the standard errors cannot drift apart between them. A4/LEARN
# has no diagnosis term because everyone enrolled was cognitively unimpaired.
#
# C2N p-tau217, and the %p-tau217 ratio built on it, are reported as half the
# detection limit in 424 of 1,179 ADNI samples. Treating that substituted value
# as a measurement biases the estimate, so every level model of those two
# measures is the left-censored regression of eMethods 5, the same one step 05
# uses: censored at the limit of 1.30 pg/mL for p-tau217, and at 100 x 1.30 /
# np-tau217 for the ratio. Every other measure is fitted by least squares.
# The column method in each table says which was used.
#
# False discovery rate: Benjamini-Hochberg within the families of eMethods 9,
# which for the levels are all biomarkers by all factors of one role within one
# cohort. The unadjusted p and the adjusted q are both kept in T1.
#
# Reads    paper1_out/adni_visit.rds, a4_screen.rds, a4_lab.rds
# Writes   paper1_out/tables/
#            T1_levels.csv             every factor against every measure, both cohorts
#            T1b_egfr_shape.csv        eGFR modelled flexibly, contrasts at 60 and 45
#            T1c_mutual_adjustment.csv the primary factors in one model
#            T1d_negative_control.csv  amyloid PET as the outcome
#            T1e_sex_stratified.csv    creatine kinase in men and women
#            T1f_lab_recency.csv       ADNI, laboratory values close to the draw
#            R1_reliability.csv        reliability of a single measurement, from
#                                      the A4 repeats; used here to correct for
#                                      measurement error and again in step 06
#
# Takes about a minute.
# =============================================================================
SCRIPT_VERSION <- "1.3.0"
cat("03_rq1_levels.R version", SCRIPT_VERSION, "\n")

if (!exists("PAPER1_VERSION")) source("00_setup.R")
if (!requireNamespace("splines", quietly = TRUE)) stop("The splines package is needed.")

adni <- readRDS(file.path(OUT_DIR, "adni_visit.rds"))
a4   <- readRDS(file.path(OUT_DIR, "a4_screen.rds"))
lab  <- readRDS(file.path(OUT_DIR, "a4_lab.rds"))

# The measures each cohort reports. ADNI has several platforms per analyte, so
# question 1 uses all of them and question 2 asks whether they agree.
ADNI_MEASURES <- c("C2N_plasma_ptau217", "Fuji_plasma_ptau217", "AlzPath_plasma_ptau217",
                   "Janssen_plasma_ptau217", "C2N_plasma_nptau217", "C2N_plasma_ptau217_ratio",
                   "Roche_plasma_ptau181", "QX_plasma_ptau181",
                   "Roche_plasma_NfL", "QX_plasma_NfL", "Roche_plasma_GFAP", "QX_plasma_GFAP",
                   "C2N_plasma_Abeta40", "Fuji_plasma_Ab40", "Roche_plasma_Ab40", "QX_plasma_Ab40",
                   "C2N_plasma_Abeta42", "Fuji_plasma_Ab42", "Roche_plasma_Ab42", "QX_plasma_Ab42",
                   "C2N_plasma_Abeta42_Abeta40", "Fuji_plasma_Ab42_Ab40",
                   "Roche_plasma_Ab42_Ab40", "QX_plasma_Ab42_Ab40")
A4_MEASURES   <- c("pTau181", "NfL", "GFAP", "Abeta40", "Abeta42", "Ab42_40")

short <- function(x) gsub("_plasma", "", x)

# ---- Which measures carry the C2N detection limit ----------------------------
# The numerator is censored at the limit itself; the ratio at a limit that moves
# with its denominator. Both take the censoring indicator from the laboratory's
# flag, which agrees exactly with the values (every flagged value is 0.65 and no
# unflagged value is at or below 1.30).
C2N_NUM   <- "C2N_plasma_ptau217"
C2N_DEN   <- "C2N_plasma_nptau217"
C2N_RATIO <- "C2N_plasma_ptau217_ratio"
C2N_CENSORED <- c(C2N_NUM, C2N_RATIO)

# The rows a model actually uses: complete cases on the outcome, the covariates
# and the cluster, and for the ratio also the denominator its limit needs. The
# censored fit is given only these rows, so that its counts of samples, of
# censored samples and of people are those of the fitted sample, whatever
# fit_tobit counts internally.
analysis_rows <- function(formula, data, cluster) {
  y    <- all.vars(formula)[1]
  vars <- unique(c(all.vars(formula), cluster,
                   if (y %in% C2N_CENSORED) "c2n_at_floor",
                   if (y == C2N_RATIO) C2N_DEN))
  d <- as.data.frame(data)
  d <- d[stats::complete.cases(d[, vars, drop = FALSE]), , drop = FALSE]
  rownames(d) <- NULL
  d
}

# One level model. Least squares with the CR1 sandwich for every measure, except
# the two C2N measures with a detection limit, which get the cluster-robust
# censored regression. Returns the columns of fit_cl plus method and n_censored.
fit_level <- function(formula, data, term, cluster) {
  y <- all.vars(formula)[1]
  if (!y %in% C2N_CENSORED)
    return(fit_cl(formula, data, term, cluster = cluster) %>%
             mutate(method = "least squares", n_censored = NA_integer_))
  d <- analysis_rows(formula, data, cluster)
  empty <- tibble(n_obs = nrow(d), n_people = dplyr::n_distinct(d[[cluster]]),
                  b = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_,
                  p = NA_real_, level = 0.95, method = "censored regression",
                  n_censored = sum(d$c2n_at_floor))
  if (nrow(d) < 10 || sum(!d$c2n_at_floor) < 10) return(empty)
  limit <- if (y == C2N_RATIO) 100 * C2N_LOD / exp(d[[C2N_DEN]]) else C2N_LOD
  r <- tryCatch(fit_tobit(formula, d, term, censor_at = limit, cluster = cluster,
                          censored = d$c2n_at_floor),
                error = function(e) NULL)
  if (is.null(r)) return(empty)
  empty %>% mutate(b = r$b, se = r$se, lo = r$lo, hi = r$hi, p = r$p)
}

# ---- The level model, used for both cohorts ---------------------------------
# One model per measure and factor: the log biomarker on the signed z-scored
# factor, with the cohort's covariates and standard errors clustered by person.
# vars limits the factors, for the laboratory recency analysis below.
levels_for <- function(data, measures, covs, id, cohort, vars = NULL) {
  factors <- FACTORS %>% filter(paste0("z_", var) %in% names(data))
  if (!is.null(vars)) factors <- factors %>% filter(var %in% vars)
  measures <- intersect(measures, names(data))
  purrr::map_dfr(measures, function(ms) purrr::map_dfr(seq_len(nrow(factors)), function(j) {
    zv <- paste0("z_", factors$var[j])
    fit_level(as.formula(sprintf("`%s` ~ %s + %s", ms, zv, covs)), data, zv, cluster = id) %>%
      mutate(cohort = cohort, measure = short(ms), measure_col = ms, factor = factors$label[j],
             var = factors$var[j], role = factors$role[j], .before = 1)
  }))
}

T1 <- bind_rows(
  levels_for(adni, ADNI_MEASURES, COVS_ADNI, "RID", "ADNI"),
  levels_for(a4,   A4_MEASURES,   COVS_A4,   "BID", "A4 screening")
) %>%
  mutate(pct_diff = pct(b), pct_lo = pct(lo), pct_hi = pct(hi)) %>%
  # eMethods 9: one family per cohort and factor role, across every biomarker.
  group_by(cohort, role) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  ungroup()

# ---- Reliability of a single measurement of each factor ---------------------
# From the repeated laboratory values in A4/LEARN. A factor measured once is a
# noisy measure of the underlying level, which attenuates every estimate above.
# The same table is used again in question 4, where it explains why the obvious
# within-person analysis fails.
say("Estimating the reliability of each factor from the A4 repeats")
# Glucose is included from version 1.2.1, now that A4/LEARN carries it.
REL_VARS <- intersect(c("egfr", "l_ck", "hb", "hba1c", "urea", "l_ast", "l_alt",
                        "albumin", "l_crp", "glucose"), names(lab))

# The occasion-to-occasion (error) variance of each factor, from the same
# random intercept and slope model that reliability() fits. It is needed for
# the conditional reliability below and is not returned by reliability().
error_variance <- function(value, time, id) {
  d <- tibble(x = as.numeric(value), t = as.numeric(time), id = id) %>%
    filter(!is.na(x), !is.na(t)) %>%
    group_by(id) %>% filter(n() >= 3) %>% ungroup()
  if (nrow(d) < 30 || dplyr::n_distinct(d$id) < 10) return(NA_real_)
  m <- tryCatch(lme4::lmer(x ~ t + (t | id), data = d,
                           control = lme4::lmerControl(calc.derivs = FALSE)),
                error = function(e) NULL)
  if (is.null(m)) return(NA_real_)
  vc <- as.data.frame(lme4::VarCorr(m))
  vc$vcov[nrow(vc)]
}

# Regression calibration with covariates in the model divides by the
# reliability of the factor given those covariates, 1 - error variance /
# var(factor | covariates), not by the marginal intraclass correlation: the
# covariates explain part of the true variation in the factor but none of its
# error, so the reliability that attenuates the coefficient is lower than the
# marginal one. var(factor | covariates) is the residual variance of the raw
# factor regressed on the model's covariates, in the rows the model used.
#
# In A4 the error variance is the one estimated from the A4 repeats. ADNI has no
# repeats, and its laboratory reports urea and albumin in different units from
# A4, so an error variance cannot be carried across in absolute terms. It is
# carried across as the same fraction of the total variance, (1 - ICC), which is
# what applying the A4 reliability to ADNI means (eMethods 8.4), giving
# 1 - (1 - ICC) x var(factor) / var(factor | covariates) in the ADNI sample.
conditional_rel <- function(data, measure_col, v, covs, id, cohort, icc, err_var) {
  if (is.na(icc) || !v %in% names(data)) return(NA_real_)
  f <- as.formula(sprintf("`%s` ~ z_%s + %s", measure_col, v, covs))
  d <- analysis_rows(f, data, id)
  d <- d[!is.na(d[[v]]), , drop = FALSE]
  if (nrow(d) < 30) return(NA_real_)
  v_cond <- summary(stats::lm(as.formula(sprintf("%s ~ %s", v, covs)), data = d))$sigma^2
  err <- if (cohort == "A4 screening") err_var else (1 - icc) * stats::var(d[[v]])
  1 - err / v_cond
}

R1 <- purrr::map_dfr(REL_VARS, function(v) {
  reliability(lab[[v]], lab$lab_time, lab$BID) %>%
    mutate(var = v, error_var = error_variance(lab[[v]], lab$lab_time, lab$BID), .before = 1)
})
# The conditional reliability depends on the analysis sample, so T1 carries the
# value used for each of its rows (rel_conditional). R1 gives it for one
# reference sample, the A4 screening p-tau181 sample, beside the unconditional
# intraclass correlation (icc_single), which step 06 and eTable 5 also use.
R1 <- R1 %>%
  mutate(rel_conditional = purrr::map_dbl(seq_len(n()), function(i)
    conditional_rel(a4, "pTau181", var[i], COVS_A4, "BID", "A4 screening",
                    icc_single[i], error_var[i]))) %>%
  left_join(FACTORS %>% select(var, label, role), by = "var") %>%
  select(var, label, role, n_people, rel_conditional, icc_single, error_var)
write_csv(R1, file.path(TAB_DIR, "R1_reliability.csv"))
print(as.data.frame(R1 %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

# The corrected estimate is what a perfectly measured factor would have given.
# The interval is scaled with it. That treats the reliability as known, which
# it is not, so the corrected interval is a little too narrow and is reported
# as a sensitivity analysis rather than as the main result. For the censored
# rows the coefficient corrected is that of the censored regression.
T1 <- T1 %>%
  left_join(R1 %>% select(var, icc_single, error_var), by = "var") %>%
  mutate(rel_conditional = purrr::pmap_dbl(
    list(cohort, measure_col, var, icc_single, error_var),
    function(co, ms, v, icc, ev) {
      if (co == "ADNI") conditional_rel(adni, ms, v, COVS_ADNI, "RID", co, icc, ev)
      else conditional_rel(a4, ms, v, COVS_A4, "BID", co, icc, ev)
    })) %>%
  mutate(b_corrected = rc_correct(b, rel_conditional),
         pct_corrected = pct(b_corrected),
         pct_corrected_lo = pct(rc_correct(lo, rel_conditional)),
         pct_corrected_hi = pct(rc_correct(hi, rel_conditional))) %>%
  select(-error_var, -measure_col)
write_csv(T1, file.path(TAB_DIR, "T1_levels.csv"))

say("Primary factors, per standard deviation, in both cohorts:")
print(as.data.frame(
  T1 %>% filter(role == "primary", measure %in% c("pTau181", "QX_ptau181", "Roche_ptau181",
                                                  "NfL", "Roche_NfL", "GFAP", "Roche_GFAP",
                                                  "C2N_ptau217")) %>%
    transmute(cohort, measure, factor = sub(" \\(.*", "", factor),
              estimate = sprintf("%.2f (%.2f to %.2f) q=%.3f", pct_diff, pct_lo, pct_hi, q)) %>%
    pivot_wider(names_from = factor, values_from = estimate)), row.names = FALSE)

# ---- eGFR is not linear ------------------------------------------------------
# A flexible curve against a straight line, and the difference from eGFR 90 at
# the two values that matter clinically.
#
# p_curved is a cluster-robust Wald test (CR1, clustered by person as in the
# level models, chi-squared on 3 degrees of freedom) that the three non-linear
# spline terms are zero. An F test from anova() would treat the three ADNI draws
# of each person as independent. The linear term plus the first three spline
# columns span the same space as the four spline columns, which is checked
# through the residual sum of squares, so testing those three columns is the
# test against a straight line.
#
# From version 1.3.0 the C2N p-tau217 spline is fitted by left-censored
# regression on every sample (spline_censored in 00_setup.R), as every other
# C2N p-tau217 estimate is; before, it used the samples above the detection
# limit only, which selects on the outcome. Every difference from eGFR 90 now
# carries a 95% interval, from the same cluster-robust variance as the test.
# These tests are reported
# without false discovery rate adjustment: they describe the shape of an
# association already tested in T1, and are not a family of their own.
egfr_shape <- function(data, measures, covs, id, cohort) {
  measures <- intersect(measures, names(data))
  purrr::map_dfr(measures, function(ms) {
    d <- as.data.frame(data) %>% filter(!is.na(.data[[ms]]), !is.na(egfr), !is.na(centiloid))
    if (ms == C2N_NUM) {
      r <- spline_censored(d, ms, covs, id, at = c(60, 45, 30))
      return(tibble(cohort = cohort, measure = short(ms), floor_handling = "left-censored regression",
                    n_obs = r$n_obs, n_people = r$n_people, p_curved = r$p_curved,
                    pct_at_60 = r$points$pct[1], pct_at_45 = r$points$pct[2],
                    pct_at_30 = r$points$pct[3],
                    pct_at_60_lo = r$points$lo[1], pct_at_60_hi = r$points$hi[1],
                    pct_at_45_lo = r$points$lo[2], pct_at_45_hi = r$points$hi[2],
                    pct_at_30_lo = r$points$lo[3], pct_at_30_hi = r$points$hi[3]))
    }
    handling <- "all samples"
    rownames(d) <- NULL
    if (n_distinct(d$egfr) < 50) return(NULL)
    f2 <- as.formula(sprintf("`%s` ~ splines::ns(egfr, knots = c(45, 60, 75)) + %s", ms, covs))
    m2 <- lm(f2, d)
    at <- function(x) { nd <- d[c(1, 1), ]; nd$egfr <- c(90, x); diff(predict(m2, nd)) }

    # The same fit with the last spline column replaced by eGFR itself.
    X2 <- model.matrix(m2)
    X2 <- X2[, !is.na(coef(m2)), drop = FALSE]
    rows <- as.integer(rownames(model.frame(m2)))
    y  <- model.response(model.frame(m2))
    ns_cols <- grep("^splines::ns\\(egfr", colnames(X2))
    X3 <- X2
    X3[, ns_cols[length(ns_cols)]] <- d$egfr[rows]
    colnames(X3)[ns_cols[length(ns_cols)]] <- "egfr_linear"
    m3 <- stats::lm.fit(X3, y)
    if (abs(sum(m3$residuals^2) - sum(residuals(m2)^2)) > 1e-6 * sum(residuals(m2)^2))
      stop("The reparameterised spline does not reproduce the spline fit for ", ms)
    test <- colnames(X3)[ns_cols[-length(ns_cols)]]
    cl <- d[[id]][rows]
    G <- length(unique(cl)); n <- nrow(X3); k <- ncol(X3)
    bread <- solve(crossprod(X3))
    U <- rowsum(X3 * m3$residuals, cl)
    V <- bread %*% crossprod(U) %*% bread * (G / (G - 1)) * ((n - 1) / (n - k))
    bt <- m3$coefficients[test]
    W  <- as.numeric(t(bt) %*% solve(V[test, test]) %*% bt)

    # The difference from eGFR 90 with its 95% interval, on the same
    # parameterisation and variance as the test. The basis is the one the
    # formula built, on the same rows, and the point estimate is checked
    # against predict() on the spline fit.
    B0 <- splines::ns(d$egfr, knots = c(45, 60, 75))
    keep_sp <- colnames(X3)[ns_cols[-length(ns_cols)]]
    ctr <- function(x) {
      cv <- stats::setNames(rep(0, ncol(X3)), colnames(X3))
      cv[keep_sp] <- (as.numeric(predict(B0, newx = x)) - as.numeric(predict(B0, newx = 90)))[seq_along(keep_sp)]
      cv["egfr_linear"] <- x - 90
      e <- sum(cv * m3$coefficients)
      s <- sqrt(as.numeric(t(cv) %*% V %*% cv))
      c(e, e - 1.96 * s, e + 1.96 * s)
    }
    c60 <- ctr(60); c45 <- ctr(45); c30 <- ctr(30)
    if (abs(c60[1] - at(60)) > 1e-6 || abs(c45[1] - at(45)) > 1e-6)
      stop("The interval contrast does not reproduce the spline fit for ", ms)

    tibble(cohort = cohort, measure = short(ms), floor_handling = handling,
           n_obs = nobs(m2), n_people = G,
           p_curved = stats::pchisq(W, df = length(test), lower.tail = FALSE),
           pct_at_60 = pct(at(60)), pct_at_45 = pct(at(45)), pct_at_30 = pct(at(30)),
           pct_at_60_lo = pct(c60[2]), pct_at_60_hi = pct(c60[3]),
           pct_at_45_lo = pct(c45[2]), pct_at_45_hi = pct(c45[3]),
           pct_at_30_lo = pct(c30[2]), pct_at_30_hi = pct(c30[3]))
  })
}
T1b <- bind_rows(
  # Every ADNI platform for the four biomarkers drawn in Figure 2 and
  # eFigure 1, so that each curve drawn has its row here.
  egfr_shape(adni, c("C2N_plasma_ptau217", "Fuji_plasma_ptau217", "AlzPath_plasma_ptau217",
                     "Janssen_plasma_ptau217",
                     "Roche_plasma_ptau181", "QX_plasma_ptau181",
                     "Roche_plasma_NfL", "QX_plasma_NfL",
                     "Roche_plasma_GFAP", "QX_plasma_GFAP"), COVS_ADNI, "RID", "ADNI"),
  egfr_shape(a4, c("pTau181", "NfL", "GFAP"), COVS_A4, "BID", "A4 screening"))
write_csv(T1b, file.path(TAB_DIR, "T1b_egfr_shape.csv"))

# ---- The primary factors in one model ---------------------------------------
mutual <- function(data, measures, covs, id, cohort) {
  vars <- FACTORS %>% filter(role == "primary", paste0("z_", var) %in% names(data)) %>% pull(var)
  zs <- paste0("z_", vars)
  measures <- intersect(measures, names(data))
  purrr::map_dfr(measures, function(ms) purrr::map_dfr(zs, function(zv) {
    f <- as.formula(sprintf("`%s` ~ %s + %s", ms, paste(zs, collapse = " + "), covs))
    fit_level(f, data, zv, cluster = id) %>%
      mutate(cohort = cohort, measure = short(ms), var = sub("^z_", "", zv), .before = 1)
  }))
}
T1c <- bind_rows(
  mutual(adni, c("Roche_plasma_ptau181", "QX_plasma_ptau181", "Roche_plasma_NfL"), COVS_ADNI, "RID", "ADNI"),
  mutual(a4, c("pTau181", "NfL", "GFAP"), COVS_A4, "BID", "A4 screening")) %>%
  mutate(pct_diff = pct(b), pct_lo = pct(lo), pct_hi = pct(hi))
write_csv(T1c, file.path(TAB_DIR, "T1c_mutual_adjustment.csv"))

# ---- Negative control --------------------------------------------------------
# A factor that acts on the assay or on clearance should not predict amyloid in
# the brain. ALT is included because it turns out not to be a clean control.
# The outcome is Centiloid, which has no detection limit, so least squares.
negative <- function(data, id, cohort) {
  factors <- FACTORS %>% filter(paste0("z_", var) %in% names(data))
  covs <- if (cohort == "ADNI") "age_at_visit + female + e4 + dx + plasma_time" else
    "age_at_visit + female + e4"
  purrr::map_dfr(seq_len(nrow(factors)), function(j) {
    zv <- paste0("z_", factors$var[j])
    fit_cl(as.formula(sprintf("centiloid ~ %s + %s", zv, covs)), data, zv, cluster = id) %>%
      mutate(cohort = cohort, factor = factors$label[j], var = factors$var[j],
             role = factors$role[j], .before = 1)
  })
}
T1d <- bind_rows(negative(adni, "RID", "ADNI"), negative(a4, "BID", "A4 screening")) %>%
  rename(centiloid_per_sd = b)
write_csv(T1d, file.path(TAB_DIR, "T1d_negative_control.csv"))

# ---- Creatine kinase in men and women ---------------------------------------
sex_split <- function(data, measures, covs, id, cohort) {
  measures <- intersect(measures, names(data))
  purrr::map_dfr(c(0, 1), function(f) purrr::map_dfr(measures, function(ms) {
    d <- data %>% filter(female == f)
    covs_nosex <- gsub(" \\+ female", "", covs)
    fit_level(as.formula(sprintf("`%s` ~ z_l_ck + %s", ms, covs_nosex)), d, "z_l_ck", cluster = id) %>%
      mutate(cohort = cohort, sex = if (f == 1) "Women" else "Men",
             measure = short(ms), .before = 1)
  }))
}
T1e <- bind_rows(
  sex_split(adni, c("Roche_plasma_ptau181", "QX_plasma_ptau181", "AlzPath_plasma_ptau217",
                    "C2N_plasma_ptau217"), COVS_ADNI, "RID", "ADNI"),
  sex_split(a4, c("pTau181", "NfL"), COVS_A4, "BID", "A4 screening")) %>%
  mutate(pct_diff = pct(b), pct_lo = pct(lo), pct_hi = pct(hi))
write_csv(T1e, file.path(TAB_DIR, "T1e_sex_stratified.csv"))

# ---- How old the ADNI laboratory values are ---------------------------------
# The ADNI panel is measured once, at screening, and the plasma series can come
# years later. This repeats the primary estimates keeping only the draws whose
# own factor value was measured within the window of that draw. The gap is
# worked out draw by draw, not from the first draw, because every draw enters
# the model: a person whose screening laboratory was a month before the first
# draw is six years from the third.
# Each factor is timed by its own measurement: eGFR by the date of the
# creatinine it is computed from, the laboratory panel by the laboratory date,
# and body mass index by the vitals nearest the draw, which are always within a
# year of it, so its restricted rows differ from 'any' only by rows without a
# vitals date.
# The 'any' row applies no restriction at all, so it is the main estimate of T1
# exactly, on the same sample.
draw_gap_years <- function(data, v) {
  switch(v,
         egfr = as.numeric(data$egfr_date - data$visit_date) / YEAR_DAYS,
         bmi  = data$bmi_gap_days / YEAR_DAYS,
         as.numeric(data$lab_date - data$visit_date) / YEAR_DAYS)
}
gap_source <- c(egfr = "creatinine date", bmi = "vitals date")
T1f_measures <- c("Roche_plasma_ptau181", "QX_plasma_ptau181", "Roche_plasma_NfL")
T1f_vars <- FACTORS %>% filter(role == "primary", paste0("z_", var) %in% names(adni)) %>% pull(var)
T1f <- purrr::map_dfr(c(1, 3, NA), function(lim) purrr::map_dfr(T1f_vars, function(v) {
  d <- if (is.na(lim)) adni else {
    g <- draw_gap_years(adni, v)
    adni[!is.na(g) & abs(g) <= lim, , drop = FALSE]
  }
  levels_for(d, T1f_measures, COVS_ADNI, "RID", "ADNI", vars = v) %>%
    select(-measure_col) %>%
    mutate(within_years = lim,
           gap_measured_from = if (is.na(lim)) "no restriction (main estimate)" else
             if (v %in% names(gap_source)) gap_source[[v]] else "laboratory date",
           .before = 1)
})) %>% mutate(pct_diff = pct(b), pct_lo = pct(lo), pct_hi = pct(hi))
write_csv(T1f, file.path(TAB_DIR, "T1f_lab_recency.csv"))

say("Written to %s", TAB_DIR)
say("T1 %d rows, T1b %d, T1c %d, T1d %d, T1e %d, T1f %d",
    nrow(T1), nrow(T1b), nrow(T1c), nrow(T1d), nrow(T1e), nrow(T1f))
})


# #############################################################################
# Section 3.2. Level models under six specifications (eTable 35)
# Developed as 14_review_checks.R version 1.1.0, section 7
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
T1   <- read_tab("T1_levels")

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

# From section 1b of 14_review_checks.R: the primary factors.
PRIMARY <- FACTORS %>% filter(role == "primary", paste0("z_", var) %in% names(adni))
# From section 5 of 14_review_checks.R: the Fujirebio columns of the cleared ratio.
FUJI_P217 <- "Fuji_plasma_ptau217"; FUJI_AB42 <- "Fuji_plasma_Ab42"

# =============================================================================
# 7. Level models: without diagnosis, with Centiloid as a spline, at the first
#    draw, and in cognitively unimpaired amyloid-negative participants
# =============================================================================
# The level model of step 03 for the four primary factors. C2N p-tau217 and the
# %p-tau217 ratio are left-censored as in step 03; every other measure is least
# squares, including the cleared ratio, here taken with the Fujirebio floor as
# supplied (eTable 32 shows that handling gives the same answer).
#   without diagnosis: diagnosis lies on a path from the factors through
#     vascular and kidney disease and from pathology, so it can act as a
#     collider; this is the check on keeping it
#   Centiloid as a natural spline (3 df), with and without diagnosis: amyloid
#     burden is related to p-tau nonlinearly, and part of what diagnosis
#     absorbs may be that
#   first draw: each participant's first draw with the measure, the draw
#     nearest the screening laboratory panel, without the time term
#   cognitively normal at the first draw with a Centiloid below 20 at the
#     draw: addresses pathology and prodromal weight loss, not collider bias,
#     since it also conditions on diagnosis
say("7. Level models: without diagnosis, Centiloid spline, first draw, CU amyloid-negative")
adni$cleared_ratio_ls <- adni[[FUJI_P217]] - adni[[FUJI_AB42]]
LEVEL_MEASURES <- c(P217, "C2N_plasma_nptau217", "C2N_plasma_ptau217_ratio", P181,
                    "Roche_plasma_NfL", "QX_plasma_NfL", "Roche_plasma_GFAP", "QX_plasma_GFAP",
                    "C2N_plasma_Abeta42_Abeta40", FUJI_AB42, "cleared_ratio_ls")
LEVEL_MEASURES <- intersect(LEVEL_MEASURES, names(adni))
C2N_RATIO <- "C2N_plasma_ptau217_ratio"; C2N_DEN <- "C2N_plasma_nptau217"
fit_level14 <- function(data, ms, zv, covs) {
  f <- stats::as.formula(sprintf("`%s` ~ %s + %s", ms, zv, covs))
  if (!ms %in% c(C2N_P217, C2N_RATIO)) {
    r <- fit_cl(f, data, zv)
    return(tibble(n_obs = r$n_obs, n_people = r$n_people, b = r$b, lo = r$lo, hi = r$hi,
                  method = "least squares"))
  }
  vars <- unique(c(all.vars(f), "RID", "c2n_at_floor", if (ms == C2N_RATIO) C2N_DEN))
  dd <- as.data.frame(data)
  dd <- dd[stats::complete.cases(dd[, vars, drop = FALSE]), , drop = FALSE]
  rownames(dd) <- NULL
  limit <- if (ms == C2N_RATIO) 100 * C2N_LOD / exp(dd[[C2N_DEN]]) else C2N_LOD
  r <- tryCatch(fit_tobit(f, dd, zv, censor_at = limit, censored = as.logical(dd$c2n_at_floor)),
                error = function(e) NULL)
  if (is.null(r)) return(tibble(n_obs = nrow(dd), n_people = n_distinct(dd$RID), b = NA_real_,
                                lo = NA_real_, hi = NA_real_, method = "censored regression"))
  tibble(n_obs = nrow(dd), n_people = n_distinct(dd$RID), b = r$b, lo = r$lo, hi = r$hi,
         method = "censored regression")
}
first_with <- function(data, ms) data %>% filter(!is.na(.data[[ms]])) %>%
  group_by(RID) %>% filter(visit_no == min(visit_no)) %>% ungroup()
cu_neg <- adni %>% filter(as.character(dx) == "CN", !is.na(centiloid), centiloid < 20)
ANALYSES <- list(
  main                 = list(data = function(ms) adni,                 covs = COVS_ADNI),
  without_diagnosis    = list(data = function(ms) adni,                 covs = covs_without("dx")),
  centiloid_spline     = list(data = function(ms) adni,                 covs = covs_spline_cl()),
  centiloid_spline_without_diagnosis = list(data = function(ms) adni,
                                            covs = covs_spline_cl(covs_without("dx"))),
  first_draw           = list(data = function(ms) first_with(adni, ms), covs = covs_without("plasma_time")),
  cu_amyloid_negative  = list(data = function(ms) cu_neg,               covs = covs_without("dx")))
R7 <- purrr::map_dfr(LEVEL_MEASURES, function(ms) purrr::map_dfr(PRIMARY$var, function(v) {
  zv <- paste0("z_", v)
  purrr::map_dfr(names(ANALYSES), function(an) {
    a <- ANALYSES[[an]]
    fit_level14(a$data(ms), ms, zv, a$covs) %>%
      mutate(analysis = an, measure = short(ms), var = v, .before = 1)
  })
}))
# The main rows reproduce T1 wherever T1 has the measure.
chk7 <- R7 %>% filter(analysis == "main") %>%
  inner_join(T1 %>% filter(cohort == "ADNI") %>% select(measure, var, b_T1 = b), by = c("measure", "var"))
if (nrow(chk7) == 0 || any(abs(chk7$b - chk7$b_T1) > 1e-6, na.rm = TRUE))
  stop("Step 14 level models do not reproduce T1_levels.csv.")
say("  check: %d main level estimates reproduce T1", nrow(chk7))
R7 <- R7 %>% mutate(pct_diff = pct(b), pct_lo = pct(lo), pct_hi = pct(hi)) %>%
  left_join(R7 %>% filter(analysis == "main") %>% transmute(measure, var, main_pct = pct(b)),
            by = c("measure", "var"))
put_rev(R7, "review_7_level_sensitivity")
print(as.data.frame(R7 %>% filter(var %in% c("egfr", "bmi")) %>%
  transmute(measure, var, analysis, n_people, estimate = fmt(pct_diff, pct_lo, pct_hi))),
  row.names = FALSE)
})
