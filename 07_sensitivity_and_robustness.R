# =============================================================================
# 07_sensitivity_and_robustness.R
# Stage 7. Sensitivity analyses and robustness checks (eMethods 11 and 12).
#   7.1  classification agreement, a visit effect shared by both series,
#        creatinine in place of eGFR, E-values (eTables 20 to 25)
#   7.2  hypertension and diabetes added as covariates (eTable 29)
#   7.3  agreement with amyloid PET per SD of each factor and by band
#        (eTables 39, 40)
# Writes paper1_out/tables/T7 to T11, review/review_6*.csv, manuscript/S19 to S25, S27.
# =============================================================================

# #############################################################################
# Section 7.1. Classification, visit effect, creatinine, E-values
# Developed as 10_robustness.R version 1.1.3, whole
# One change: the visit-effect model (eTable 25) is fitted for eGFR only.
# #############################################################################
local({
# =============================================================================
# 10_robustness.R
# Four checks that answer the four ways a reviewer can argue the main results
# are an artefact rather than a finding.
#
#   1. Classification. A percentage difference on the log scale does not tell
#      anyone how often two assays would disagree about a patient. This turns
#      the between-platform differences into the proportion of samples that two
#      p-tau217 platforms place on opposite sides of an amyloid threshold, and
#      asks whether that proportion moves with each factor.
#
#   2. A shared disturbance at one visit. Creatinine and p-tau217 are measured
#      at the same visits, and the joint growth model has no term for anything
#      that moves both at a single visit, such as an acute illness or a sample
#      handled differently. Such a disturbance would be absorbed into the
#      correlation between the two rates of change and inflate it. The model is
#      refitted with a random intercept at the visit shared by both series, so
#      that the correlation is estimated free of it.
#
#   3. Creatinine in place of eGFR. eGFR is computed from creatinine using age
#      and sex, and the models also adjust for age and sex, so the exposure and
#      two covariates are linked by construction. Repeating the level models on
#      creatinine itself removes that link.
#
#   4. E-values. The level associations are between people and cannot be freed
#      of confounding by adjustment. An E-value states how strongly an
#      unmeasured confounder would have to be associated with both the factor
#      and the biomarker to explain the association away, which is a quantity a
#      reader can judge rather than a caveat they must accept.
#
# Nothing here changes any earlier result. Every check is reported beside the
# estimate it is checking, so a reader can see both.
#
# Reads    paper1_out/adni_visit.rds, a4_screen.rds, a4_long.rds, a4_lab.rds
#          paper1_out/tables/T1_levels.csv, T5b_joint_growth.csv,
#          T5d_span_sensitivity.csv (from step 06)
# Writes   paper1_out/tables/
#            T7_classification.csv      discordance and what moves it
#            T7b_thresholds.csv         the threshold used for each platform
#            T8_evalues.csv             E-values for the level associations
#            T9_creatinine.csv          creatinine beside eGFR
#            T10_occasion_check.csv     the joint model with a visit effect
#          paper1_out/manuscript/
#            S19.csv to S25.csv         the same, formatted; S25 is the span
#                                       sensitivity of the joint model, which
#                                       step 06 computes
#
# C2N p-tau217 has 424 of 1,179 samples below its detection limit, reported as
# half the limit. Wherever this step fits a level model to it, in checks 3 and
# 4, it uses the censored regression, censored at the limit itself, as step 05
# does, so that no substituted value is treated as a measurement.
#
# This step formats its own supplementary tables rather than adding to step 09,
# so that it can be run and rerun without touching anything step 09 produces.
#
# Needs    nlme for check 2. Takes about fifteen minutes, almost all of it in
#          check 2, which fits two mixed models per factor.
# =============================================================================
SCRIPT_VERSION <- "1.1.3"
cat("10_robustness.R version", SCRIPT_VERSION, "\n")

if (!exists("PAPER1_VERSION")) source("00_setup.R")

MS_DIR <- file.path(OUT_DIR, "manuscript")
dir.create(MS_DIR, showWarnings = FALSE, recursive = TRUE)

put <- function(x, name, what) {
  if (is.null(x) || nrow(x) == 0) { say("  %s skipped, nothing to write", name); return(invisible(NULL)) }
  write_csv(x, file.path(MS_DIR, paste0(name, ".csv")))
  say("  %-8s %-58s %3d rows", paste0(name, ".csv"), what, nrow(x))
}
short <- function(x) gsub("_plasma", "", x)
short_factor <- function(x) sub(" \\(.*", "", x)
fmt <- function(e, lo, hi) sprintf("%.2f (%.2f to %.2f)", e, lo, hi)

adni <- readRDS(file.path(OUT_DIR, "adni_visit.rds"))
a4   <- readRDS(file.path(OUT_DIR, "a4_screen.rds"))

PRIMARY <- FACTORS %>% filter(role == "primary")

# =============================================================================
# 1. Classification: how often do two p-tau217 platforms disagree about a
#    person, and does the disagreement track a physiological factor?
# =============================================================================
# Each platform is given the threshold that best separates amyloid PET positive
# from negative in this sample, by Youden's index. A threshold chosen this way
# flatters every platform equally, which is what is wanted here: the question is
# not which platform classifies best, it is whether two platforms that have each
# been given their best threshold still disagree, and whether the disagreement
# is patterned by something physiological rather than random.
say("1. Classification agreement between p-tau217 platforms")

CENTILOID_CUTS <- c(20, 24)   # amyloid PET positive, two thresholds in common use

PTAU217 <- intersect(c("C2N_plasma_ptau217", "Fuji_plasma_ptau217",
                       "AlzPath_plasma_ptau217", "Janssen_plasma_ptau217"),
                     names(adni))

# The value of x that maximises sensitivity plus specificity against pos.
youden_threshold <- function(x, pos) {
  ok <- !is.na(x) & !is.na(pos)
  x <- x[ok]; pos <- pos[ok]
  if (length(x) < 20 || length(unique(pos)) < 2) return(NA_real_)
  cand <- sort(unique(x))
  j <- vapply(cand, function(k) mean(x[pos] >= k) - mean(x[!pos] >= k), numeric(1))
  cand[which.max(j)]
}

classification_for <- function(d, cut, sample_label) {
  d <- d[!is.na(d$centiloid), , drop = FALSE]
  if (nrow(d) < 50) return(NULL)
  d$amy <- d$centiloid >= cut
  ms <- PTAU217[vapply(PTAU217, function(m) sum(!is.na(d[[m]])) >= 50, logical(1))]
  if (length(ms) < 2) return(NULL)
  ths <- vapply(ms, function(m) youden_threshold(d[[m]], d$amy), numeric(1))

  thresholds <<- bind_rows(thresholds, tibble(
    sample = sample_label, centiloid_cut = cut, platform = short(ms),
    threshold_pg_ml = exp(unname(ths)),
    n_samples = vapply(ms, function(m) sum(!is.na(d[[m]])), integer(1)),
    called_positive_pct = vapply(ms, function(m)
      100 * mean(d[[m]] >= ths[[m]], na.rm = TRUE), numeric(1)),
    pet_positive_pct = 100 * mean(d$amy)))

  purrr::map_dfr(combn(ms, 2, simplify = FALSE), function(pr) {
    a <- d[[pr[1]]] >= ths[[pr[1]]]
    b <- d[[pr[2]]] >= ths[[pr[2]]]
    ok <- !is.na(a) & !is.na(b)
    if (sum(ok) < 50) return(NULL)
    dd <- d[ok, , drop = FALSE]
    a <- a[ok]; b <- b[ok]
    dd$disagree <- as.numeric(a != b)
    # Positive on the first platform only counts as 1, on the second only as -1,
    # so a positive coefficient means the factor pushes the first platform up
    # relative to the second.
    dd$net <- as.numeric(a & !b) - as.numeric(b & !a)
    # Two samples are involved and both are reported. The disagreement
    # percentage is over every sample that both platforms measured (the pair
    # sample). Each per-SD model then also needs the factor and the
    # covariates, so it uses fewer samples and people (the model sample).
    # Before version 1.1.0 S19 paired the pair's samples with the model's
    # people, which described neither.
    base <- tibble(sample = sample_label, centiloid_cut = cut,
                   pair = sprintf("%s minus %s", short(pr[1]), short(pr[2])),
                   n_samples = nrow(dd), n_people_pair = n_distinct(dd$RID),
                   disagree_pct = 100 * mean(dd$disagree))
    purrr::map_dfr(seq_len(nrow(PRIMARY)), function(j) {
      zv <- paste0("z_", PRIMARY$var[j])
      if (!zv %in% names(dd)) return(NULL)
      r1 <- fit_cl(as.formula(sprintf("disagree ~ %s + %s", zv, COVS_ADNI)), dd, zv, cluster = "RID")
      r2 <- fit_cl(as.formula(sprintf("net ~ %s + %s", zv, COVS_ADNI)), dd, zv, cluster = "RID")
      base %>% mutate(
        factor = PRIMARY$label[j], var = PRIMARY$var[j],
        n_samples_model = r1$n_obs, n_people = r1$n_people,
        # The outcome is a proportion, so the coefficient is already a
        # difference in probability. It is reported in percentage points.
        disagree_pp = 100 * r1$b, disagree_lo = 100 * r1$lo, disagree_hi = 100 * r1$hi,
        disagree_p = r1$p,
        net_pp = 100 * r2$b, net_lo = 100 * r2$lo, net_hi = 100 * r2$hi, net_p = r2$p)
    })
  })
}

thresholds <- tibble()
adni_nofloor <- if ("c2n_at_floor" %in% names(adni))
  adni[!coalesce(adni$c2n_at_floor, FALSE), , drop = FALSE] else NULL

T7 <- bind_rows(
  purrr::map_dfr(CENTILOID_CUTS, function(k) classification_for(adni, k, "all samples")),
  if (!is.null(adni_nofloor) && nrow(adni_nofloor) >= 50)
    purrr::map_dfr(CENTILOID_CUTS, function(k)
      classification_for(adni_nofloor, k, "C2N floor excluded")) else NULL
)
if (nrow(T7) > 0) {
  # False discovery rate control within each sample and threshold, so that the
  # sensitivity runs cannot affect the adjusted values of the main one.
  T7 <- T7 %>% group_by(sample, centiloid_cut) %>%
    mutate(disagree_q = p.adjust(disagree_p, "BH"),
           net_q = p.adjust(net_p, "BH")) %>% ungroup()
  write_csv(T7, file.path(TAB_DIR, "T7_classification.csv"))
  write_csv(thresholds, file.path(TAB_DIR, "T7b_thresholds.csv"))
  say("Thresholds and how often each platform calls a sample positive:")
  print(as.data.frame(thresholds %>% filter(sample == "all samples", centiloid_cut == 20) %>%
    transmute(platform, threshold = round(threshold_pg_ml, 2), n_samples,
              called_positive = round(called_positive_pct, 2),
              pet_positive = round(pet_positive_pct, 2))), row.names = FALSE)
  say("How often two platforms disagree, and what moves the disagreement:")
  print(as.data.frame(T7 %>%
    filter(sample == "all samples", centiloid_cut == 20, var == "l_ck") %>%
    transmute(pair, n = n_samples, disagree = round(disagree_pct, 2),
              per_sd_ck = fmt(disagree_pp, disagree_lo, disagree_hi))), row.names = FALSE)
} else say("  no classification results, so this check is skipped")

# =============================================================================
# 2. A disturbance shared by both series at one visit
# =============================================================================
say("2. The joint growth model with a visit-level effect shared by both series")

L   <- if (file.exists(file.path(OUT_DIR, "a4_long.rds")))
  readRDS(file.path(OUT_DIR, "a4_long.rds")) else NULL
lab <- if (file.exists(file.path(OUT_DIR, "a4_lab.rds")))
  readRDS(file.path(OUT_DIR, "a4_lab.rds")) else NULL

TIME_CENTRE <- 2.5      # as in step 06
MIN_VALUES <- 2
MIN_SPAN_YEARS <- 1
# A laboratory value and a blood draw this close together are treated as the
# same visit. Five weeks is wide enough to catch a panel drawn either side of a
# study visit and narrow enough that two real visits are not merged.
OCC_TOLERANCE <- 0.10
# eGFR only: the creatine kinase and haemoglobin rows were removed from eTable 25
# on 16 Nov 2026 (decisions log 8.3).
OCC_VARS <- strsplit(Sys.getenv("PAPER1_OCC_VARS", unset = "egfr"), ",")[[1]]

T10 <- NULL
if (!is.null(L) && !is.null(lab) && requireNamespace("nlme", quietly = TRUE)) {

  eligible_people <- function(v) {
    f <- lab %>% filter(!is.na(.data[[v]])) %>%
      group_by(BID) %>%
      summarise(n = n(), span = max(lab_time) - min(lab_time), .groups = "drop") %>%
      filter(n >= MIN_VALUES, span >= MIN_SPAN_YEARS)
    b <- L %>% filter(!is.na(pTau217)) %>% count(BID) %>% filter(n >= MIN_VALUES)
    intersect(f$BID, b$BID)
  }
  # The same person-level covariates as the joint model of step 06 (from its
  # version 1.4.0): age at the first draw and Centiloid standardised, sex and
  # APOE4, each on the level and the rate of change of both series.
  JOINT_COVS <- c("age_z", "female", "e4", "cl_z")
  person_covs <- L %>% distinct(BID, age0, female, e4, centiloid) %>%
    group_by(BID) %>% slice(1) %>% ungroup()
  JOINT_FIXED <- as.formula(paste("y ~ 0 +", paste(c("ix", "ip", "ix:t", "ip:t",
    as.vector(outer(c("ix", "ip", "ix:t", "ip:t"), JOINT_COVS, paste, sep = ":"))),
    collapse = " + ")))
  joint_data <- function(v) {
    keep <- eligible_people(v)
    cv <- person_covs %>% filter(BID %in% keep) %>%
      filter(!is.na(age0), !is.na(female), !is.na(e4), !is.na(centiloid)) %>%
      mutate(age_z = as.numeric(scale(age0)), cl_z = as.numeric(scale(centiloid))) %>%
      select(BID, all_of(JOINT_COVS))
    keep <- intersect(keep, cv$BID)
    bio <- L %>% filter(!is.na(pTau217), BID %in% keep) %>%
      transmute(BID, t = plasma_time, y = pTau217, o = "bio")
    fac <- lab %>% filter(!is.na(.data[[v]]), BID %in% keep) %>%
      transmute(BID, t = lab_time, y = FACTOR_SIGN[[v]] * .data[[v]], o = "fac")
    bind_rows(bio, fac) %>%
      group_by(o) %>% mutate(y = as.numeric(scale(y))) %>% ungroup() %>%
      left_join(cv, by = "BID") %>%
      mutate(t = t - TIME_CENTRE,
             ip = as.integer(o == "bio"), ix = as.integer(o == "fac")) %>%
      arrange(BID, o, t)
  }
  # The covariance of the person-level random effects. getVarCov cannot return
  # one when a model has more than one level of nesting, so the matrix is taken
  # from the fitted structure itself, which works whatever the model's shape.
  person_cov <- function(m, group = "BID") {
    if (is.null(m)) return(NULL)
    rs <- m$modelStruct$reStruct
    pd <- if (group %in% names(rs)) rs[[group]] else rs[[length(rs)]]
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
  ctrl <- nlme::lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200,
                           returnObject = TRUE, apVar = FALSE)
  # With returnObject = TRUE nlme hands back a fit even when the optimiser
  # stopped without converging, and says so only in a warning. The warnings
  # are therefore caught rather than suppressed, and a fit that did not
  # converge is marked, so that it is not reported as if it had (step 06 does
  # the same). The transient "singular precision matrix" warnings are not a
  # sign of non-convergence and are ignored.
  NONCONVERGENCE_10 <- "convergence error code|without convergence"
  lme_checked <- function(...) {
    warned <- character()
    m <- withCallingHandlers(tryCatch(nlme::lme(...), error = function(e) NULL),
      warning = function(w) {
        warned <<- c(warned, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    if (!is.null(m)) attr(m, "converged") <- !any(grepl(NONCONVERGENCE_10, warned))
    m
  }
  converged_10 <- function(m) !is.null(m) && isTRUE(attr(m, "converged"))
  fit_plain <- function(d)
    lme_checked(JOINT_FIXED,
                random = list(BID = nlme::pdSymm(~ 0 + ix + ip + ix:t + ip:t)),
                weights = nlme::varIdent(form = ~ 1 | o),
                data = d, control = ctrl, method = "REML")
  # The same model with one addition: the residuals of the two series are
  # allowed to correlate within a visit. Anything that moves creatinine and
  # p-tau217 together at a single visit is then carried by that correlation
  # instead of by the two rates of change. The estimated correlation is itself
  # the answer to whether the concern is real, so it is reported.
  #
  # A shared random intercept at the visit expresses the same idea and is kept
  # as a fallback, because the two parameterisations do not always converge on
  # the same data.
  fit_with_visit <- function(d, optimiser = "optim") {
    ct <- nlme::lmeControl(opt = optimiser, maxIter = 200, msMaxIter = 200,
                           returnObject = TRUE, apVar = FALSE)
    # A fit that did not converge is passed over in the same way as one that
    # failed, so the fallback is tried rather than an unconverged answer used.
    m <- lme_checked(JOINT_FIXED,
                     random = list(BID = nlme::pdSymm(~ 0 + ix + ip + ix:t + ip:t)),
                     correlation = nlme::corCompSymm(form = ~ 1 | BID/occ),
                     weights = nlme::varIdent(form = ~ 1 | o),
                     data = d, control = ct, method = "REML")
    if (converged_10(m)) return(list(fit = m, how = "correlated residuals within visit"))
    m <- lme_checked(JOINT_FIXED,
                     random = list(BID = nlme::pdSymm(~ 0 + ix + ip + ix:t + ip:t),
                                   occ = nlme::pdIdent(~ 1)),
                     weights = nlme::varIdent(form = ~ 1 | o),
                     data = d, control = ct, method = "REML")
    if (converged_10(m)) return(list(fit = m, how = "shared visit intercept"))
    list(fit = NULL, how = "did not fit or did not converge")
  }
  within_visit_rho <- function(m) {
    if (is.null(m) || is.null(m$modelStruct$corStruct)) return(NA_real_)
    tryCatch(unname(coef(m$modelStruct$corStruct, unconstrained = FALSE))[1],
             error = function(e) NA_real_)
  }

  T10 <- purrr::map_dfr(intersect(OCC_VARS, FACTORS$var), function(v) {
    if (!v %in% names(lab)) return(NULL)
    d <- joint_data(v)
    if (n_distinct(d$BID) < 100) return(NULL)
    d$occ <- factor(paste(d$BID, round(d$t / OCC_TOLERANCE)))
    d <- d %>% arrange(BID, occ)
    shared <- d %>% group_by(occ) %>% summarise(both = n_distinct(o) == 2, .groups = "drop")
    t0 <- Sys.time()
    m0 <- fit_plain(d)
    r1 <- fit_with_visit(d)
    # A second optimiser, only if the first did not produce a usable answer.
    if (is.na(slope_correlation(r1$fit))) r1 <- fit_with_visit(d, "nlminb")
    say("  %s: %d of %d visits carry both series, fitted by %s, %.0f seconds",
        v, sum(shared$both), nrow(shared), r1$how,
        as.numeric(difftime(Sys.time(), t0, units = "secs")))
    tibble(var = v, n_people = n_distinct(d$BID),
           visits = nrow(shared), visits_with_both = sum(shared$both),
           correlation = if (converged_10(m0)) slope_correlation(m0) else NA_real_,
           converged = converged_10(m0),
           correlation_with_visit_effect = slope_correlation(r1$fit),
           within_visit_residual_correlation = within_visit_rho(r1$fit),
           fitted_by = r1$how,
           seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  })
  if (!is.null(T10) && nrow(T10) > 0) {
    write_csv(T10, file.path(TAB_DIR, "T10_occasion_check.csv"))
    print(as.data.frame(T10 %>% transmute(
      factor = var, people = n_people, visits_with_both,
      without = round(correlation, 2),
      with_visit_effect = round(correlation_with_visit_effect, 2),
      within_visit_r = round(within_visit_residual_correlation, 2), fitted_by)),
      row.names = FALSE)
  }
} else say("  A4 longitudinal data or nlme not available, so this check is skipped")

# =============================================================================
# 3. Creatinine in place of eGFR
# =============================================================================
say("3. The level models on creatinine itself rather than eGFR")

add_creatinine <- function(d) {
  if (!"creatinine" %in% names(d)) return(d)
  d$l_creatinine <- safe_log(d$creatinine)
  # Higher creatinine means worse kidney function, which is the direction eGFR
  # is signed in, so no sign flip is needed for the two to be comparable.
  d$z_l_creatinine <- zf(d$l_creatinine)
  d
}
adni <- add_creatinine(adni)
a4   <- add_creatinine(a4)

ADNI_MAIN <- intersect(c("C2N_plasma_ptau217", "AlzPath_plasma_ptau217",
                         "Roche_plasma_ptau181", "QX_plasma_ptau181",
                         "Roche_plasma_NfL", "QX_plasma_NfL",
                         "Roche_plasma_GFAP", "QX_plasma_GFAP",
                         "C2N_plasma_Abeta42_Abeta40"), names(adni))
A4_MAIN <- intersect(c("pTau181", "NfL", "GFAP", "Ab42_40"), names(a4))

# C2N p-tau217 below its limit is reported as half the limit, and 424 of 1,179
# samples are. A linear model treats that substituted value as a measurement,
# which eMethods 5 rules out, so this measure is fitted by censored regression,
# censored at the limit itself, exactly as step 05 fits it. Every other
# measure here is quantified throughout and keeps the linear model.
CENSORED_MEASURES <- intersect("C2N_plasma_ptau217", c(ADNI_MAIN, A4_MAIN))

# The rows a level model actually uses: the outcome, the exposure, every
# covariate and the person identifier all present. People are counted and the
# outcome's standard deviation is taken on these rows, so that both describe
# the sample that produced the estimate. Passing only these rows to the
# censored regression also keeps its own counts right.
analysis_rows <- function(formula, d, id) {
  d <- as.data.frame(d)
  vars <- unique(c(all.vars(formula), id))
  d[stats::complete.cases(d[, vars, drop = FALSE]), , drop = FALSE]
}

# The standard deviation of the outcome in the analysis sample. For a censored
# outcome the coefficient is on the scale of the uncensored (latent) value, and
# the standard deviation of the observed values would be that of a variable
# with 424 samples replaced by one constant. The matching quantity is the
# latent standard deviation, which is the scale of a censored regression with
# no predictors, fitted on the same rows.
outcome_sd_of <- function(y, censored) {
  if (!censored) return(stats::sd(y))
  lim <- log(C2N_LOD)
  obs <- y > lim + 1e-8
  fit <- survival::survreg(survival::Surv(ifelse(obs, y, lim), obs, type = "left") ~ 1,
                           dist = "gaussian")
  fit$scale
}

# One level model, by the method the measure needs, with its analysis sample.
level_fit <- function(d, ms, zv, covs, id) {
  f  <- as.formula(sprintf("`%s` ~ %s + %s", ms, zv, covs))
  dd <- analysis_rows(f, d, id)
  cens <- ms %in% CENSORED_MEASURES
  r <- if (cens) fit_tobit(f, dd, zv, censor_at = C2N_LOD, cluster = id) else
    fit_cl(f, dd, zv, cluster = id)
  tibble(method = if (cens) "censored regression" else "linear regression",
         n_obs = nrow(dd), n_people = n_distinct(dd[[id]]),
         b = r$b, lo = r$lo, hi = r$hi,
         outcome_sd = outcome_sd_of(dd[[ms]], cens),
         outcome_sd_basis = if (cens) "latent, censored model, analysis sample" else
           "observed, analysis sample")
}

creatinine_for <- function(d, measures, covs, id, cohort) {
  if (!"z_l_creatinine" %in% names(d)) return(NULL)
  purrr::map_dfr(measures, function(ms) {
    r <- level_fit(d, ms, "z_l_creatinine", covs, id)
    out <- tibble(cohort = cohort, measure = short(ms), method = r$method,
                  n_people = r$n_people, creatinine_pct = pct(r$b),
                  creatinine_lo = pct(r$lo), creatinine_hi = pct(r$hi),
                  egfr_pct = NA_real_, egfr_lo = NA_real_, egfr_hi = NA_real_)
    # A censored creatinine estimate has to be read beside a censored eGFR
    # estimate, so that one is fitted here, by the same method. The linear
    # eGFR estimates are taken from T1 below.
    if (ms %in% CENSORED_MEASURES && "z_egfr" %in% names(d)) {
      e <- level_fit(d, ms, "z_egfr", covs, id)
      out <- out %>% mutate(egfr_pct = pct(e$b), egfr_lo = pct(e$lo), egfr_hi = pct(e$hi))
    }
    out
  })
}
T9 <- bind_rows(
  creatinine_for(adni, ADNI_MAIN, COVS_ADNI, "RID", "ADNI"),
  creatinine_for(a4,   A4_MAIN,   COVS_A4,   "BID", "A4 screening"))

# Beside the eGFR estimate the paper reports, so the two can be compared.
T1r <- if (file.exists(file.path(TAB_DIR, "T1_levels.csv")))
  read_csv(file.path(TAB_DIR, "T1_levels.csv"), show_col_types = FALSE) else NULL
if (!is.null(T1r) && nrow(T9) > 0)
  T9 <- T9 %>% left_join(
    T1r %>% filter(var == "egfr") %>%
      transmute(cohort, measure, .e = pct_diff, .elo = pct_lo, .ehi = pct_hi),
    by = c("cohort", "measure")) %>%
    mutate(egfr_pct = coalesce(egfr_pct, .e), egfr_lo = coalesce(egfr_lo, .elo),
           egfr_hi = coalesce(egfr_hi, .ehi)) %>%
    select(-.e, -.elo, -.ehi)
if (nrow(T9) > 0) {
  write_csv(T9, file.path(TAB_DIR, "T9_creatinine.csv"))
  print(as.data.frame(T9 %>% transmute(
    cohort, measure, method, people = n_people,
    creatinine = fmt(creatinine_pct, creatinine_lo, creatinine_hi),
    egfr = fmt(egfr_pct, egfr_lo, egfr_hi))),
    row.names = FALSE)
}

# =============================================================================
# 4. E-values for the level associations
# =============================================================================
# How strongly would an unmeasured confounder have to be associated with both
# the factor and the biomarker, over and above the covariates already in the
# model, to move the estimate to the null? The outcome is continuous, so the
# estimate is first put on an approximate risk ratio scale by standardising it
# against the spread of the outcome.
#
# eMethods 12.2 takes that spread in the analysis sample. Before version 1.1.0
# it was taken over the whole file, which includes rows the model never saw.
# Each model is therefore refitted here on its own rows, which gives the
# standard deviation on the right sample and, for C2N p-tau217, the censored
# estimate. For every other measure the refit is the T1 model itself, and the
# two are checked against each other.
say("4. E-values for the level associations")

evalue <- function(rr) {
  rr <- ifelse(!is.na(rr) & rr < 1, 1 / rr, rr)
  ifelse(is.na(rr), NA_real_, rr + sqrt(rr * (rr - 1)))
}

T8 <- NULL
if (!is.null(T1r)) {
  MAIN <- bind_rows(
    tibble(cohort = "ADNI", ms = ADNI_MAIN, covs = COVS_ADNI, id = "RID"),
    tibble(cohort = "A4 screening", ms = A4_MAIN, covs = COVS_A4, id = "BID")) %>%
    mutate(measure = short(ms))
  rows <- T1r %>% filter(role == "primary") %>%
    inner_join(MAIN, by = c("cohort", "measure"))
  fits <- purrr::map_dfr(seq_len(nrow(rows)), function(i) {
    r <- rows[i, ]
    d <- if (r$cohort == "ADNI") adni else a4
    level_fit(d, r$ms, paste0("z_", r$var), r$covs, r$id) %>%
      rename_with(~ paste0("e_", .x))
  })
  rows <- bind_cols(rows, fits)
  # For a linear model the refit must reproduce T1 exactly; if it does not, T1
  # was built from different data and the E-values would not match it.
  lin <- rows %>% filter(e_method == "linear regression")
  if (nrow(lin) && (max(abs(lin$e_b - lin$b), na.rm = TRUE) > 1e-8 ||
                    any(lin$e_n_obs != lin$n_obs)))
    say("  WARNING: the refitted level models do not reproduce T1_levels.csv. Rerun step 03.")
  cen <- rows %>% filter(e_method == "censored regression")
  if (nrow(cen) && max(abs(cen$e_b - cen$b), na.rm = TRUE) > 1e-6)
    say("  note: T1_levels.csv holds an estimate for %s that is not the censored one; the E-values use the censored estimate.",
        paste(unique(cen$measure), collapse = ", "))
  T8 <- rows %>%
    mutate(rr    = exp(0.91 * e_b / e_outcome_sd),
           rr_lo = exp(0.91 * e_lo / e_outcome_sd),
           rr_hi = exp(0.91 * e_hi / e_outcome_sd),
           e_value = evalue(rr),
           # The limit of the interval nearer the null. An interval that spans
           # the null has an E-value of 1, meaning no confounder is needed.
           e_value_interval = dplyr::case_when(
             is.na(rr_lo) | is.na(rr_hi) ~ NA_real_,
             rr_lo <= 1 & rr_hi >= 1 ~ 1,
             rr_lo > 1 ~ evalue(rr_lo),
             TRUE ~ evalue(rr_hi))) %>%
    transmute(cohort, measure, factor = short_factor(factor), role,
              method = e_method, n_people = e_n_people, n_obs = e_n_obs,
              estimate = pct(e_b), lo = pct(e_lo), hi = pct(e_hi),
              outcome_sd = e_outcome_sd, outcome_sd_basis = e_outcome_sd_basis,
              risk_ratio = rr, e_value, e_value_interval)
  write_csv(T8, file.path(TAB_DIR, "T8_evalues.csv"))
  print(as.data.frame(T8 %>% filter(factor == "eGFR") %>%
    transmute(cohort, measure, method, estimate = fmt(estimate, lo, hi),
              outcome_sd = round(outcome_sd, 3),
              e_value = round(e_value, 2),
              e_value_interval = round(e_value_interval, 2))), row.names = FALSE)
  say("  eGFR E-values range from %.2f to %.2f; all primary level associations from %.2f to %.2f.",
      min(T8$e_value[T8$factor == "eGFR"], na.rm = TRUE),
      max(T8$e_value[T8$factor == "eGFR"], na.rm = TRUE),
      min(T8$e_value, na.rm = TRUE), max(T8$e_value, na.rm = TRUE))
} else say("  T1_levels.csv is not there, so the E-values are skipped")

# =============================================================================
# The supplementary tables for this step
# =============================================================================
say("Building the supplementary tables for this step")

if (nrow(T7) > 0) {
  put(T7 %>% filter(sample == "all samples", centiloid_cut == 20) %>%
        transmute(pair, pair_samples = n_samples, pair_people = n_people_pair,
                  disagree_pct = round(disagree_pct, 2),
                  factor = short_factor(factor),
                  # The sample of the per-SD models, which need the factor and
                  # the covariates as well as both platforms.
                  model_samples = n_samples_model, model_people = n_people,
                  disagreement_per_sd = fmt(disagree_pp, disagree_lo, disagree_hi),
                  disagreement_q = signif(disagree_q, 2),
                  direction_per_sd = fmt(net_pp, net_lo, net_hi),
                  direction_q = signif(net_q, 2)),
      "S19", "how often two p-tau217 platforms disagree about a person")
  put(T7 %>% filter(!(sample == "all samples" & centiloid_cut == 20)) %>%
        transmute(sample, centiloid_cut, pair, pair_samples = n_samples,
                  disagree_pct = round(disagree_pct, 2), factor = short_factor(factor),
                  model_samples = n_samples_model, model_people = n_people,
                  disagreement_per_sd = fmt(disagree_pp, disagree_lo, disagree_hi)),
      "S20", "classification disagreement, other thresholds and samples")
}
if (!is.null(T8))
  put(T8 %>% transmute(cohort, measure, factor, method, people = n_people,
                       estimate = fmt(estimate, lo, hi),
                       e_value = round(e_value, 2),
                       e_value_interval = round(e_value_interval, 2)),
      "S21", "E-values for the level associations")
if (nrow(T9) > 0)
  put(T9 %>% transmute(cohort, measure, method, people = n_people,
                       creatinine = fmt(creatinine_pct, creatinine_lo, creatinine_hi),
                       egfr = if ("egfr_pct" %in% names(T9))
                         fmt(egfr_pct, egfr_lo, egfr_hi) else NA_character_),
      "S22", "levels on creatinine itself beside the eGFR estimate")
# The three ways of handling the C2N detection limit, for the ratio and its two
# components. Computed in step 05; formatted here so that Table 4 has a
# supplementary table to point to.
T4b <- if (file.exists(file.path(TAB_DIR, "T4b_ratio_floor.csv")))
  read_csv(file.path(TAB_DIR, "T4b_ratio_floor.csv"), show_col_types = FALSE) else NULL
if (!is.null(T4b) && nrow(T4b) > 0)
  # Each handling carries its estimate with the 95% interval and its false
  # discovery rate adjusted value, so that the claim that a conclusion holds
  # under every handling can be checked from the table itself.
  put(T4b %>% transmute(factor = short_factor(factor),
                        measure = short(measure),
                        samples = n_all, samples_after_exclusion = n_after_exclusion,
                        as_delivered = fmt(pct_as_delivered, lo_as_delivered, hi_as_delivered),
                        q_as_delivered = signif(q_as_delivered, 2),
                        floor_excluded = fmt(pct_floor_excluded, lo_floor_excluded, hi_floor_excluded),
                        q_floor_excluded = signif(q_floor_excluded, 2),
                        censored = fmt(pct_censored, lo_censored, hi_censored),
                        q_censored = signif(q_censored, 2)),
      "S24", "the detection limit and the ratio, three ways of handling it")

if (!is.null(T10) && nrow(T10) > 0)
  put(T10 %>% left_join(FACTORS %>% select(var, label), by = "var") %>%
        transmute(factor = short_factor(label), people = n_people,
                  visits, visits_with_both,
                  correlation = round(correlation, 2),
                  with_visit_effect = round(correlation_with_visit_effect, 2),
                  within_visit_residual_correlation =
                    round(within_visit_residual_correlation, 2),
                  fitted_by),
      "S23", "change on change with a visit effect shared by both series")

# The joint model under the two other span rules, computed in step 06. The
# main rule is the first block, so the table reads as the main result and the
# two departures from it.
T5d <- if (file.exists(file.path(TAB_DIR, "T5d_span_sensitivity.csv")))
  read_csv(file.path(TAB_DIR, "T5d_span_sensitivity.csv"), show_col_types = FALSE) else NULL
if (!is.null(T5d) && nrow(T5d) > 0)
  put(T5d %>% left_join(FACTORS %>% select(var, label), by = "var") %>%
        transmute(rule = description, factor = short_factor(label),
                  people = n_people, factor_values = n_factor_values,
                  ptau217_draws = n_draws,
                  correlation = sprintf("%.2f", slope_correlation),
                  interval = ifelse(is.na(lo), "", sprintf("%.2f to %.2f", lo, hi)),
                  bootstrap_draws = ifelse(is.na(draws), "", as.character(draws)),
                  usable_draws = ifelse(is.na(usable), "", as.character(usable))),
      "S25", "change on change under two other span rules, with the out-of-range p-tau217 draws excluded, and with a quadratic term for time")

say("Written to %s and %s", TAB_DIR, MS_DIR)
say("Nothing in this step changes an earlier result. Each check is reported")
say("beside the estimate it is checking.")
})


# #############################################################################
# Section 7.2. Hypertension and diabetes (eTable 29)
# Developed as 12_comorbidity_and_labels.R version 1.0.4, sections 1 to 2c
# One change: the joint growth model is fitted for eGFR only.
# #############################################################################
local({
# =============================================================================
# 12_comorbidity_and_labels.R
# Paper 1: the comorbidity sensitivity analysis.
#
# (Sections 2d, 3 and 4 of this script are now sections 5.2 and 9.4.)
#
#   1  Hypertension and diabetes as additional covariates. The level
#      associations are comparisons between people, and hypertension and
#      diabetes both lower eGFR and are linked to higher NfL and GFAP through
#      small-vessel disease. The change analysis has the same exposure: faster
#      kidney decline and a faster rise in p-tau217 could both reflect vascular
#      ageing. This step refits
#        a  the level models of the four primary factors with the main
#           biomarkers (the E-value set of step 10),
#        b  the two creatine kinase differences, C2N minus ALZpath for
#           p-tau217 and Roche minus Quanterix for p-tau181,
#        c  the joint growth model of step 06, for eGFR (the creatine kinase
#           and haemoglobin rows were removed on 16 Nov 2026),
#      each twice on the same rows: without and with hypertension and diabetes.
#      Comparing the two on the same rows isolates the adjustment from the
#      change of sample. The main estimate of the paper is reported beside
#      them, from the table of the step that produced it.
#
#      This is a sensitivity analysis and replaces nothing. Hypertension lies
#      partly on the path from kidney function to the biomarkers, so adjusting
#      for it can remove part of a genuine kidney effect as well as confounding.
#      A confounder acting on the circulating concentration cancels from a
#      same-aliquot difference, so part b is expected to move little; it is
#      run because it answers the comorbidity reading of the subsamples.
#
# Definitions
#   ADNI hypertension  recorded in the medical history (RECMHIST), or an
#                      antihypertensive in the medication file, dated on or
#                      before the draw or undated, or blood pressure of 140/90
#                      or above at the vitals visit nearest the draw, within a
#                      year. 0 when the history and the medication record both
#                      exist and neither, nor blood pressure, is positive.
#   ADNI diabetes      recorded in the medical history, or a glucose-lowering
#                      drug, on the same dating rule. 0 when both records exist
#                      and neither is positive. ADNI glucose is not used, as the
#                      release does not say whether it is fasting.
#   A4/LEARN           the release has no medical history or medication file.
#                      Hypertension is screening blood pressure of 140/90 or
#                      above (visit 001) and diabetes HbA1c of 6.5% or above,
#                      from the laboratory value nearest the screening draw
#                      (a4_screen) or the earliest laboratory visit (joint
#                      model). Treated hypertension with normal blood pressure
#                      is therefore counted as absent in A4/LEARN.
#   Every matched history and medication string is written to checks/ so that
#   the definitions can be inspected, as step 01 does for statins.
#
# Reads    paper1_out/adni_visit.rds, a4_screen.rds, a4_long.rds, a4_lab.rds
#          paper1_out/tables/T1_levels.csv, T2_platform_difference.csv,
#            T5b_joint_growth.csv
#          ADNI RECMHIST.csv, RECCMEDS.csv, VITALS.csv; A4 Raw Data/vitals.csv
#          All raw files read only.
# Writes   paper1_out/tables/T11_comorbidity.csv
#          paper1_out/checks/comorbidity_build.csv, adni_hypertension_strings.csv,
#            adni_diabetes_strings.csv, adni_antihypertensive_strings.csv,
#            adni_glucose_lowering_strings.csv
#          paper1_out/manuscript/S27.csv  the comorbidity sensitivity analysis
#
# The joint model is fitted without a bootstrap by default: the question is
# whether the correlation moves, which the point estimate answers. To add an
# interval for the adjusted eGFR correlation, set before sourcing the launcher
#      Sys.setenv(PAPER1_COMORB_BOOT = 1000)
# which takes two hours or more (the unadjusted 1,000 draws took 117 minutes).
#
# Needs    survival and nlme. Takes a few minutes without the bootstrap.
# =============================================================================
# Changes in 1.0.4: the between-platform verdicts are relabelled "equivalent"
# and "differ", as in the manuscript; an earlier bootstrap interval for the
# adjusted eGFR correlation is kept when the step is rerun without the
# bootstrap. No estimate changes.
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

adni <- readRDS(file.path(OUT_DIR, "adni_visit.rds"))
a4   <- readRDS(file.path(OUT_DIR, "a4_screen.rds"))
L    <- readRDS(file.path(OUT_DIR, "a4_long.rds"))
lab  <- readRDS(file.path(OUT_DIR, "a4_lab.rds"))
PRIMARY <- FACTORS %>% filter(role == "primary")

# =============================================================================
# 1. Hypertension and diabetes
# =============================================================================
say("1. Building hypertension and diabetes")

# ---- ADNI: recorded medical history -----------------------------------------
# RECMHIST.csv is registered in 00_setup.R (ADNI_FILES$mhist).
mh_raw  <- read_source(adni_file("mhist"), "history")
mh_date <- pick_col(mh_raw, c("VISDATE", "EXAMDATE", "USERDATE"))
mh <- mh_raw %>%
  transmute(RID = as.integer(RID),
            mh_date = parse_date_ymd(.data[[mh_date]]),
            desc = str_to_lower(str_squish(coalesce(MHDESC, "")))) %>%
  filter(!is.na(RID)) %>%
  semi_join(adni %>% distinct(RID), by = "RID")

HTN_PATTERN     <- "hypertens|\\bhtn\\b|high blood pressure|elevated blood pressure|\\bhbp\\b"
NOT_HTN_PATTERN <- paste0("pre-?hypertens|pre hypertens|borderline|pulmonary|ocular|intraocular|intracranial|portal|white[- ]coat|",
                          "gestational|pregnan|eclamp|family hist|hypotens")
# Borderline hypertension is excluded, as borderline diabetes is below, so that
# the two conditions are defined the same way (version 1.0.3).
DM_PATTERN      <- "diabet|\\bdm\\b|\\bt2dm\\b|\\bt1dm\\b|\\bniddm\\b|\\biddm\\b|\\baodm\\b"
NOT_DM_PATTERN  <- paste0("insipidus|pre-?diabet|pre diabet|borderline|gestational|",
                          "family hist|risk of|impaired|intoleran")
mh <- mh %>%
  mutate(is_htn = str_detect(desc, HTN_PATTERN) & !str_detect(desc, NOT_HTN_PATTERN),
         is_dm  = str_detect(desc, DM_PATTERN)  & !str_detect(desc, NOT_DM_PATTERN))

# ---- ADNI: medication --------------------------------------------------------
# Generic and common brand names, listed rather than matched by suffix, as in
# step 01. Topical minoxidil and ophthalmic timolol are not antihypertensives
# and are not listed. Alpha blockers (doxazosin, terazosin, prazosin) are left
# out because in older men they are mostly prescribed for prostate symptoms.
# Short brand names are matched as whole words, so that "actos" does not match
# "lactose". Beta blockers and some diuretics are also prescribed for
# other reasons (tremor, arrhythmia, heart failure), which this definition
# cannot separate; the strings are written out for inspection.
ANTIHTN_PATTERN <- paste0(
  # ACE inhibitors, angiotensin receptor blockers, and their combinations
  "lisinopril|enalapril|ramipril|benazepril|quinapril|fosinopril|captopril|perindopril|",
  "trandolapril|moexipril|zestril|prinivil|zestoretic|vasotec|altace|lotensin|lotrel|accupril|",
  "monopril|aceon|quinopril|linsinopril|",
  "losartan|valsartan|irbesartan|candesartan|olmesartan|telmisartan|azilsartan|eprosartan|",
  "cozaar|hyzaar|diovan|exforge|avapro|avalide|benicar|\\bazor\\b|micardis|atacand|edarbi|temisartan|",
  # diuretics
  "furosemide|frusemide|lasix|bumetanide|bumex|torsemide|torasemide|demadex|",
  "hydrochlorothiazide|hctz|chlorthalidone|chlortalidone|indapamide|metolazone|chlorothiazide|",
  "spironolactone|aldactone|eplerenone|inspra|amiloride|triamterene|dyazide|maxzide|microzide|",
  "thiazide|thizide|thazide|triazide|athyzide|hydrochlorathiazide|hydrochlorothyazide|hydrocholoro|",
  # calcium channel blockers
  "amlodipine|norvasc|nifedipine|procardia|adalat|felodipine|plendil|isradipine|nicardipine|",
  "diltiazem|cardizem|tiazac|verapamil|\\bcalan\\b|isoptin|",
  # beta blockers
  "metoprolol|toprol|lopressor|atenolol|tenormin|carvedilol|\\bcoreg\\b|bisoprolol|\\bziac\\b|nebivolol|",
  "bystolic|labetalol|propranolol|inderal|nadolol|",
  # others
  "clonidine|catapres|hydralazine|methyldopa")
GLUCOSE_PATTERN <- paste0(
  "metformin|glucophage|glipizide|glucotrol|glyburide|glibenclamide|glimepiride|amaryl|",
  "sitagliptin|januvia|janumet|saxagliptin|onglyza|linagliptin|tradjenta|alogliptin|",
  "pioglitazone|\\bactos\\b|rosiglitazone|avandia|repaglinide|prandin|nateglinide|starlix|acarbose|",
  "insulin|lantus|levemir|humalog|novolog|humulin|novolin|toujeo|tresiba|basaglar|apidra|",
  "liraglutide|victoza|dulaglutide|trulicity|semaglutide|ozempic|rybelsus|exenatide|byetta|bydureon|",
  "empagliflozin|jardiance|dapagliflozin|farxiga|canagliflozin|invokana|invokamet|ertugliflozin")

med_raw <- read_source(adni_file("meds"), "medication")
m_mdate <- pick_col(med_raw, c("VISDATE", "EXAMDATE", "USERDATE"))
med <- med_raw %>%
  transmute(RID = as.integer(RID),
            med_date = parse_date_ymd(.data[[m_mdate]]),
            med = str_to_lower(coalesce(CMMED, ""))) %>%
  filter(!is.na(RID)) %>%
  semi_join(adni %>% distinct(RID), by = "RID") %>%
  mutate(is_antihtn = str_detect(med, ANTIHTN_PATTERN),
         is_glucose = str_detect(med, GLUCOSE_PATTERN))

strings_out <- function(d, flag, text, file) {
  s <- d %>% filter(.data[[flag]]) %>% count(.data[[text]], name = "records") %>%
    arrange(desc(records))
  names(s)[1] <- "text"
  write_csv(s, file.path(CHECK_DIR, file))
  s
}
s_htn  <- strings_out(mh,  "is_htn",     "desc", "adni_hypertension_strings.csv")
s_dm   <- strings_out(mh,  "is_dm",      "desc", "adni_diabetes_strings.csv")
s_ahtn <- strings_out(med, "is_antihtn", "med",  "adni_antihypertensive_strings.csv")
s_glu  <- strings_out(med, "is_glucose", "med",  "adni_glucose_lowering_strings.csv")
say("  ADNI strings matched: hypertension %d, diabetes %d, antihypertensive %d, glucose-lowering %d (checks/)",
    nrow(s_htn), nrow(s_dm), nrow(s_ahtn), nrow(s_glu))

# ---- ADNI: blood pressure at the vitals visit nearest each draw -------------
vt_raw <- read_source(adni_file("vitals"), "vitals")
v_date <- pick_col(vt_raw, c("EXAMDATE", "VISDATE"))
bp <- vt_raw %>%
  transmute(RID = as.integer(RID),
            exam_date = parse_date_ymd(.data[[v_date]]),
            sbp = suppressWarnings(as.numeric(VSBPSYS)),
            dbp = suppressWarnings(as.numeric(VSBPDIA))) %>%
  mutate(sbp = if_else(sbp >= 60 & sbp <= 260, sbp, NA_real_),
         dbp = if_else(dbp >= 30 & dbp <= 150, dbp, NA_real_)) %>%
  filter(!is.na(RID), !is.na(exam_date), !is.na(sbp) | !is.na(dbp)) %>%
  semi_join(adni %>% distinct(RID), by = "RID") %>%
  distinct(RID, exam_date, .keep_all = TRUE) %>%
  nearest_visit("exam_date", adni %>% select(RID, visit_no, visit_date), window = 365) %>%
  transmute(RID, visit_no, bp_high = as.integer(coalesce(sbp >= 140, FALSE) | coalesce(dbp >= 90, FALSE)))

# ---- ADNI: at each draw ------------------------------------------------------
# A history or medication record counts at a draw when it is dated on or before
# that draw, or undated, the rule of step 01.
at_draw <- function(d, flag, date_col, name) {
  d %>% filter(.data[[flag]]) %>%
    inner_join(adni %>% select(RID, visit_no, visit_date), by = "RID",
               relationship = "many-to-many") %>%
    group_by(RID, visit_no) %>%
    summarise(!!name := as.integer(any(is.na(.data[[date_col]]) | .data[[date_col]] <= visit_date)),
              .groups = "drop")
}
has_hist <- mh  %>% distinct(RID) %>% mutate(hist_recorded = 1L)
has_med  <- med %>% distinct(RID) %>% mutate(med_recorded = 1L)
adni_cm <- adni %>%
  select(RID, visit_no) %>%
  left_join(at_draw(mh,  "is_htn",     "mh_date",  "htn_hist"), by = c("RID", "visit_no")) %>%
  left_join(at_draw(mh,  "is_dm",      "mh_date",  "dm_hist"),  by = c("RID", "visit_no")) %>%
  left_join(at_draw(med, "is_antihtn", "med_date", "htn_med"),  by = c("RID", "visit_no")) %>%
  left_join(at_draw(med, "is_glucose", "med_date", "dm_med"),   by = c("RID", "visit_no")) %>%
  left_join(bp, by = c("RID", "visit_no")) %>%
  left_join(has_hist, by = "RID") %>%
  left_join(has_med,  by = "RID") %>%
  mutate(across(c(htn_hist, dm_hist, htn_med, dm_med), ~ coalesce(.x, 0L)),
         both_records = !is.na(hist_recorded) & !is.na(med_recorded),
         htn = case_when(htn_hist == 1 | htn_med == 1 | coalesce(bp_high, 0L) == 1 ~ 1L,
                         both_records ~ 0L, TRUE ~ NA_integer_),
         dm  = case_when(dm_hist == 1 | dm_med == 1 ~ 1L,
                         both_records ~ 0L, TRUE ~ NA_integer_)) %>%
  select(RID, visit_no, htn, dm, htn_hist, htn_med, bp_high, dm_hist, dm_med)
adni <- adni %>% left_join(adni_cm, by = c("RID", "visit_no"))

# ---- A4/LEARN: screening blood pressure and HbA1c ----------------------------
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
a4 <- a4 %>%
  left_join(a4_vt, by = "BID") %>%
  mutate(dm = if_else(is.na(hba1c), NA_integer_, as.integer(hba1c >= 6.5)))
# For the joint model, one value per person: screening blood pressure, and the
# HbA1c at the person's earliest laboratory visit that has one.
a4_person_cm <- lab %>%
  filter(!is.na(hba1c)) %>%
  group_by(BID) %>% slice_min(day, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(BID, dm = as.integer(hba1c >= 6.5)) %>%
  full_join(a4_vt, by = "BID")

# ---- What was built ----------------------------------------------------------
first_draw <- adni %>% group_by(RID) %>% slice_min(visit_no, n = 1, with_ties = FALSE) %>% ungroup()
prev <- function(x) sprintf("%d of %d (%.2f%%)", sum(x == 1, na.rm = TRUE), sum(!is.na(x)),
                            100 * mean(x == 1, na.rm = TRUE))
chk <- tibble::tribble(
  ~item, ~value,
  "ADNI people, hypertension at the first draw", prev(first_draw$htn),
  "  from the medical history", prev(first_draw$htn_hist),
  "  from an antihypertensive", prev(first_draw$htn_med),
  "  from blood pressure 140/90 or above", prev(first_draw$bp_high),
  "ADNI people, diabetes at the first draw", prev(first_draw$dm),
  "  from the medical history", prev(first_draw$dm_hist),
  "  from a glucose-lowering drug", prev(first_draw$dm_med),
  "ADNI draws, hypertension", prev(adni$htn),
  "ADNI draws, diabetes", prev(adni$dm),
  "A4/LEARN screening, hypertension by blood pressure", prev(a4$htn),
  "A4/LEARN screening, diabetes by HbA1c", prev(a4$dm),
  "A4/LEARN people with screening blood pressure and an HbA1c",
    as.character(sum(!is.na(a4_person_cm$htn) & !is.na(a4_person_cm$dm)))
)
print(as.data.frame(chk), row.names = FALSE)
write_csv(chk, file.path(CHECK_DIR, "comorbidity_build.csv"))

CM_COVS <- "htn + dm"

# =============================================================================
# 2a. Level models, with and without hypertension and diabetes
# =============================================================================
say("2a. Level models")
ADNI_MAIN <- intersect(c("C2N_plasma_ptau217", "AlzPath_plasma_ptau217",
                         "Roche_plasma_ptau181", "QX_plasma_ptau181",
                         "Roche_plasma_NfL", "QX_plasma_NfL",
                         "Roche_plasma_GFAP", "QX_plasma_GFAP",
                         "C2N_plasma_Abeta42_Abeta40"), names(adni))
A4_MAIN <- intersect(c("pTau181", "NfL", "GFAP", "Ab42_40"), names(a4))
C2N_NUM <- "C2N_plasma_ptau217"

# One level model on given rows. C2N p-tau217 is the censored regression of
# step 03, censored at the limit with the laboratory's flag; everything else
# least squares with the CR1 sandwich.
level_on <- function(d, ms, zv, covs, id) {
  f <- as.formula(sprintf("`%s` ~ %s + %s", ms, zv, covs))
  if (ms != C2N_NUM) return(fit_cl(f, d, zv, cluster = id) %>% mutate(method = "least squares"))
  r <- tryCatch(fit_tobit(f, d, zv, censor_at = C2N_LOD, cluster = id,
                          censored = coalesce(as.logical(d$c2n_at_floor), FALSE)),
                error = function(e) NULL)
  if (is.null(r)) return(tibble(n_obs = nrow(d), n_people = NA_integer_, b = NA_real_,
                                lo = NA_real_, hi = NA_real_, p = NA_real_,
                                method = "censored regression"))
  r %>% mutate(method = "censored regression")
}
level_pair <- function(data, ms, v, covs, id, cohort) {
  zv <- paste0("z_", v)
  need <- unique(c(ms, zv, all.vars(as.formula(paste("~", covs))), "htn", "dm", id,
                   if (ms == C2N_NUM) "c2n_at_floor"))
  d <- as.data.frame(data)
  d <- d[stats::complete.cases(d[, intersect(need, names(d)), drop = FALSE]), , drop = FALSE]
  rownames(d) <- NULL
  if (nrow(d) < 30) return(NULL)
  r0 <- level_on(d, ms, zv, covs, id)
  r1 <- level_on(d, ms, zv, paste(covs, "+", CM_COVS), id)
  tibble(analysis = "level", cohort = cohort, analyte = NA_character_, measure = short(ms), var = v,
         method = r0$method, n_obs = nrow(d), n_people = n_distinct(d[[id]]),
         pct0 = pct(r0$b), lo0 = pct(r0$lo), hi0 = pct(r0$hi),
         pct1 = pct(r1$b), lo1 = pct(r1$lo), hi1 = pct(r1$hi))
}
LV <- bind_rows(
  purrr::map_dfr(ADNI_MAIN, function(ms) purrr::map_dfr(PRIMARY$var, function(v)
    if (paste0("z_", v) %in% names(adni)) level_pair(adni, ms, v, COVS_ADNI, "RID", "ADNI"))),
  purrr::map_dfr(A4_MAIN, function(ms) purrr::map_dfr(PRIMARY$var, function(v)
    if (paste0("z_", v) %in% names(a4)) level_pair(a4, ms, v, COVS_A4, "BID", "A4 screening"))))
T1r <- read_csv(file.path(TAB_DIR, "T1_levels.csv"), show_col_types = FALSE)
LV <- LV %>% left_join(T1r %>% transmute(cohort, measure, var, main_pct = pct_diff,
                                         main_lo = pct_lo, main_hi = pct_hi, main_n = n_people),
                       by = c("cohort", "measure", "var"))

# =============================================================================
# 2b. The two creatine kinase differences
# =============================================================================
say("2b. Creatine kinase differences between platforms")
# The same fit as step 04: C2N minus the other platform, left-censored at
# log(1.30) minus the other platform's value where C2N is below the limit.
diff_on <- function(d, pr, zv, covs) {
  f <- as.formula(sprintf("diff ~ %s + %s", zv, covs))
  if (pr[1] != C2N_NUM) return(fit_cl(f, d, zv))
  tryCatch(fit_tobit(f, d, zv, censor_at = C2N_LOD / exp(d[[pr[2]]]),
                     censored = coalesce(as.logical(d$c2n_at_floor), FALSE)),
           error = function(e) tibble(b = NA_real_, lo = NA_real_, hi = NA_real_))
}
DIFF_PAIRS <- list(c("C2N_plasma_ptau217", "AlzPath_plasma_ptau217"),
                   c("Roche_plasma_ptau181", "QX_plasma_ptau181"))
DF <- purrr::map_dfr(DIFF_PAIRS, function(pr) {
  if (!all(pr %in% names(adni))) return(NULL)
  zv <- "z_l_ck"
  d <- as.data.frame(adni)
  d$diff <- d[[pr[1]]] - d[[pr[2]]]
  need <- unique(c(pr, "diff", zv, all.vars(as.formula(paste("~", COVS_ADNI))), "htn", "dm", "RID",
                   if (pr[1] == C2N_NUM) "c2n_at_floor"))
  d <- d[stats::complete.cases(d[, need, drop = FALSE]), , drop = FALSE]
  rownames(d) <- NULL
  r0 <- diff_on(d, pr, zv, COVS_ADNI)
  r1 <- diff_on(d, pr, zv, paste(COVS_ADNI, "+", CM_COVS))
  tibble(analysis = "difference", cohort = "ADNI",
         analyte = if (pr[1] == C2N_NUM) "pTau217" else "pTau181",
         measure = paste(short(pr[1]), "minus", short(pr[2])), var = "l_ck",
         method = if (pr[1] == C2N_NUM) "censored regression" else "least squares",
         n_obs = nrow(d), n_people = n_distinct(d$RID),
         pct0 = pct(r0$b), lo0 = pct(r0$lo), hi0 = pct(r0$hi),
         pct1 = pct(r1$b), lo1 = pct(r1$lo), hi1 = pct(r1$hi))
})
T2r <- read_csv(file.path(TAB_DIR, "T2_platform_difference.csv"), show_col_types = FALSE)
DF <- DF %>% left_join(T2r %>% filter(var == "l_ck") %>%
                         transmute(measure = pair, main_pct = pct_diff, main_lo = pct_lo,
                                   main_hi = pct_hi, main_n = n_people),
                       by = "measure")

# =============================================================================
# 2c. The joint growth model
# =============================================================================
# The model of step 06, section 3, restated here so that step 06 and its
# bootstraps do not have to be rerun. Covariates enter the level and the rate
# of change of both series. The unadjusted fit on every eligible person is
# checked against T5b_joint_growth.csv before anything else is reported.
say("2c. The joint growth model")
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
person_covs <- L %>% distinct(BID, age0, female, e4, centiloid) %>%
  group_by(BID) %>% slice(1) %>% ungroup() %>%
  left_join(a4_person_cm, by = "BID")
joint_fixed <- function(extra = character()) {
  base <- c("ix", "ip", "ix:t", "ip:t")
  adj  <- as.vector(outer(base, c(JOINT_COVS, extra), paste, sep = ":"))
  as.formula(paste("y ~ 0 +", paste(c(base, adj), collapse = " + ")))
}
joint_data <- function(v, need_cm) {
  keep <- eligible_people(v)
  cv <- person_covs %>% filter(BID %in% keep) %>%
    filter(!is.na(age0), !is.na(female), !is.na(e4), !is.na(centiloid))
  if (need_cm) cv <- cv %>% filter(!is.na(htn), !is.na(dm))
  cv <- cv %>% mutate(age_z = as.numeric(scale(age0)), cl_z = as.numeric(scale(centiloid))) %>%
    select(BID, all_of(JOINT_COVS), if (need_cm) c("htn", "dm"))
  keep <- intersect(keep, cv$BID)
  bio <- L %>% filter(!is.na(pTau217), BID %in% keep) %>%
    transmute(BID, t = plasma_time, y = pTau217, o = "bio")
  fac <- lab %>% filter(!is.na(.data[[v]]), BID %in% keep) %>%
    transmute(BID, t = lab_time, y = FACTOR_SIGN[[v]] * .data[[v]], o = "fac")
  bind_rows(bio, fac) %>%
    group_by(o) %>% mutate(y = as.numeric(scale(y))) %>% ungroup() %>%
    left_join(cv, by = "BID") %>%
    mutate(t = t - TIME_CENTRE, ip = as.integer(o == "bio"), ix = as.integer(o == "fac")) %>%
    arrange(BID, o, t)
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

T5b <- read_csv(file.path(TAB_DIR, "T5b_joint_growth.csv"), show_col_types = FALSE)
check_d <- joint_data("egfr", need_cm = FALSE)
check_r <- slope_correlation(fit_joint(check_d))
main_r  <- T5b$slope_correlation[T5b$var == "egfr" & T5b$group == "everyone"]
say("  check: unadjusted eGFR correlation %.4f here, %.4f in T5b_joint_growth.csv",
    check_r, main_r)
if (length(main_r) != 1 || is.na(check_r) || abs(check_r - main_r) > 1e-4)
  stop("The restated joint model does not reproduce step 06. Nothing further is reported.")

JM <- purrr::map_dfr(intersect("egfr", names(lab)), function(v) {
  d  <- joint_data(v, need_cm = TRUE)
  m0 <- fit_joint(d)
  m1 <- fit_joint(d, extra = c("htn", "dm"))
  m1b <- fit_joint(d, extra = c("htn", "dm"), optimiser = "nlminb")
  tibble(analysis = "joint growth model", cohort = "A4", analyte = "pTau217",
         measure = "p-tau217 with the factor", var = v,
         method = "bivariate linear mixed model", n_obs = sum(d$o == "fac"),
         n_draws = sum(d$o == "bio"), n_people = n_distinct(d$BID),
         r0 = slope_correlation(m0), r1 = slope_correlation(m1),
         r1_other_optimiser = slope_correlation(m1b),
         converged0 = converged(m0), converged1 = converged(m1),
         main_r = T5b$slope_correlation[T5b$var == v & T5b$group == "everyone"][1],
         main_n = T5b$n_people[T5b$var == v & T5b$group == "everyone"][1])
})

BOOT <- suppressWarnings(as.integer(Sys.getenv("PAPER1_COMORB_BOOT", "0")))
if (is.na(BOOT)) BOOT <- 0
JM$lo1 <- NA_real_; JM$hi1 <- NA_real_; JM$usable1 <- NA_integer_
if (BOOT > 0) {
  say("  bootstrapping the adjusted eGFR correlation, %d draws. This is the slow part.", BOOT)
  d <- joint_data("egfr", need_cm = TRUE)
  people <- unique(d$BID)
  set.seed(1)
  t0 <- Sys.time()
  rs <- vapply(seq_len(BOOT), function(i) {
    if (i %% 50 == 0) say("  egfr adjusted: %d of %d, %.0f minutes so far", i, BOOT,
                          as.numeric(difftime(Sys.time(), t0, units = "mins")))
    ids <- sample(people, length(people), replace = TRUE)
    dd <- tibble(BID = ids, .new = seq_along(ids)) %>%
      inner_join(d, by = "BID", relationship = "many-to-many") %>%
      mutate(BID = .new) %>% select(-.new)
    m <- tryCatch(fit_joint(dd, extra = c("htn", "dm")), error = function(e) NULL)
    if (!converged(m)) NA_real_ else slope_correlation(m)
  }, numeric(1))
  JM <- JM %>% mutate(lo1 = if_else(var == "egfr", unname(quantile(rs, 0.025, na.rm = TRUE)), lo1),
                      hi1 = if_else(var == "egfr", unname(quantile(rs, 0.975, na.rm = TRUE)), hi1),
                      usable1 = if_else(var == "egfr", sum(!is.na(rs)), usable1))
}

# An interval bootstrapped in an earlier run is kept when this run does not
# bootstrap, provided the adjusted fit is the same one: same people and the same
# correlation to 1e-8. Without this a quick rerun, for the labels alone, would
# silently drop the interval (version 1.0.4).
if (BOOT == 0 && file.exists(file.path(TAB_DIR, "T11_comorbidity.csv"))) {
  prev <- read_csv(file.path(TAB_DIR, "T11_comorbidity.csv"), show_col_types = FALSE) %>%
    filter(analysis == "joint growth model", !is.na(lo1)) %>%
    select(var, n_people_prev = n_people, r1_prev = r1, lo1_prev = lo1, hi1_prev = hi1,
           usable1_prev = usable1)
  if (nrow(prev) > 0) {
    JM <- JM %>% left_join(prev, by = "var") %>%
      mutate(same = !is.na(r1_prev) & n_people == n_people_prev & abs(r1 - r1_prev) < 1e-8,
             lo1 = if_else(same, lo1_prev, lo1), hi1 = if_else(same, hi1_prev, hi1),
             usable1 = if_else(same, as.integer(usable1_prev), usable1))
    kept <- JM$var[JM$same]
    JM <- JM %>% select(-ends_with("_prev"), -same)
    if (length(kept)) say("  bootstrap not run; kept the earlier interval for %s", paste(kept, collapse = ", "))
    else say("  bootstrap not run; the earlier interval belongs to a different fit and is not reused")
  }
}

T11 <- bind_rows(LV, DF, JM) %>%
  left_join(FACTORS %>% select(var, label), by = "var")
write_csv(T11, file.path(TAB_DIR, "T11_comorbidity.csv"))

say("How much the estimates move with hypertension and diabetes added, on the same rows:")
print(as.data.frame(T11 %>% filter(analysis != "joint growth model") %>%
  transmute(analysis, cohort, measure, factor = short_factor(label), people = n_people,
            without = fmt(pct0, lo0, hi0), with = fmt(pct1, lo1, hi1))), row.names = FALSE)
print(as.data.frame(JM %>% transmute(factor = var, people = n_people,
                                     without = round(r0, 3), with = round(r1, 3),
                                     other_optimiser = round(r1_other_optimiser, 3),
                                     paper = round(main_r, 3), converged1)), row.names = FALSE)

# The formatted table, eTable 29 in the supplement.
S27 <- bind_rows(
  T11 %>% filter(analysis != "joint growth model") %>%
    transmute(analysis, cohort, analyte = coalesce(analyte, ""), measure, factor = short_factor(label), method,
              people_main = as.character(main_n), main_estimate = fmt(main_pct, main_lo, main_hi),
              people = as.character(n_people), without_adjustment = fmt(pct0, lo0, hi0),
              with_adjustment = fmt(pct1, lo1, hi1)),
  T11 %>% filter(analysis == "joint growth model") %>%
    transmute(analysis, cohort, analyte, measure, factor = short_factor(label), method,
              people_main = as.character(main_n), main_estimate = sprintf("%.2f", main_r),
              people = as.character(n_people), without_adjustment = sprintf("%.2f", r0),
              with_adjustment = ifelse(is.na(lo1), sprintf("%.2f", r1),
                                       sprintf("%.2f (%.2f to %.2f)", r1, lo1, hi1))))
put(S27, "S27", "hypertension and diabetes added, on the same rows")
})


# #############################################################################
# Section 7.3. Agreement with amyloid PET by physiological factor (eTables 39, 40)
# Developed as 14_review_checks.R version 1.1.0, sections 6 and 6b
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
T7b  <- read_tab("T7b_thresholds")

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
# From section 5 of 14_review_checks.R: the cleared ratio and its categories.
FUJI_P217 <- "Fuji_plasma_ptau217"; FUJI_AB42 <- "Fuji_plasma_Ab42"
CUT_NEGATIVE <- 0.00370; CUT_POSITIVE <- 0.00738
cleared_category <- function(r) {
  r5 <- round(r, 5)
  dplyr::case_when(is.na(r5) ~ NA_character_, r5 <= CUT_NEGATIVE ~ "negative",
                   r5 >= CUT_POSITIVE ~ "positive", TRUE ~ "indeterminate")
}
# From section 4 of 14_review_checks.R: each participant's first draw, for the
# creatine kinase quarters (race and ethnicity, joined there, are not used here).
p_adni <- adni %>% arrange(RID, visit_no) %>% group_by(RID) %>% slice(1) %>% ungroup()

# =============================================================================
# 6. Agreement with amyloid PET per SD of each factor, and by band
# =============================================================================
# Each p-tau217 platform at its own threshold from step 10 (T7b, all samples)
# and the cleared ratio at its two cutoffs, against amyloid PET at a Centiloid
# of 20 and of 24. The main estimate is a linear probability model with
# standard errors clustered by participant, as in step 10: among PET-positive
# samples, a positive call on the factor; among PET-negative samples, a
# negative call; each adjusted for Centiloid as a natural spline (3 df), age
# and sex. The coefficient is the change in sensitivity or specificity, in
# percentage points, per SD of the factor at a fixed amyloid burden. One false
# discovery rate family per Centiloid cut. The thresholds were set in the
# whole sample and are held fixed, so the estimates are descriptive.
# Bands, descriptive: eGFR below 60 against 60 or above; body mass index 30 or
# above against below 30; creatine kinase in the top quarter of participants
# of the same sex against the rest. Intervals resample participants.
say("6. Agreement with amyloid PET per SD of each factor, and by band")
TESTS <- c(setNames(P217, short(P217)), cleared = "cleared")
test_data <- function(tn, cut) {
  d <- adni %>% filter(!is.na(centiloid))
  if (tn == "cleared") {
    d <- d %>% filter(!is.na(.data[[FUJI_P217]]), !is.na(.data[[FUJI_AB42]])) %>%
      mutate(cat = cleared_category(exp(.data[[FUJI_P217]] - .data[[FUJI_AB42]])),
             call_pos = cat == "positive", call_neg = cat == "negative")
    attr(d, "threshold") <- NA_real_
  } else {
    thr <- T7b$threshold_pg_ml[T7b$sample == "all samples" & T7b$centiloid_cut == cut &
                                 T7b$platform == tn]
    if (length(thr) != 1) stop("No threshold in T7b for ", tn, " at Centiloid ", cut)
    d <- d %>% filter(!is.na(.data[[TESTS[[tn]]]])) %>%
      mutate(call_pos = .data[[TESTS[[tn]]]] >= log(thr) - 1e-9, call_neg = !call_pos)
    attr(d, "threshold") <- thr
  }
  d$pet <- d$centiloid >= cut
  d$y_pos <- as.numeric(d$call_pos); d$y_neg <- as.numeric(d$call_neg)
  d
}
test_label <- function(tn) if (tn == "cleared") "FDA-cleared ratio, two cutoffs" else tn
R6 <- purrr::map_dfr(c(20, 24), function(cut) purrr::map_dfr(names(TESTS), function(tn) {
  d <- test_data(tn, cut)
  purrr::map_dfr(PRIMARY$var, function(v) {
    zv <- paste0("z_", v)
    rhs <- sprintf("%s + splines::ns(centiloid, df = 3) + age_at_visit + female", zv)
    sp <- fit_cl(stats::as.formula(paste("y_pos ~", rhs)), d[d$pet, , drop = FALSE], zv)
    sn <- fit_cl(stats::as.formula(paste("y_neg ~", rhs)), d[!d$pet, , drop = FALSE], zv)
    tibble(test = test_label(tn), threshold_pg_ml = attr(d, "threshold"), centiloid_cut = cut, var = v,
           measure = c("sensitivity", "specificity"),
           n_obs = c(sp$n_obs, sn$n_obs), n_people = c(sp$n_people, sn$n_people),
           pp_per_sd = 100 * c(sp$b, sn$b), lo = 100 * c(sp$lo, sn$lo), hi = 100 * c(sp$hi, sn$hi),
           p = c(sp$p, sn$p))
  })
})) %>%
  group_by(centiloid_cut) %>% mutate(q = p.adjust(p, "BH")) %>% ungroup()
put_rev(R6, "review_6_accuracy_per_sd")
print(as.data.frame(R6 %>% filter(centiloid_cut == 20) %>%
  transmute(test, var, measure, n_obs, per_sd = fmt(pp_per_sd, lo, hi), q = signif(q, 2))),
  row.names = FALSE)

ck_q75 <- p_adni %>% filter(!is.na(ck)) %>% group_by(female) %>%
  summarise(q75 = stats::quantile(ck, 0.75), .groups = "drop")
say("  creatine kinase top quarter: men %.2f U/L, women %.2f U/L",
    ck_q75$q75[ck_q75$female == 0], ck_q75$q75[ck_q75$female == 1])
ck_top <- function(d) {
  q <- ck_q75$q75[match(d$female, ck_q75$female)]
  ifelse(is.na(d$ck) | is.na(q), NA, d$ck >= q)
}
BANDS <- list(
  egfr = list(label = "eGFR", f = function(d) ifelse(is.na(d$egfr), NA, d$egfr < 60),
              yes = "below 60", no = "60 or above"),
  bmi  = list(label = "Body mass index", f = function(d) ifelse(is.na(d$bmi), NA, d$bmi >= 30),
              yes = "30 or above", no = "below 30"),
  ck   = list(label = "Creatine kinase", f = ck_top,
              yes = "top quarter within sex", no = "lower three quarters within sex"))
acc_stats <- function(d) {
  a <- d[d$band, , drop = FALSE]; b <- d[!d$band, , drop = FALSE]
  sens <- function(x) 100 * mean(x$call_pos[x$pet])
  spec <- function(x) 100 * mean(x$call_neg[!x$pet])
  ind  <- function(x) 100 * mean(!x$call_pos & !x$call_neg)
  c(sens_band = sens(a), sens_rest = sens(b), sens_difference = sens(a) - sens(b),
    spec_band = spec(a), spec_rest = spec(b), spec_difference = spec(a) - spec(b),
    indeterminate_band = ind(a), indeterminate_rest = ind(b))
}
R6b <- purrr::map_dfr(names(TESTS), function(tn) purrr::map_dfr(c(20, 24), function(cut) {
  d <- test_data(tn, cut)
  purrr::map_dfr(names(BANDS), function(bk) {
    db <- d; db$band <- BANDS[[bk]]$f(db); db <- db[!is.na(db$band), , drop = FALSE]
    boot_pct(db, acc_stats) %>%
      mutate(test = test_label(tn), threshold_pg_ml = attr(d, "threshold"), centiloid_cut = cut,
             band_factor = BANDS[[bk]]$label, band = BANDS[[bk]]$yes, rest = BANDS[[bk]]$no,
             samples_band = sum(db$band), people_band = n_distinct(db$RID[db$band]),
             pet_positive_band = sum(db$band & db$pet),
             samples_rest = sum(!db$band), pet_positive_rest = sum(!db$band & db$pet), .before = 1)
  })
}))
R6b <- R6b %>% filter(!(test != "FDA-cleared ratio, two cutoffs" & grepl("^indeterminate", stat)))
put_rev(R6b, "review_6b_accuracy_by_band")
})
