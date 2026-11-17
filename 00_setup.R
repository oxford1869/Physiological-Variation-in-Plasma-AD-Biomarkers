# =============================================================================
# 00_setup.R
# Paper 1: settings, paths and the estimation methods, in one place. Every
# other file in this folder sources this file first. Nothing here reads data.
#
# The pipeline is self-contained. It reads the raw ADNI and A4/LEARN releases
# itself. Both raw folders are read only and nothing is ever written to them.
#
# Everything the pipeline writes goes to one output folder, paper1_out/ beside
# this folder unless PAPER1_OUT_DIR names another.
# =============================================================================
PAPER1_VERSION <- "1.2.0"

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(lubridate)
  library(purrr)
  library(tibble)
})

# ---- Paths -------------------------------------------------------------------
# The launcher moves R to this folder, so the working directory is the code
# folder. Everything else is worked out from there, and every path can be
# overridden with an environment variable if the folders ever move.
if (!file.exists("00_setup.R"))
  stop("Run this through run_all.R, which moves R to the folder the code is in.")

CODE_DIR <- normalizePath(getwd())
OUT_ROOT <- Sys.getenv("PAPER1_OUT_ROOT", unset = dirname(CODE_DIR))
OUT_DIR  <- Sys.getenv("PAPER1_OUT_DIR", unset = file.path(OUT_ROOT, "paper1_out"))
TAB_DIR  <- file.path(OUT_DIR, "tables")
CHECK_DIR <- file.path(OUT_DIR, "checks")
for (p in c(OUT_DIR, TAB_DIR, CHECK_DIR)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

# The ADNI release is read from PAPER1_ADNI_DIR or, without it, from the folder
# ADNI beside the code folder (OUT_ROOT/ADNI). RECCMEDS.csv must be the full
# download (79,873 rows in the release used for the paper; an earlier download
# of 73,806 rows finds 388 statin records where the current file finds 733).
# More than one folder can be given, separated by semicolons, in the environment
# variable PAPER1_ADNI_DIR. Each file is then taken from the first folder that
# has it, and if the copies differ in size the larger one is used and the choice
# is printed, because an ADNI file only gains rows between releases.
ADNI_DIRS <- Sys.getenv("PAPER1_ADNI_DIR", unset = "")
ADNI_DIRS <- if (nzchar(ADNI_DIRS)) strsplit(ADNI_DIRS, ";", fixed = TRUE)[[1]] else
  file.path(OUT_ROOT, "ADNI")
ADNI_DIRS <- ADNI_DIRS[dir.exists(ADNI_DIRS)]
if (length(ADNI_DIRS) == 0)
  stop("The ADNI folder was not found. Expected ", file.path(OUT_ROOT, "ADNI"))
ADNI_DIR <- ADNI_DIRS[1]
A4_DIR <- Sys.getenv("PAPER1_A4_DIR", unset = file.path(OUT_ROOT, "A4&learn"))

ADNI_FILES <- list(
  plasma  = "FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.csv",
  merge   = "ADNIMERGE_DATA.csv",
  amyloid = "UCBERKELEY_AMY_6MM.csv",
  lab     = "LABDATA.csv",
  vitals  = "VITALS.csv",
  demog   = "PTDEMOG.csv",
  meds    = "RECCMEDS.csv",
  # the recorded medical history, for hypertension and diabetes (section 7.2)
  mhist   = "RECMHIST.csv"
)
# ADNI file names differ in case between releases, so resolve each one on disk.
adni_file <- function(key) {
  want <- ADNI_FILES[[key]]
  if (is.null(want)) stop("No ADNI file is registered under '", key, "'")
  pat <- paste0("^", gsub("\\.", "\\\\.", want), "$")
  hit <- unlist(lapply(ADNI_DIRS, function(d)
    list.files(d, pattern = pat, ignore.case = TRUE, full.names = TRUE)))
  if (length(hit) == 0) stop("File not found in any ADNI folder: ", want)
  if (length(hit) > 1) {
    sz <- file.info(hit)$size
    if (diff(range(sz)) > 0) {
      keep <- which.max(sz)
      say("  note: %s is in %d folders and the copies differ. Using the larger one,",
          want, length(hit))
      say("        %.1f MB in %s, rather than %.1f MB in %s.",
          sz[keep] / 1024^2, dirname(hit[keep]),
          sz[-keep][1] / 1024^2, dirname(hit[-keep][1]))
      hit <- hit[keep]
    }
  }
  hit[1]
}

# ---- Constants ---------------------------------------------------------------
YEAR_DAYS      <- 365.25
VISIT_GAP_DAYS <- 7      # plasma dates this close together are the same sample

# C2N p-tau217 below the limit of detection is reported as half that limit.
# The censored regression censors at the limit itself, not at the reported value.
C2N_LOD        <- 1.30
C2N_REPORTED   <- 0.65

# Equivalence margin for the between-platform difference, as a percentage per
# standard deviation of the factor. Tested against a 90% interval, so that the
# one-sided tests are at 5%.
EQUIV_MARGIN_PCT <- 5

# ADNI screening laboratory panel. RCT392 (creatinine) is handled separately
# because it is measured at every visit, not only at screening.
ADNI_LAB_CODES <- c(ck = "RCT14", ast = "RCT5", alt = "RCT4", urea = "RCT6",
                    hb = "HMT40", albumin = "RCT13", glucose = "RCT11")

# A4/LEARN uses the same test codes as ADNI, so the panel is selected by code
# rather than by the printed test name. Names differ between releases and a
# changed name fails silently; a code does not.
A4_LAB_CODES <- c(creatinine = "RCT392", ck = "RCT14", ast = "RCT5", alt = "RCT4",
                  urea = "RCT6", hb = "HMT40", albumin = "RCT13",
                  hba1c = "HMT370", crp = "SCT1528",
                  # A4/LEARN reports glucose under its own code: random glucose,
                  # RCT1669. ADNI's RCT11 does not occur in the A4 file.
                  glucose = "RCT1669")

# The A4/LEARN p-tau217 assay reports "<LLOQ" below 0.131 U/mL and ">ULOQ"
# above 1.43 U/mL, but the file also carries the unquantified reading itself in
# ORRESRAW. Those readings are used, so that no draw is dropped for being low,
# and the flags allow a sensitivity analysis that excludes them.
A4_PTAU217_LLOQ <- 0.131
A4_PTAU217_ULOQ <- 1.43

# Creatinine is reported in micromoles per litre in A4/LEARN and in milligrams
# per decilitre in ADNI. The eGFR equations take mg/dL.
UMOL_PER_MGDL <- 88.4

A4_FILES <- list(
  roche    = file.path("External Data", "biomarker_Plasma_Roche_Results.csv"),
  ptau217  = file.path("External Data", "biomarker_pTau217.csv"),
  lab      = file.path("External Data", "clrm_lab.csv"),
  subjinfo = file.path("Derived Data", "SUBJINFO.csv"),
  vitals   = file.path("Raw Data", "vitals.csv")
)
a4_file <- function(key) {
  want <- A4_FILES[[key]]
  if (is.null(want)) stop("No A4 file is registered under '", key, "'")
  p <- file.path(A4_DIR, want)
  if (!file.exists(p)) stop("File not found: ", p)
  p
}

# ---- The factors -------------------------------------------------------------
# One definition for both cohorts, so that a factor is z-scored and signed the
# same way everywhere. sign is the direction of interest: a positive estimate
# always means the direction named in the label. Scripts subset this table to
# the factors their cohort has.
FACTORS <- tibble::tribble(
  ~var,      ~label,                                  ~role,
  "egfr",    "eGFR (per SD lower)",                   "primary",
  "bmi",     "Body mass index (per SD higher)",       "primary",
  "l_ck",    "Creatine kinase (per SD higher, log)",  "primary",
  "hb",      "Haemoglobin (per SD lower)",            "primary",
  "urea",    "Urea (per SD higher)",                  "exploratory",
  "l_ast",   "AST (per SD higher, log)",              "exploratory",
  "l_alt",   "ALT (per SD higher, log)",              "exploratory",
  "hba1c",   "HbA1c (per SD higher)",                 "exploratory",
  "l_crp",   "C-reactive protein (per SD higher, log)", "exploratory",
  "albumin", "Albumin (per SD lower)",                "exploratory",
  "glucose", "Glucose (per SD higher)",               "exploratory"
)
FACTOR_SIGN <- c(egfr = -1, bmi = 1, l_ck = 1, hb = -1, urea = 1, l_ast = 1,
                 l_alt = 1, hba1c = 1, l_crp = 1, albumin = -1, glucose = 1)

# Covariates used by every level and difference model. A4/LEARN has no
# diagnosis, because everyone enrolled was cognitively unimpaired.
COVS_ADNI <- "age_at_visit + female + e4 + centiloid + dx + plasma_time"
COVS_A4   <- "age_at_visit + female + e4 + centiloid"
COVS      <- COVS_ADNI

# ---- Small utilities ---------------------------------------------------------
say <- function(...) cat(sprintf(...), "\n")

# A coefficient on the log scale as a percentage difference.
pct <- function(x) 100 * (exp(x) - 1)

zf  <- function(x) as.numeric(scale(x))

safe_log <- function(x) ifelse(!is.na(x) & x > 0, log(x), NA_real_)

save_p1 <- function(obj, name) {
  attr(obj, "paper1_version") <- PAPER1_VERSION
  attr(obj, "built_at") <- Sys.time()
  saveRDS(obj, file.path(OUT_DIR, paste0(name, ".rds")))
  invisible(obj)
}

# Read a CSV with every column as text, so nothing is coerced on the way in.
# Column names are upper-cased because the ADNI exports are inconsistent.
read_chr <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  x <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()),
                       show_col_types = FALSE, progress = FALSE)
  names(x) <- toupper(names(x))
  x
}

# Read a raw file and record where it came from. The record is written to
# checks/ by the build scripts, so that a set of results can always be traced
# back to the exact files it was built from. Two copies of an ADNI release are
# not always the same, and that has to be visible rather than silent.
.p1_sources <- new.env(parent = emptyenv())
.p1_sources$rows <- list()

read_source <- function(path, label) {
  x    <- read_chr(path)
  info <- file.info(path)
  .p1_sources$rows[[length(.p1_sources$rows) + 1]] <- tibble(
    label = label, file = basename(path), folder = dirname(path),
    rows = nrow(x), columns = ncol(x),
    modified = format(info$mtime, "%Y-%m-%d"),
    size_mb = round(info$size / 1024^2, 2))
  say("  %-9s %7d rows, %3d columns, modified %s", label, nrow(x), ncol(x),
      format(info$mtime, "%Y-%m-%d"))
  x
}

source_table <- function() bind_rows(.p1_sources$rows)

# The first of several candidate column names that the file actually has.
pick_col <- function(df, candidates, required = TRUE, label = candidates[1]) {
  hit <- candidates[toupper(candidates) %in% names(df)]
  if (length(hit) == 0) {
    if (required) stop(sprintf("None of [%s] found for %s. Columns are: %s",
                               paste(candidates, collapse = ", "), label,
                               paste(names(df), collapse = ", ")))
    return(NA_character_)
  }
  toupper(hit[1])
}

# The plasma file gives dates as day/month/year; the others as year-month-day.
parse_date_dmy <- function(x) as.Date(lubridate::parse_date_time(x, orders = c("dmy", "ymd", "mdy"), quiet = TRUE))
parse_date_ymd <- function(x) as.Date(lubridate::parse_date_time(x, orders = c("ymd", "dmy", "mdy"), quiet = TRUE))

# Screening and follow-up codes collapse to baseline for joins to ADNIMERGE.
recode_viscode <- function(v) {
  v <- as.character(v)
  dplyr::case_when(v %in% c("sc", "f", "scmri") ~ "bl", TRUE ~ v)
}

# eGFR from serum creatinine (mg/dL), age (years) and sex ("Female" / "Male").
# CKD-EPI 2021 is the equation used throughout; the 2009 one is kept for
# comparison with earlier work.
egfr_ckdepi_2021 <- function(scr, age, sex) {
  female <- sex == "Female"
  k      <- ifelse(female, 0.7, 0.9)
  alpha  <- ifelse(female, -0.241, -0.302)
  142 * pmin(scr / k, 1)^alpha * pmax(scr / k, 1)^-1.200 * 0.9938^age * ifelse(female, 1.012, 1)
}
egfr_ckdepi_2009 <- function(scr, age, sex) {
  female <- sex == "Female"
  k      <- ifelse(female, 0.7, 0.9)
  alpha  <- ifelse(female, -0.329, -0.411)
  141 * pmin(scr / k, 1)^alpha * pmax(scr / k, 1)^-1.209 * 0.993^age * ifelse(female, 1.018, 1)
}

# For each plasma visit (blood draw), find that person's record in target (an
# id and a date, for example an amyloid scan or a vitals visit) closest in time.
# The match runs from the draw to the record, as eMethods 7 describes: every
# draw takes its own nearest record, so two draws close together can share one
# record, and neither loses a record that lies inside its own window. (Version
# 1.0.0 matched the other way, each record to its nearest draw, which left a
# draw without a scan when a neighbouring draw was nearer to that scan.)
# window, in days, limits the match; draws with no record inside it are left
# out. Returns one row per matched draw: the record's columns with visit_no,
# visit_date and gap_days, the record date minus the plasma date, so negative
# means before the blood draw.
nearest_visit <- function(target, date_col, visits, id = "RID", window = Inf) {
  stopifnot(all(c(id, "visit_no", "visit_date") %in% names(visits)))
  target$.row <- seq_len(nrow(target))
  target$.d   <- target[[date_col]]
  matched <- visits %>%
    select(all_of(id), visit_no, visit_date) %>%
    inner_join(target %>% select(.row, all_of(id), .d), by = id,
               relationship = "many-to-many") %>%
    mutate(gap_days = as.integer(.d - visit_date)) %>%
    filter(abs(gap_days) <= window) %>%
    group_by(across(all_of(id)), visit_no) %>%
    slice_min(abs(gap_days), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(.row, visit_no, visit_date, gap_days)
  matched %>%
    left_join(target %>% select(-all_of(id)), by = ".row") %>%
    left_join(target %>% select(.row, all_of(id)), by = ".row") %>%
    select(-.row, -.d)
}

# ---- Estimation --------------------------------------------------------------

# Linear model with standard errors clustered by person (CR1). Used for every
# level and every between-platform difference, in both cohorts. level sets the
# interval: 0.95 for reporting, 0.90 for the equivalence test.
fit_cl <- function(formula, data, term, level = 0.95, cluster = "RID") {
  data <- as.data.frame(data)
  rownames(data) <- NULL
  mf  <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  idx <- as.integer(rownames(mf))
  if (length(idx) < 10) return(tibble(n_obs = length(idx), n_people = NA_integer_,
                                      b = NA_real_, se = NA_real_, lo = NA_real_,
                                      hi = NA_real_, p = NA_real_, level = level))
  cl <- data[[cluster]][idx]
  m  <- stats::lm(formula, data = data[idx, , drop = FALSE])
  if (!term %in% names(stats::coef(m)) || is.na(stats::coef(m)[[term]]))
    return(tibble(n_obs = length(idx), n_people = dplyr::n_distinct(cl), b = NA_real_,
                  se = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, level = level))
  # A covariate can be constant inside a subsample, for example sex in a
  # sex-stratified model or time at the first visit. lm sets those coefficients
  # to NA; they are dropped here so that the sandwich can be formed.
  cf   <- stats::coef(m)
  keep <- !is.na(cf)
  X <- stats::model.matrix(m)[, keep, drop = FALSE]
  e <- stats::residuals(m)
  G <- length(unique(cl)); n <- nrow(X); k <- ncol(X)
  bread <- solve(crossprod(X))
  meat  <- Reduce(`+`, lapply(split(seq_len(n), cl), function(i) {
    s <- crossprod(X[i, , drop = FALSE], e[i]); s %*% t(s)
  }))
  V  <- bread %*% meat %*% bread * (G / (G - 1)) * ((n - 1) / (n - k))
  b  <- stats::coef(m)[[term]]
  se <- sqrt(V[term, term])
  z  <- stats::qnorm(1 - (1 - level) / 2)
  tibble(n_obs = n, n_people = G, b = b, se = se,
         lo = b - z * se, hi = b + z * se, p = 2 * stats::pnorm(-abs(b / se)),
         level = level)
}

# Equivalence. Takes an estimate and its standard error on the log scale and
# asks whether the 90% interval lies inside a margin of margin_pct per standard
# deviation. Equivalent means the platforms agree to within the margin; not
# equivalent means only that agreement has not been shown, which is not the
# same as a difference.
equivalence <- function(b, se, margin_pct = EQUIV_MARGIN_PCT) {
  z  <- stats::qnorm(0.95)
  m  <- log(1 + margin_pct / 100)
  lo <- b - z * se
  hi <- b + z * se
  tibble(margin_pct = margin_pct,
         lo90 = pct(lo), hi90 = pct(hi),
         equivalent = !is.na(lo) & lo > -m & hi < m)
}

# Left-censored (Tobit) regression, for an outcome with a detection limit.
# y is the log value, observed above the limit and censored at log(limit).
# The limit itself is the censoring point, not the value the laboratory
# reports for a censored sample.
#
# censor_at is the limit on the original scale. It may be one number, or one
# number per row: the limit of a ratio moves with its denominator, so the
# %p-tau217 ratio is censored at a point that differs from sample to sample.
# censored, when given, names the rows that are censored, for a measure whose
# limit is known from a flag rather than from the value itself.
fit_tobit <- function(formula, data, term, censor_at, cluster = "RID", censored = NULL) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("The survival package is needed for the censored regression.")
  y_name <- all.vars(formula)[1]
  d <- as.data.frame(data)
  if (length(censor_at) > 1) stopifnot(length(censor_at) == nrow(d))
  d$.limit <- log(censor_at)
  if (!is.null(censored)) {
    stopifnot(length(censored) == nrow(d))
    d$.cens <- censored
  } else {
    d$.cens <- NA
  }
  ok <- !is.na(d[[y_name]]) & !is.na(d$.limit)
  d <- d[ok, , drop = FALSE]
  d$.obs <- if (all(is.na(d$.cens))) d[[y_name]] > d$.limit + 1e-8 else !d$.cens
  d$.y   <- ifelse(d$.obs, d[[y_name]], d$.limit)
  rhs <- paste(deparse(formula[[3]]), collapse = " ")
  f   <- stats::as.formula(sprintf("survival::Surv(.y, .obs, type = 'left') ~ %s + cluster(%s)",
                                   rhs, cluster))
  fit <- survival::survreg(f, data = d, dist = "gaussian")
  # survreg drops rows with a missing covariate, factor or cluster itself, so
  # the counts are taken on the rows the fit actually used, not on d. Version
  # 1.0.0 counted d, which overstated n_obs and n_censored whenever a covariate
  # was missing. The estimates are unaffected.
  used <- rep(TRUE, nrow(d))
  if (!is.null(fit$na.action)) used[as.integer(fit$na.action)] <- FALSE
  du <- d[used, , drop = FALSE]
  stopifnot(nrow(du) == NROW(fit$y))
  s   <- summary(fit)$table
  if (!term %in% rownames(s))
    return(tibble(n_obs = nrow(du), n_censored = sum(!du$.obs), b = NA_real_, se = NA_real_,
                  lo = NA_real_, hi = NA_real_, p = NA_real_))
  b <- s[term, 1]; se <- s[term, 2]
  tibble(n_obs = nrow(du), n_censored = sum(!du$.obs), n_people = dplyr::n_distinct(du[[cluster]]),
         b = b, se = se, lo = b - 1.96 * se, hi = b + 1.96 * se,
         p = 2 * stats::pnorm(-abs(b / se)))
}

# Slope of y on x when both are measured with error. lambda is the ratio of the
# error variance of y to that of x. With lambda = 1 this is the orthogonal
# (total least squares) slope; with the platform-specific error variances from
# Paper 2 it is the Deming slope.
deming_slope <- function(x, y, lambda = 1) {
  ok <- stats::complete.cases(x, y); x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(NA_real_)
  sxx <- stats::var(x); syy <- stats::var(y); sxy <- stats::cov(x, y)
  (syy - lambda * sxx + sqrt((syy - lambda * sxx)^2 + 4 * lambda * sxy^2)) / (2 * sxy)
}

# Reliability of a single measurement of a factor, and of a within-person
# deviation from that person's own mean, from a series of repeated values.
# Needs at least three values per person. Returns:
#   icc_single  between-person variance over total variance
#   rel_within  true within-person variation over observed within-person
#               variation, given k observations spanning a time variance of vt
reliability <- function(value, time, id, k = NULL, vt = NULL) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("The lme4 package is needed for the reliability calculation.")
  d <- tibble(x = as.numeric(value), t = as.numeric(time), id = id) %>%
    filter(!is.na(x), !is.na(t)) %>%
    group_by(id) %>% filter(n() >= 3) %>% ungroup()
  if (nrow(d) < 30 || dplyr::n_distinct(d$id) < 10)
    return(tibble(n_people = dplyr::n_distinct(d$id), icc_single = NA_real_, rel_within = NA_real_))
  m <- tryCatch(lme4::lmer(x ~ t + (t | id), data = d,
                           control = lme4::lmerControl(calc.derivs = FALSE)),
                error = function(e) NULL)
  if (is.null(m)) return(tibble(n_people = dplyr::n_distinct(d$id),
                                icc_single = NA_real_, rel_within = NA_real_))
  vc   <- as.data.frame(lme4::VarCorr(m))
  s2i  <- vc$vcov[1]                 # between-person variance of the level
  s2s  <- vc$vcov[2]                 # between-person variance of the rate of change
  s2e  <- vc$vcov[nrow(vc)]          # occasion to occasion variance
  own  <- d %>% group_by(id) %>% summarise(vt = mean((t - mean(t))^2), k = n(), .groups = "drop")
  if (is.null(vt)) vt <- mean(own$vt)
  if (is.null(k))  k  <- mean(own$k)
  true_w <- s2s * vt
  err_w  <- s2e * (k - 1) / k
  tibble(n_people = dplyr::n_distinct(d$id),
         icc_single = s2i / (s2i + s2e),
         rel_within = true_w / (true_w + err_w))
}

# Regression calibration. An estimate for a factor measured with error is
# attenuated by the reliability of that measurement, so dividing by it gives
# the estimate that a perfectly measured factor would have produced.
rc_correct <- function(b, reliability) b / reliability

# The eGFR spline for C2N p-tau217 by left-censored regression (version 1.2.0).
# Until 1.2.0 the spline for C2N p-tau217 was fitted by least squares on the
# samples above the detection limit, which selects on the outcome, the handling
# eMethods 5 rejects. survreg accepts spline columns like any other covariate,
# so the curve is now fitted with every sample, a value below the limit
# contributing the information that it lies below log(1.30).
#
# The spline is parameterised as in step 03: the first k - 1 natural spline
# columns and eGFR itself, which span the same functions as the k spline
# columns. The test against a straight line is then the Wald test of the k - 1
# spline columns, with the cluster-robust variance survreg gives when a cluster
# term is present. That the two parameterisations fit identically is checked
# through the log-likelihood. Each point on the curve is the difference from an
# eGFR of reference, so the covariates cancel and only the eGFR columns enter.
#   d        rows with the measure on the log scale and c2n_at_floor
#   at       eGFR values at which the difference from reference is returned
# Returns a list: points (egfr, pct, lo, hi), p_curved, n_obs, n_people,
# n_censored.
spline_censored <- function(d, measure, covs, id, at, knots = c(45, 60, 75),
                            reference = 90) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("The survival package is needed for the censored spline.")
  d <- as.data.frame(d)
  d <- d[!is.na(d[[measure]]) & !is.na(d$egfr) & !is.na(d$centiloid), , drop = FALSE]
  rownames(d) <- NULL
  # The basis is built on these rows, as in step 03 and step 11, so its
  # boundary knots are the same.
  B <- splines::ns(d$egfr, knots = knots, Boundary.knots = range(d$egfr))
  k <- ncol(B)
  sp_all <- paste0(".sp", seq_len(k))
  for (i in seq_len(k)) d[[sp_all[i]]] <- B[, i]
  d$.lin <- d$egfr
  sp <- sp_all[-k]
  rhs  <- paste(c(sp, ".lin", covs), collapse = " + ")
  rhs0 <- paste(c(sp_all, covs), collapse = " + ")
  vars <- unique(c(measure, all.vars(stats::as.formula(paste("~", rhs))), id, "c2n_at_floor"))
  d <- d[stats::complete.cases(d[, vars, drop = FALSE]), , drop = FALSE]
  rownames(d) <- NULL
  d$.obs <- !as.logical(d$c2n_at_floor)
  d$.y   <- ifelse(d$.obs, d[[measure]], log(C2N_LOD))
  fit_one <- function(r) survival::survreg(
    stats::as.formula(sprintf("survival::Surv(.y, .obs, type = 'left') ~ %s + cluster(%s)", r, id)),
    data = d, dist = "gaussian")
  fit  <- fit_one(rhs)
  fit0 <- fit_one(rhs0)
  ll <- c(fit$loglik[2], fit0$loglik[2])
  if (abs(ll[1] - ll[2]) > 1e-6 * abs(ll[2]))
    stop("The reparameterised censored spline does not reproduce the spline fit for ", measure)
  b <- stats::coef(fit)
  if (anyNA(b[c(sp, ".lin")])) stop("A spline coefficient is not estimable for ", measure)
  V <- stats::vcov(fit)[names(b), names(b)]
  W <- as.numeric(t(b[sp]) %*% solve(V[sp, sp]) %*% b[sp])
  B_ref <- as.numeric(stats::predict(B, newx = reference))
  pts <- purrr::map_dfr(at, function(x) {
    cv <- stats::setNames(rep(0, length(b)), names(b))
    cv[sp] <- (as.numeric(stats::predict(B, newx = x)) - B_ref)[-k]
    cv[".lin"] <- x - reference
    e <- sum(cv * b)
    s <- sqrt(as.numeric(t(cv) %*% V %*% cv))
    tibble(egfr = x, pct = pct(e), lo = pct(e - 1.96 * s), hi = pct(e + 1.96 * s))
  })
  list(points = pts, p_curved = stats::pchisq(W, df = length(sp), lower.tail = FALSE),
       n_obs = nrow(d), n_people = dplyr::n_distinct(d[[id]]), n_censored = sum(!d$.obs))
}

say("00_setup.R loaded, Paper 1 version %s", PAPER1_VERSION)
say("  code    %s", CODE_DIR)
say("  ADNI    %s", paste(ADNI_DIRS, collapse = "\n          "))
say("  A4      %s", A4_DIR)
say("  output  %s", OUT_DIR)
