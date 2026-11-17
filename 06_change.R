# =============================================================================
# 06_change.R
# Stage 6. Research question 4: is the rate of change in a factor related to the
# rate of change in p-tau217? A4 only.
#   6.1  the hybrid model, the reliability of change, the joint growth model and
#        its bootstrap, the span rules (Table 5, eTables 18, 19, 27)
#   6.2  the joint growth model on a clinical scale and by trial arm (eTable 19)
# Writes paper1_out/tables/T5*.csv, R2_change_reliability.csv, review/review_3*.csv.
# =============================================================================

# #############################################################################
# Section 6.1. Change on change
# Developed as 06_rq4_change.R version 1.4.1, whole
# One change: the span sensitivity (eTable 27) is fitted for eGFR only.
# #############################################################################
local({
# =============================================================================
# 06_rq4_change.R
# Research question 4: is the rate of change in a factor related to the rate of
# change in the biomarker?
#
# Four published papers ask for this analysis and none has done it. The obvious
# way to do it fails, and the reason it fails is itself a result worth
# reporting, so the script does both: the analysis that fails, with the
# measurement that explains why, and then the analysis that works.
#
#   1  The hybrid model. Each factor is split into a person's own average and
#      that person's departure from it at each draw. The within-person estimate
#      is far smaller than the between-person one, which looks like a finding
#      and is not.
#   2  The reliability of a within-person change, from the repeated laboratory
#      values. Most of what looks like change in a factor is noise, which is
#      what makes the within-person estimate in step 1 too small.
#   3  The joint growth model. Both series are modelled together, each person
#      having a true level and a true rate of change in both, and the quantity
#      of interest is the correlation between the two rates of change. Using
#      every repeated value rather than the concurrent pairs is what makes the
#      noise average out.
#
# Intervals for the correlation come from a bootstrap that resamples people and
# refits, because the model's own standard errors are not always usable. The
# bootstrap is off by default because it is slow. To run it:
#      Sys.setenv(PAPER1_BOOT = 500)
# before sourcing the launcher. PAPER1_BOOT_VARS chooses which factors, and
# PAPER1_BOOT_SEED fixes the draws so the interval can be reproduced. Factors
# bootstrapped earlier are kept, so one factor can be rerun with more draws.
# Each draw refits the model in about two seconds, so 500 draws take about 16
# minutes per group and about 50 minutes for the whole sample and the two trial
# arms together. A draw counts as usable only if its fit converged; draws that
# fail and draws that stop without converging are counted separately.
#
# Reads    paper1_out/a4_long.rds, a4_lab.rds, tables/R1_reliability.csv
# Writes   paper1_out/tables/
#            T5_hybrid_within_between.csv  the analysis that fails
#            R2_change_reliability.csv     why it fails
#            T5b_joint_growth.csv          the analysis that works
#            T5c_joint_bootstrap.csv       intervals, when the bootstrap is run
#            T5d_span_sensitivity.csv      the joint model under two other
#                                          eligibility rules and with the
#                                          p-tau217 draws outside the assay's
#                                          quantification range excluded
#                                          (formatted as S25 by step 10)
#
# Changes in 1.2.0
#   the joint model is refitted with the p-tau217 draws outside the
#     quantification range excluded, the check eMethods 5 states, as a fourth
#     block of T5d
#
# Changes in 1.3.0
#   the joint model is refitted with a quadratic term for time in each series
#     (fixed effects only, so the estimand is still the correlation between the
#     true linear rates of change at the centre of follow-up), the check of
#     non-linear trajectories in the analysis plan, as a fifth block of T5d
#
# Changes in 1.4.0
#   the joint model adjusts both series, level and rate of change, for age at
#     the first draw, sex, APOE4 and Centiloid, the covariates of eMethods 7.
#     Without them part of the slope correlation could come from age, which
#     speeds both the fall in kidney function and the rise in p-tau217. The
#     estimand is now the correlation between the rates of change that these
#     covariates do not explain. People missing any of them are left out.
#
# Changes in 1.4.1
#   the hybrid model and the reliability of change use the trial participants
#     only, as eMethods 3 states that the change analysis is confined to A4;
#     before 1.4.1 LEARN participants with two draws could enter them
#
# The span sensitivity has its own eGFR bootstrap, which runs by default with
# 1,000 draws (about 25 minutes for each rule). PAPER1_N_BOOT sets the number
# of draws, and 0 leaves the interval out.
#
# Needs    lme4 and nlme. Takes about four minutes without either bootstrap
#          (PAPER1_BOOT unset and PAPER1_N_BOOT = 0),
#          because every model is fitted twice, with two different optimisers,
#          to check that the answer is a real optimum.
# =============================================================================
SCRIPT_VERSION <- "1.4.1"
cat("06_rq4_change.R version", SCRIPT_VERSION, "\n")

if (!exists("PAPER1_VERSION")) source("00_setup.R")
for (pkg in c("lme4", "nlme")) if (!requireNamespace(pkg, quietly = TRUE))
  stop("The ", pkg, " package is needed for question 4.")

L   <- readRDS(file.path(OUT_DIR, "a4_long.rds"))
lab <- readRDS(file.path(OUT_DIR, "a4_lab.rds"))
BOOT <- suppressWarnings(as.integer(Sys.getenv("PAPER1_BOOT", "0")))
if (is.na(BOOT)) BOOT <- 0

FAC <- FACTORS %>% filter(var %in% names(L), var %in% names(lab))
say("Draws: %d in %d people. Factors available: %s",
    nrow(L), n_distinct(L$BID), paste(FAC$var, collapse = ", "))

# ---- 1. The hybrid model, which is the analysis that fails -------------------
# Between is the difference between people; within is the same person at a
# draw when the factor is higher than their own average. If a change in the
# factor moved the biomarker, the two would agree.
# The change analysis is confined to A4 (eMethods 3). The joint model enforces
# that through its span rule; the hybrid model and the reliability of change
# enforce it here.
L_A4 <- L %>% filter(arm != "LEARN")
hybrid <- function(v, data = L_A4) {
  d <- data %>% filter(!is.na(.data[[v]]), !is.na(pTau217)) %>%
    group_by(BID) %>% filter(n() >= 2) %>% ungroup()
  if (nrow(d) < 100) return(NULL)
  s <- sd(d[[v]])
  sign <- FACTOR_SIGN[[v]]
  d <- d %>% group_by(BID) %>%
    mutate(between = mean(.data[[v]]), within = .data[[v]] - between) %>%
    ungroup() %>%
    mutate(between = sign * (between - mean(between)) / s, within = sign * within / s)
  arms <- length(unique(d$arm)) > 1
  f <- as.formula(paste("pTau217 ~ between + within + plasma_time + age0 + female + e4 +",
                        "centiloid +", if (arms) "arm +" else "", "(1 | BID)"))
  m <- lme4::lmer(f, data = d, REML = TRUE)
  co <- summary(m)$coefficients
  V  <- as.matrix(vcov(m))
  dif <- co["within", 1] - co["between", 1]
  sed <- sqrt(V["within", "within"] + V["between", "between"] - 2 * V["within", "between"])
  tibble(var = v, n_people = n_distinct(d$BID), n_draws = nrow(d),
         between_pct = pct(co["between", 1]),
         between_lo = pct(co["between", 1] - 1.96 * co["between", 2]),
         between_hi = pct(co["between", 1] + 1.96 * co["between", 2]),
         within_pct = pct(co["within", 1]),
         within_lo = pct(co["within", 1] - 1.96 * co["within", 2]),
         within_hi = pct(co["within", 1] + 1.96 * co["within", 2]),
         p_difference = 2 * pnorm(-abs(dif / sed)),
         within_b = co["within", 1])
}
T5 <- purrr::map_dfr(FAC$var, hybrid)

# ---- 2. Why it fails ---------------------------------------------------------
# The reliability of a within-person change depends on how many values a person
# has and over how long, so it is estimated for this sample rather than in
# general. Below about 0.20 there is almost nothing left to detect.
say("Estimating how much of a within-person change is real")
R2 <- purrr::map_dfr(FAC$var, function(v) {
  d <- L_A4 %>% filter(!is.na(.data[[v]]), !is.na(pTau217)) %>%
    group_by(BID) %>% filter(n() >= 2) %>% ungroup()
  own <- d %>% group_by(BID) %>%
    summarise(vt = mean((plasma_time - mean(plasma_time))^2), k = n(),
              span = max(plasma_time) - min(plasma_time), .groups = "drop")
  # span_years is the mean of last minus first plasma draw. It used to be
  # 2 x sqrt(vt), which equals the span only for a person with two draws and
  # understates it for the majority who have more. The quantity the
  # reliability actually uses, the mean time variance of the draws, is kept
  # beside it so that the calculation can be repeated from the table.
  reliability(lab[[v]], lab$lab_time, lab$BID,
              k = mean(own$k), vt = mean(own$vt)) %>%
    mutate(var = v, span_years = mean(own$span), time_variance = mean(own$vt),
           draws_per_person = mean(own$k), .before = 1)
})
# Correcting an estimate for measurement error means dividing by the
# reliability, which is only sensible while the reliability is not tiny. Below
# 0.20 the correction multiplies noise by five or more, so it is left out.
T5 <- T5 %>%
  left_join(R2 %>% select(var, rel_within), by = "var") %>%
  mutate(within_corrected_pct = if_else(rel_within >= 0.20,
                                        pct(rc_correct(within_b, rel_within)), NA_real_))
R2 <- R2 %>% left_join(FACTORS %>% select(var, label), by = "var") %>%
  select(var, label, n_people, draws_per_person, span_years, time_variance,
         icc_single, rel_within)
write_csv(T5, file.path(TAB_DIR, "T5_hybrid_within_between.csv"))
write_csv(R2, file.path(TAB_DIR, "R2_change_reliability.csv"))

say("The within-person estimate against the between-person estimate:")
print(as.data.frame(T5 %>% transmute(
  factor = var, people = n_people, draws = n_draws,
  between = sprintf("%.2f (%.2f to %.2f)", between_pct, between_lo, between_hi),
  within = sprintf("%.2f (%.2f to %.2f)", within_pct, within_lo, within_hi),
  reliability = round(rel_within, 2),
  within_corrected = round(within_corrected_pct, 2))), row.names = FALSE)

# ---- 3. The joint growth model, which is the analysis that works -------------
# Both series for one person are modelled together. Each person has a true
# level and a true rate of change in each, and the four are allowed to
# correlate. The quantity of interest is the correlation between the two rates
# of change: does p-tau217 rise faster in the people whose kidney function
# falls faster? Each series is standardised so the two are on the same scale,
# and time is centred at 2.5 years so the level and the rate are not confounded.
TIME_CENTRE <- 2.5

# A person can only contribute to a change-on-change analysis if both series
# are measured more than once, over a long enough span for a rate of change to
# mean anything. LEARN participants have a median of one creatinine value, so
# this rule excludes that cohort, which is as it should be: without repeats
# the model would be estimating their rate of change from nothing.
MIN_VALUES <- 2
MIN_SPAN_YEARS <- 1

# The main analysis applies the span rule to the laboratory series only, as
# eMethods 3 describes. The two span sensitivities in section 5 change that
# rule, which is why the span required of each series is an argument: lab_span
# for the laboratory series and bio_span for the p-tau217 series.
# quantified_only keeps only the p-tau217 draws inside the assay's
# quantification range (0.131 to 1.43 U/mL), for the check that using the
# instrument readings outside it does not produce the result. Eligibility is
# then judged on the draws that remain.
bio_draws <- function(quantified_only = FALSE)
  if (quantified_only) L %>% filter(quantified) else L
eligible_people <- function(v, lab_span = MIN_SPAN_YEARS, bio_span = 0, quantified_only = FALSE) {
  f <- lab %>% filter(!is.na(.data[[v]])) %>%
    group_by(BID) %>%
    summarise(n = n(), span = max(lab_time) - min(lab_time), .groups = "drop") %>%
    filter(n >= MIN_VALUES, span >= lab_span)
  b <- bio_draws(quantified_only) %>% filter(!is.na(pTau217)) %>% group_by(BID) %>%
    summarise(n = n(), span = max(plasma_time) - min(plasma_time), .groups = "drop") %>%
    filter(n >= MIN_VALUES, span >= bio_span)
  intersect(f$BID, b$BID)
}

# Person-level covariates of the joint model (eMethods 7): age at the first
# draw and Centiloid standardised over the people in the model, sex and APOE4
# as 0/1. Each enters the level and the rate of change of both series.
JOINT_COVS <- c("age_z", "female", "e4", "cl_z")
person_covs <- L %>% distinct(BID, age0, female, e4, centiloid) %>%
  group_by(BID) %>% slice(1) %>% ungroup()
joint_fixed <- function(quadratic = FALSE) {
  base <- c("ix", "ip", "ix:t", "ip:t", if (quadratic) c("ix:I(t^2)", "ip:I(t^2)"))
  adj  <- as.vector(outer(c("ix", "ip", "ix:t", "ip:t"), JOINT_COVS, paste, sep = ":"))
  as.formula(paste("y ~ 0 +", paste(c(base, adj), collapse = " + ")))
}

joint_data <- function(v, people = NULL, lab_span = MIN_SPAN_YEARS, bio_span = 0,
                       quantified_only = FALSE) {
  keep <- eligible_people(v, lab_span, bio_span, quantified_only)
  if (!is.null(people)) keep <- intersect(keep, people)
  cv <- person_covs %>% filter(BID %in% keep) %>%
    filter(!is.na(age0), !is.na(female), !is.na(e4), !is.na(centiloid)) %>%
    mutate(age_z = as.numeric(scale(age0)), cl_z = as.numeric(scale(centiloid))) %>%
    select(BID, all_of(JOINT_COVS))
  keep <- intersect(keep, cv$BID)
  bio <- bio_draws(quantified_only) %>% filter(!is.na(pTau217), BID %in% keep) %>%
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

# The fit warns about singular precision matrices while it searches, because a
# person with only a few values gives a singular block. Those warnings are
# transient, so rather than reporting them the script checks the fit the way
# that matters: it refits with a second optimiser and compares. Two optimisers
# reaching the same likelihood and the same correlation is evidence that the
# answer is a real optimum and not an artefact of the search.
#
# With returnObject = TRUE, nlme returns a fit even when the optimiser stops
# without converging, and says so only in a warning ("convergence error code"
# from optim or nlminb, or "reached without convergence" when maxIter runs
# out). Those warnings are the only record of it, so they are caught here
# rather than suppressed, and the fit carries attr(, "converged"). A bootstrap
# draw is usable only if its fit converged, as eMethods 8.10 states; before
# version 1.1.0 a draw that stopped at the iteration limit was counted too.
NONCONVERGENCE <- "convergence error code|without convergence"
# quadratic = TRUE adds a fixed quadratic term for time in each series. The
# random effects are unchanged, so the slope correlation keeps its meaning: the
# correlation between the true rates of change at the centre of follow-up.
fit_joint <- function(d, optimiser = "optim", quadratic = FALSE) {
  fixed <- joint_fixed(quadratic)
  ctrl <- nlme::lmeControl(opt = optimiser, maxIter = 200, msMaxIter = 200,
                           returnObject = TRUE, apVar = FALSE)
  warned <- character()
  m <- withCallingHandlers(tryCatch(
    nlme::lme(fixed,
              random = list(BID = nlme::pdSymm(~ 0 + ix + ip + ix:t + ip:t)),
              weights = nlme::varIdent(form = ~ 1 | o),
              data = d, control = ctrl, method = "REML"),
    error = function(e) NULL),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (!is.null(m)) attr(m, "converged") <- !any(grepl(NONCONVERGENCE, warned))
  m
}
converged <- function(m) !is.null(m) && isTRUE(attr(m, "converged"))

slope_correlation <- function(m) {
  if (is.null(m)) return(NA_real_)
  vc <- as.matrix(nlme::getVarCov(m))
  nm <- colnames(vc)
  i <- grep("ix:t|t:ix", nm)[1]
  j <- grep("ip:t|t:ip", nm)[1]
  if (is.na(i) || is.na(j)) return(NA_real_)
  vc[i, j] / sqrt(vc[i, i] * vc[j, j])
}

say("Fitting the joint growth model for each factor")
groups <- list("everyone" = NULL)
if ("arm" %in% names(L)) for (a in sort(unique(L$arm)))
  groups[[a]] <- L$BID[L$arm == a]

T5b <- purrr::map_dfr(FAC$var, function(v) purrr::map_dfr(names(groups), function(g) {
  d <- joint_data(v, groups[[g]])
  np <- n_distinct(d$BID)
  if (np < 100) return(NULL)
  m1 <- fit_joint(d, "optim")
  m2 <- fit_joint(d, "nlminb")
  r1 <- slope_correlation(m1)
  r2 <- slope_correlation(m2)
  tibble(var = v, group = g, n_people = np,
         n_factor_values = sum(d$o == "fac"), n_draws = sum(d$o == "bio"),
         slope_correlation = r1, correlation_other_optimiser = r2,
         converged = converged(m1), converged_other_optimiser = converged(m2),
         agree = !is.na(r1) && !is.na(r2) && abs(r1 - r2) < 0.05,
         logLik = if (is.null(m1)) NA_real_ else as.numeric(logLik(m1)))
}))
write_csv(T5b, file.path(TAB_DIR, "T5b_joint_growth.csv"))

say("Each factor is oriented in the direction named in its label, so a positive")
say("correlation means p-tau217 rises faster in the people whose factor moves")
say("further in that direction. For eGFR that is faster decline in kidney function.")
print(as.data.frame(T5b %>% filter(group == "everyone") %>%
  transmute(factor = var, people = n_people, values = n_factor_values, draws = n_draws,
            correlation = round(slope_correlation, 2), converged, agree)), row.names = FALSE)
print(as.data.frame(T5b %>% filter(var == "egfr") %>%
  transmute(group, people = n_people, correlation = round(slope_correlation, 2),
            agree)), row.names = FALSE)
excl <- L %>% distinct(BID, arm) %>% filter(!BID %in% eligible_people("egfr")) %>% count(arm)
if (nrow(excl) > 0) {
  say("People with too few repeated creatinine values to contribute:")
  print(as.data.frame(excl), row.names = FALSE)
}

# ---- 4. Intervals by bootstrap ----------------------------------------------
# People are resampled with replacement and the model refitted, because the
# model's own standard errors for a variance are not always usable. A method of
# moments shortcut was tried and rejected: it can return a negative variance.
BOOT_VARS <- strsplit(Sys.getenv("PAPER1_BOOT_VARS", "egfr"), ",")[[1]]
# The interval is a percentile of resampled values, so it moves a little from
# one run to the next. A published interval has to be reproducible, so each
# factor and group is seeded separately: rerunning one factor gives the same
# answer as before and does not disturb the others.
BOOT_SEED <- suppressWarnings(as.integer(Sys.getenv("PAPER1_BOOT_SEED", "1")))
if (is.na(BOOT_SEED)) BOOT_SEED <- 1
# Resampling replaces whole people, and each drawn person gets a new
# identifier so that a person drawn twice counts as two people. A draw is
# usable only if its fit converged: one that errors and one that stops without
# converging are both counted, separately, and left out of the interval.
boot_one <- function(d, people, quadratic = FALSE) {
  ids <- sample(people, length(people), replace = TRUE)
  dd <- tibble(BID = ids, .new = seq_along(ids)) %>%
    inner_join(d, by = "BID", relationship = "many-to-many") %>%
    mutate(BID = .new) %>% select(-.new)
  m <- fit_joint(dd, quadratic = quadratic)
  if (is.null(m)) return(c(r = NA_real_, status = 1))
  if (!converged(m)) return(c(r = NA_real_, status = 2))
  r <- slope_correlation(m)
  c(r = r, status = if (is.na(r)) 1 else 0)
}
# The draws for one data set, seeded so that the interval can be reproduced.
# label is used only in the progress messages.
boot_draws <- function(d, B, label, quadratic = FALSE) {
  people <- unique(d$BID)
  set.seed(BOOT_SEED)
  t0 <- Sys.time()
  res <- vapply(seq_len(B), function(i) {
    if (i %% 50 == 0) say("  %s: %d of %d, %.0f minutes so far", label, i, B,
                          as.numeric(difftime(Sys.time(), t0, units = "mins")))
    tryCatch(boot_one(d, people, quadratic), error = function(e) c(r = NA_real_, status = 1))
  }, numeric(2))
  rs <- res["r", ]
  # The median of the resampled correlations is reported beside the
  # estimate itself. Resampling people with replacement can shift a
  # variance component, and a median far from the estimate is the warning
  # that the interval should not be read as symmetric around it.
  tibble(draws = B, usable = sum(res["status", ] == 0),
         failed = sum(res["status", ] == 1),
         not_converged = sum(res["status", ] == 2),
         bootstrap_median = unname(median(rs, na.rm = TRUE)),
         lo = unname(quantile(rs, 0.025, na.rm = TRUE)),
         hi = unname(quantile(rs, 0.975, na.rm = TRUE)),
         seed = BOOT_SEED,
         minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
}

if (BOOT > 0) {
  say("Bootstrapping the correlation, %d draws for %s. This is the slow part.",
      BOOT, paste(BOOT_VARS, collapse = ", "))
  T5c <- purrr::map_dfr(intersect(BOOT_VARS, FAC$var), function(v)
    purrr::map_dfr(names(groups), function(g) {
      d <- joint_data(v, groups[[g]])
      if (n_distinct(d$BID) < 100) return(NULL)
      bd <- boot_draws(d, BOOT, paste(v, g, sep = ", "))
      tibble(var = v, group = g, draws = bd$draws, usable = bd$usable,
             failed = bd$failed, not_converged = bd$not_converged,
             correlation = slope_correlation(fit_joint(d))) %>%
        bind_cols(bd %>% select(bootstrap_median, lo, hi, seed, minutes))
    }))
  # Factors bootstrapped in an earlier run are kept, so a single factor can be
  # rerun with more draws without losing the others.
  boot_file <- file.path(TAB_DIR, "T5c_joint_bootstrap.csv")
  if (file.exists(boot_file)) {
    old <- read_csv(boot_file, show_col_types = FALSE)
    if (all(c("var", "group") %in% names(old)))
      T5c <- bind_rows(anti_join(old, T5c, by = c("var", "group")), T5c) %>%
        arrange(match(var, FAC$var), match(group, names(groups)))
  }
  write_csv(T5c, boot_file)
  print(as.data.frame(T5c %>% mutate(across(c(correlation, bootstrap_median, lo, hi),
                                                    ~ round(.x, 3)))), row.names = FALSE)
  # A row kept from a run before version 1.1.0 has no convergence columns, and
  # its usable count means only that the fit did not fail. Say so rather than
  # let it pass as a count of converged draws.
  if ("not_converged" %in% names(T5c) && any(is.na(T5c$not_converged)))
    say("  note: rows without not_converged come from an earlier version, whose usable count did not check convergence.")
} else {
  say("The bootstrap was not run, so the correlations have no intervals yet.")
  say("  To run it:  Sys.setenv(PAPER1_BOOT = 500)  before sourcing the launcher.")
  say("  It takes roughly an hour. Add PAPER1_BOOT_VARS to bootstrap more than eGFR.")
}

# ---- 5. The span rule, two sensitivity analyses ------------------------------
# The main analysis asks for at least two laboratory values spanning at least a
# year, and at least two p-tau217 draws at any span (eMethods 3). The span
# rule is a judgement, so the joint model is refitted under the two obvious
# alternatives: no span rule at all, and the one-year rule on both series.
# The first shows whether the rule matters; the second is the stricter reading
# of "a long enough span for a rate of change to mean anything".
#
# The eGFR correlation gets an interval by the same bootstrap as section 4,
# with the same seed. The number of draws is N_BOOT, 1,000 unless
# PAPER1_N_BOOT says otherwise; 0 leaves the intervals out. The main rule's
# interval is not recomputed here: it is taken from T5c, so it is the one
# the paper reports.
N_BOOT <- suppressWarnings(as.integer(Sys.getenv("PAPER1_N_BOOT", "1000")))
if (is.na(N_BOOT) || N_BOOT < 0) N_BOOT <- 1000
SPAN_RULES <- tibble::tribble(
  ~rule,        ~lab_span,      ~bio_span,      ~quantified_only, ~quadratic,
  "main",       MIN_SPAN_YEARS, 0,              FALSE,            FALSE,
  "no_span",    0,              0,              FALSE,            FALSE,
  "span_both",  MIN_SPAN_YEARS, MIN_SPAN_YEARS, FALSE,            FALSE,
  "quantified", MIN_SPAN_YEARS, 0,              TRUE,             FALSE,
  "quadratic",  MIN_SPAN_YEARS, 0,              FALSE,            TRUE)
SPAN_DESCRIPTION <- c(
  main       = "one-year span on the laboratory series only (main analysis)",
  no_span    = "no span rule on either series",
  span_both  = "one-year span on both series",
  quantified = "main rule, p-tau217 draws outside the quantification range excluded",
  quadratic  = "main rule, quadratic term for time in each series")
if (!"quantified" %in% names(L)) {
  say("  a4_long.rds has no quantified column, so the out-of-range check is skipped; rerun step 02.")
  SPAN_RULES <- SPAN_RULES %>% filter(!quantified_only)
}
# eGFR only: the creatine kinase and haemoglobin rows were removed from eTable 27
# on 16 Nov 2026 (decisions log 8.3).
SPAN_VARS <- intersect("egfr", FAC$var)
SPAN_BOOT_VAR <- "egfr"

say("Refitting the joint model under two other span rules")
boot_file <- file.path(TAB_DIR, "T5c_joint_bootstrap.csv")
T5c_now <- if (file.exists(boot_file)) read_csv(boot_file, show_col_types = FALSE) else NULL
T5d <- purrr::map_dfr(seq_len(nrow(SPAN_RULES)), function(k) {
  ru <- SPAN_RULES[k, ]
  purrr::map_dfr(SPAN_VARS, function(v) {
    d <- joint_data(v, NULL, ru$lab_span, ru$bio_span, ru$quantified_only)
    m <- fit_joint(d, quadratic = ru$quadratic)
    out <- tibble(rule = ru$rule, description = SPAN_DESCRIPTION[[ru$rule]],
                  lab_span_years = ru$lab_span, ptau217_span_years = ru$bio_span,
                  var = v, n_people = n_distinct(d$BID),
                  n_factor_values = sum(d$o == "fac"), n_draws = sum(d$o == "bio"),
                  slope_correlation = slope_correlation(m), converged = converged(m),
                  draws = NA_integer_, usable = NA_integer_, not_converged = NA_integer_,
                  lo = NA_real_, hi = NA_real_)
    if (v != SPAN_BOOT_VAR) return(out)
    if (ru$rule == "main") {
      b <- if (!is.null(T5c_now)) T5c_now %>% filter(var == v, group == "everyone") else NULL
      if (!is.null(b) && nrow(b) == 1)
        out <- out %>% mutate(draws = as.integer(b$draws), usable = as.integer(b$usable),
                              not_converged = if ("not_converged" %in% names(b))
                                as.integer(b$not_converged) else NA_integer_,
                              lo = b$lo, hi = b$hi)
    } else if (N_BOOT > 0) {
      bd <- boot_draws(d, N_BOOT, paste(v, ru$rule, sep = ", "), quadratic = ru$quadratic)
      out <- out %>% mutate(draws = bd$draws, usable = bd$usable,
                            not_converged = bd$not_converged, lo = bd$lo, hi = bd$hi)
    }
    out
  })
})
write_csv(T5d, file.path(TAB_DIR, "T5d_span_sensitivity.csv"))
print(as.data.frame(T5d %>% transmute(rule, factor = var, people = n_people,
                                      correlation = round(slope_correlation, 3),
                                      converged, draws, usable,
                                      lo = round(lo, 3), hi = round(hi, 3))),
      row.names = FALSE)

say("Written to %s", TAB_DIR)
})


# #############################################################################
# Section 6.2. The joint growth model on a clinical scale and by trial arm (eTable 19)
# Developed as 14_review_checks.R version 1.1.0, section 3 and its bootstraps
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


L    <- readRDS(file.path(OUT_DIR, "a4_long.rds"))
lab  <- readRDS(file.path(OUT_DIR, "a4_lab.rds"))
read_tab <- function(name) read_csv(file.path(TAB_DIR, paste0(name, ".csv")), show_col_types = FALSE)
T5b  <- read_tab("T5b_joint_growth")
T5c  <- if (file.exists(file.path(TAB_DIR, "T5c_joint_bootstrap.csv"))) read_tab("T5c_joint_bootstrap") else NULL

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
# 3. The joint growth model: clinical scale and trial arm
# =============================================================================
# The model of step 06, section 3, restated as step 12 does, and checked against
# T5b, overall and in each arm, before anything is reported.
#   The regression of the rate of p-tau217 on the rate of eGFR, from the same
#   covariance matrix as the correlation: cov / var of the eGFR rate on the
#   standardised scale, converted back with the standard deviations used to
#   standardise each series. It rescales the correlation between participants'
#   true rates; it is not the effect of one person's decline. eGFR is signed so
#   that a positive rate is faster decline, so the result is the percentage
#   change per year in p-tau217 per 1 mL/min/1.73 m2 per year faster decline.
#   The difference in the correlation between the solanezumab and placebo arms
#   (step 06 fits each arm separately; eTable 19), with an interval from
#   resampling participants within each arm.
#   The model with trial arm on the level and rate of both series is fitted
#   once, for its arm-by-time coefficients: whether solanezumab changed either
#   rate.
say("3. The joint growth model: clinical scale and trial arm")
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
# people restricts the model to a group, as step 06 does for each arm: the
# covariates and both series are then standardised within that group.
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
# The correlation, the regression of the p-tau217 rate on the eGFR rate in
# percent per year per mL/min/1.73 m2 per year, the standard deviation of the
# true eGFR rate in mL/min/1.73 m2 per year, and the regression per standard
# deviation of the eGFR rate.
slope_stats <- function(m, sd_bio, sd_fac) {
  na <- c(r = NA_real_, pct_per_unit = NA_real_, sd_egfr_rate = NA_real_, pct_per_sd = NA_real_)
  if (is.null(m)) return(na)
  vc <- as.matrix(nlme::getVarCov(m)); nm <- colnames(vc)
  i <- grep("ix:t|t:ix", nm)[1]; j <- grep("ip:t|t:ip", nm)[1]
  if (is.na(i) || is.na(j)) return(na)
  beta_raw <- vc[i, j] / vc[i, i] * sd_bio / sd_fac
  sd_rate  <- sqrt(vc[i, i]) * sd_fac
  c(r = vc[i, j] / sqrt(vc[i, i] * vc[j, j]), pct_per_unit = pct(beta_raw),
    sd_egfr_rate = sd_rate, pct_per_sd = pct(beta_raw * sd_rate))
}

d3 <- joint_data("egfr")
sd_bio <- attr(d3, "sd_bio"); sd_fac <- attr(d3, "sd_fac")
arm_tab <- d3 %>% distinct(BID, arm) %>% count(arm)
say("  people in the joint model by arm: %s", paste(arm_tab$arm, arm_tab$n, collapse = ", "))
ARM_TERMS <- c("sola", if (any(arm_tab$arm == "LEARN")) "learn")
m_main <- fit_joint(d3)
s_main <- slope_stats(m_main, sd_bio, sd_fac)
ref_r <- T5b$slope_correlation[T5b$var == "egfr" & T5b$group == "everyone"]
say("  check: eGFR correlation %.4f here, %.4f in T5b_joint_growth.csv", s_main[["r"]], ref_r)
if (length(ref_r) != 1 || is.na(s_main[["r"]]) || abs(s_main[["r"]] - ref_r) > 1e-4)
  stop("The restated joint model does not reproduce step 06. Nothing further is reported.")
d_pl <- joint_data("egfr", L$BID[L$arm == "Placebo"])
d_so <- joint_data("egfr", L$BID[L$arm == "Solanezumab"])
r_pl <- slope_correlation(fit_joint(d_pl))
r_so <- slope_correlation(fit_joint(d_so))
ref_pl <- T5b$slope_correlation[T5b$var == "egfr" & T5b$group == "Placebo"]
ref_so <- T5b$slope_correlation[T5b$var == "egfr" & T5b$group == "Solanezumab"]
say("  check: placebo %.4f here, %.4f in T5b; solanezumab %.4f here, %.4f in T5b",
    r_pl, ref_pl, r_so, ref_so)
if (length(ref_pl) != 1 || length(ref_so) != 1 || is.na(r_pl) || is.na(r_so) ||
    abs(r_pl - ref_pl) > 1e-4 || abs(r_so - ref_so) > 1e-4)
  stop("The restated arm models do not reproduce step 06. Nothing further is reported.")

# The arm-by-time coefficients, standardised units per year, with their model
# intervals: how much faster each series changed in the solanezumab arm.
m_arm <- fit_joint(d3, extra = ARM_TERMS)
arm_fx <- if (is.null(m_arm)) tibble() else {
  tt <- summary(m_arm)$tTable
  rn <- rownames(tt)
  pick <- rn[grepl("sola", rn) & grepl("(^|:)t(:|$)", rn)]
  tibble(term = pick, estimate = tt[pick, "Value"], se = tt[pick, "Std.Error"]) %>%
    mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se,
           series = if_else(grepl("ip", term), "p-tau217", "eGFR, signed so higher is faster decline"),
           converged = converged(m_arm), r_with_arm = slope_correlation(m_arm))
}
if (nrow(arm_fx)) put_rev(arm_fx, "review_3b_arm_by_time")

# The bootstraps, as in step 06: participants resampled with replacement, each
# drawn participant a new identifier, seed 1, and only converged fits used.
resample_people <- function(d, people) {
  ids <- sample(people, length(people), replace = TRUE)
  tibble(BID = ids, .new = seq_along(ids)) %>%
    inner_join(d, by = "BID", relationship = "many-to-many") %>%
    mutate(BID = .new) %>% select(-.new)
}
boot_main <- function(d, B) {
  people <- unique(d$BID)
  set.seed(1)
  t0 <- Sys.time()
  res <- vapply(seq_len(B), function(b) {
    if (b %% 50 == 0) say("  main model: %d of %d, %.0f minutes so far", b, B,
                          as.numeric(difftime(Sys.time(), t0, units = "mins")))
    dd <- resample_people(d, people)
    m <- tryCatch(fit_joint(dd), error = function(e) NULL)
    if (!converged(m)) return(c(r = NA_real_, pct_per_unit = NA_real_,
                                sd_egfr_rate = NA_real_, pct_per_sd = NA_real_))
    slope_stats(m, sd_bio, sd_fac)
  }, numeric(4))
  list(draws = t(res), usable = sum(!is.na(res["r", ])))
}
# Resampling within each arm keeps the arm sizes fixed. A draw is usable only
# if both arm fits converged.
boot_arms <- function(dp, ds, B) {
  pp <- unique(dp$BID); ps <- unique(ds$BID)
  set.seed(1)
  t0 <- Sys.time()
  res <- vapply(seq_len(B), function(b) {
    if (b %% 50 == 0) say("  arm difference: %d of %d, %.0f minutes so far", b, B,
                          as.numeric(difftime(Sys.time(), t0, units = "mins")))
    mp <- tryCatch(fit_joint(resample_people(dp, pp)), error = function(e) NULL)
    ms <- tryCatch(fit_joint(resample_people(ds, ps)), error = function(e) NULL)
    if (!converged(mp) || !converged(ms)) return(c(placebo = NA_real_, solanezumab = NA_real_,
                                                   difference = NA_real_))
    rp <- slope_correlation(mp); rs <- slope_correlation(ms)
    c(placebo = rp, solanezumab = rs, difference = rs - rp)
  }, numeric(3))
  list(draws = t(res), usable = sum(!is.na(res["difference", ])))
}
boot_interval <- function(bd, stat) {
  if (is.null(bd)) return(c(NA_real_, NA_real_))
  unname(stats::quantile(bd$draws[, stat], c(0.025, 0.975), na.rm = TRUE))
}
build_R3 <- function(bd_main, bd_arm) {
  main <- tibble(model = "main",
                 quantity = c("slope correlation",
                              "p-tau217 change, % per year, per 1 mL/min/1.73 m2 per year faster eGFR decline",
                              "standard deviation of the true eGFR rate, mL/min/1.73 m2 per year",
                              "p-tau217 change, % per year, per SD faster eGFR decline"),
                 stat = c("r", "pct_per_unit", "sd_egfr_rate", "pct_per_sd"),
                 estimate = unname(s_main[c("r", "pct_per_unit", "sd_egfr_rate", "pct_per_sd")]),
                 n_people = n_distinct(d3$BID),
                 draws = if (is.null(bd_main)) 0L else as.integer(REVIEW_BOOT),
                 usable = if (is.null(bd_main)) NA_integer_ else bd_main$usable)
  main$lo <- vapply(main$stat, function(s) boot_interval(bd_main, s)[1], numeric(1))
  main$hi <- vapply(main$stat, function(s) boot_interval(bd_main, s)[2], numeric(1))
  arms <- tibble(model = "by arm",
                 quantity = c("slope correlation, placebo", "slope correlation, solanezumab",
                              "difference, solanezumab minus placebo"),
                 stat = c("placebo", "solanezumab", "difference"),
                 estimate = c(r_pl, r_so, r_so - r_pl),
                 n_people = c(n_distinct(d_pl$BID), n_distinct(d_so$BID),
                              n_distinct(d_pl$BID) + n_distinct(d_so$BID)),
                 draws = if (is.null(bd_arm)) 0L else as.integer(REVIEW_ARM_BOOT),
                 usable = if (is.null(bd_arm)) NA_integer_ else bd_arm$usable)
  arms$lo <- vapply(arms$stat, function(s) boot_interval(bd_arm, s)[1], numeric(1))
  arms$hi <- vapply(arms$stat, function(s) boot_interval(bd_arm, s)[2], numeric(1))
  bind_rows(main, arms) %>% mutate(sd_ptau217_log = sd_bio, sd_egfr = sd_fac)
}
R3 <- build_R3(NULL, NULL)
put_rev(R3, "review_3_joint_model")
print(as.data.frame(R3 %>% transmute(model, quantity, estimate = round(estimate, 3))), row.names = FALSE)

# =============================================================================
# 3, continued. The bootstraps of the joint growth models, the slow part
# =============================================================================
# Run last, so that a problem in the other sections shows within minutes.
# review_3_joint_model.csv is rewritten with the intervals.
say("3, continued. Bootstrapping: %d draws for the main model, %d for the arm difference",
    REVIEW_BOOT, REVIEW_ARM_BOOT)
bd_main <- if (REVIEW_BOOT > 0) boot_main(d3, REVIEW_BOOT) else NULL
if (!is.null(bd_main) && REVIEW_BOOT == T5c$draws[T5c$var == "egfr" & T5c$group == "everyone"][1])
  say("  check: correlation interval %.3f to %.3f here, %.3f to %.3f in T5c",
      boot_interval(bd_main, "r")[1], boot_interval(bd_main, "r")[2],
      T5c$lo[T5c$var == "egfr" & T5c$group == "everyone"][1],
      T5c$hi[T5c$var == "egfr" & T5c$group == "everyone"][1])
bd_arm <- if (REVIEW_ARM_BOOT > 0) boot_arms(d_pl, d_so, REVIEW_ARM_BOOT) else NULL
R3 <- build_R3(bd_main, bd_arm)
put_rev(R3, "review_3_joint_model")
print(as.data.frame(R3 %>% transmute(model, quantity, estimate = round(estimate, 3),
                                     lo = round(lo, 3), hi = round(hi, 3), usable)), row.names = FALSE)
})
