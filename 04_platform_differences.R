# =============================================================================
# 04_platform_differences.R
# Stage 4. Research question 2: does a factor change what is in the blood, or
# what a particular assay reads? Differences between platforms on the same aliquot.
#   4.1  every platform pair and factor, the subsamples, the MIMIC model and the
#        C2N detection limit (Table 3, Figure 1, eTables 10 to 12, 14 to 17)
#   4.2  creatine kinase interactions, the C2N differences as a difference of
#        coefficients, and the date of the draw (eTables 36 to 38)
# Writes paper1_out/tables/T2*.csv, T3, T6, T6b, review/review_1*, review_2*.
# =============================================================================

# #############################################################################
# Section 4.1. Differences between platforms
# Developed as 04_rq2_platform.R version 1.4.1, whole
# #############################################################################
local({
SCRIPT_VERSION <- "1.4.1"
cat("04_rq2_platform.R version", SCRIPT_VERSION, "\n")
# =============================================================================
# 04_rq2_platform.R
# Research question 2, the centrepiece: does a factor change what is in the
# blood, or what a particular assay reads?
#
# Every platform measured the same aliquot. A factor that changes the
# circulating concentration moves all platforms of an analyte together. A
# factor that acts on one assay opens a difference between platforms. The
# difference is taken on the log scale, so it is a proportional difference;
# standardised differences are not used because they attenuate the noisier
# assay and can manufacture a difference where there is none.
#
# Three things are reported for every platform pair and every factor.
#   the difference itself, with standard errors clustered by person
#   a verdict with three levels (eMethods 8.6 and 9):
#     "platforms agree to within the margin" when the 90% interval lies inside
#       a margin of 5% per standard deviation;
#     "platforms differ" when the false-discovery-rate adjusted q is below 0.05
#       and the pair is not equivalent;
#     "no difference detected" otherwise. A difference that is merely not
#       significant is not read as agreement
#   a slope-corrected difference, for assays that are not proportional to each
#     other, computed across a range of plausible error-variance ratios rather
#     than assuming one
#
# The latent concentration model is the omnibus version of the same question:
# the factor acts on a latent concentration, and a direct effect of the factor
# on the platforms beyond it is the test of an assay-specific effect. Only the
# omnibus test is reported, for the reason given where the model is fitted.
#
# Reads    paper1_out/adni_visit.rds
# Writes   paper1_out/tables/
#            T2_platform_difference.csv        every pair, every factor
#            T2b_equivalence_summary.csv       verdict counts by analyte and factor
#            T2c_sensitivity.csv               statins, sex, first visit, recency
#            T3_latent_direct_effects.csv      the latent model, omnibus test
#            T6_floor_sensitivity.csv          the C2N detection limit
#            T6b_difference_floor_excluded.csv C2N pairs, floor values set aside
#
# Needs    lavaan for the latent model, survival for the censored regression
# Takes about two minutes.
#
# Changes in 1.2.0
#   verdicts use the adjusted q, not the unadjusted p, and every sensitivity
#     table is adjusted within itself (eMethods 9)
#   Deming slopes and slope-corrected differences for pairs with C2N p-tau217
#     use only samples above the C2N detection limit (eMethods 5)
#   statin subsamples use statin status at each draw when the data carry it;
#     "first visit" is each person's first draw with both platforms of the pair
#   the latent model reports only the omnibus test, which is identified
#   T6 counts the samples the censored regression actually used
#
# Changes in 1.4.1: a stale comment corrected; no change to any result.
#
# Changes in 1.4.0
#   a difference involving C2N p-tau217 is fitted by censored regression, as
#     every other C2N p-tau217 estimate is (eMethods 5): a C2N value below the
#     limit makes the difference left-censored at log(1.30) minus the other
#     platform's log value. Before 1.4.0 the difference used the reported 0.65
#     as if it were a measurement. T6b, the floor-excluded check, stays least
#     squares on the samples above the limit.
#   the laboratory-recency subsamples measure the gap from each draw to the
#     laboratory date, as eMethods 11 states, not from the first draw
#   the latent model for p-tau217 no longer reads the C2N value 0.65 as a
#     measurement. lavaan cannot fit a censored indicator with this estimator,
#     so the main model leaves C2N out (Fujirebio, ALZpath, Janssen) and a
#     sensitivity model keeps C2N with its values below the limit set to
#     missing. Other analytes are unchanged, as none has a detection floor.
# =============================================================================

if (!exists("PAPER1_VERSION")) source("00_setup.R")
HAVE_LAVAAN <- requireNamespace("lavaan", quietly = TRUE)
if (!HAVE_LAVAAN) say("lavaan is not installed, so the latent model is skipped.")

d <- readRDS(file.path(OUT_DIR, "adni_visit.rds"))
short <- function(x) gsub("_plasma", "", x)

ANALYTES <- list(
  pTau217 = c("C2N_plasma_ptau217", "Fuji_plasma_ptau217", "AlzPath_plasma_ptau217",
              "Janssen_plasma_ptau217"),
  pTau181 = c("Roche_plasma_ptau181", "QX_plasma_ptau181"),
  NfL     = c("Roche_plasma_NfL", "QX_plasma_NfL"),
  GFAP    = c("Roche_plasma_GFAP", "QX_plasma_GFAP"),
  Abeta40 = c("C2N_plasma_Abeta40", "Fuji_plasma_Ab40", "Roche_plasma_Ab40", "QX_plasma_Ab40"),
  Abeta42 = c("C2N_plasma_Abeta42", "Fuji_plasma_Ab42", "Roche_plasma_Ab42", "QX_plasma_Ab42"),
  Ab42_40 = c("C2N_plasma_Abeta42_Abeta40", "Fuji_plasma_Ab42_Ab40",
              "Roche_plasma_Ab42_Ab40", "QX_plasma_Ab42_Ab40"))
ANALYTES <- lapply(ANALYTES, function(m) m[m %in% names(d)])
FAC <- FACTORS %>% filter(paste0("z_", var) %in% names(d))

# ---- The C2N detection limit, for the slope ----------------------------------
# C2N reports p-tau217 below its detection limit as one constant (0.65 pg/mL).
# eMethods 5 says such a value is never treated as a measurement. A Deming slope
# fitted through several hundred copies of one constant is pulled towards
# whatever line passes through that cluster, and in these data it moved every
# C2N slope by 0.2 to 0.36. So for any pair that includes C2N p-tau217 the slope,
# and the slope-corrected difference that uses it, are estimated only on the
# samples above the limit. Pairs without C2N p-tau217 use every sample with both
# values, as before. The uncorrected difference is fitted by censored
# regression (fit_diff_censored below), and T6b repeats it with the floor
# values set aside.
C2N_P217 <- "C2N_plasma_ptau217"
c2n_below_limit <- function(data) {
  if ("c2n_at_floor" %in% names(data)) coalesce(as.logical(data$c2n_at_floor), FALSE)
  else !is.na(data[[C2N_P217]]) & data[[C2N_P217]] <= log(C2N_LOD) + 1e-8
}

# ---- The slope correction ----------------------------------------------------
# Two assays that are not proportional to each other can show a difference that
# is an artefact of scale rather than of the factor. Correcting for that needs a
# Deming slope, which depends on the ratio of the two platforms' measurement
# error variances. That ratio is not known here, so rather than assume a value
# the difference is corrected across a range of plausible ratios and reported
# for each. A ratio of 1.00 is the orthogonal regression, which assumes the two
# assays are equally noisy.
LAMBDAS <- c(0.50, 1.00, 2.00)
lam_tag <- function(l) sprintf("%03d", round(l * 100))

# ---- The verdict and the false discovery rate --------------------------------
# eMethods 9: Benjamini-Hochberg within each family, and a pair is said to
# differ only when its adjusted q is below 0.05. The unadjusted p stays in the
# table beside q, because the adjustment depends on where the family boundary is
# drawn. A pair can be both equivalent and, with a very narrow interval,
# significantly different from zero; equivalence takes precedence, because a
# difference shown to lie inside the margin is by definition too small to matter.
VERDICT_AGREE  <- "platforms agree to within the margin"
VERDICT_DIFFER <- "platforms differ"
VERDICT_NONE   <- "no difference detected"
add_verdict <- function(tab, ...) {
  tab %>%
    group_by(...) %>%
    mutate(p_fdr = p.adjust(p, "BH")) %>%
    ungroup() %>%
    mutate(verdict = case_when(equivalent ~ VERDICT_AGREE,
                               !is.na(p_fdr) & p_fdr < 0.05 ~ VERDICT_DIFFER,
                               TRUE ~ VERDICT_NONE))
}

# ---- The between-platform difference -----------------------------------------
# pair_rows, when given, is a function(data, pair) that picks the rows for one
# pair. It is needed where the subsample depends on the pair, as the first draw
# with both platforms does.
# extra_covs adds covariates to every model of the table, for the statin
# adjustment; the Deming slope is unaffected, as it does not use covariates.
# A difference involving C2N p-tau217 is left-censored wherever the C2N value is
# below the limit. C2N is the first platform of each of its pairs (ANALYTES lists
# it first and combn keeps the order), so the difference is C2N minus the other
# platform and its censoring point is log(1.30) minus the other platform's log
# value, which fit_tobit takes as 1.30 / exp(other) on the original scale. The
# fit is given only the complete cases, so that its counts are those of the
# fitted sample. censor_c2n = FALSE is for data with the floor values removed.
fit_diff_censored <- function(dp, pr, zv, covs) {
  stopifnot(pr[1] == C2N_P217)
  f <- as.formula(sprintf("diff ~ %s + %s", zv, covs))
  vars <- unique(c(all.vars(f), "RID", pr[2]))
  dd <- as.data.frame(dp)
  dd <- dd[stats::complete.cases(dd[, vars, drop = FALSE]), , drop = FALSE]
  rownames(dd) <- NULL
  cens <- c2n_below_limit(dd)
  empty <- tibble(n_obs = nrow(dd), n_people = dplyr::n_distinct(dd$RID),
                  b = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_,
                  p = NA_real_, level = 0.95, n_censored = sum(cens),
                  method = "censored regression")
  if (nrow(dd) < 10 || sum(!cens) < 10) return(empty)
  r <- tryCatch(fit_tobit(f, dd, zv, censor_at = C2N_LOD / exp(dd[[pr[2]]]),
                          censored = cens),
                error = function(e) NULL)
  if (is.null(r) || is.na(r$b)) return(empty)
  empty %>% mutate(b = r$b, se = r$se, lo = r$lo, hi = r$hi, p = r$p)
}

difference_table <- function(data, analytes, factors, label = "main", pair_rows = NULL,
                             extra_covs = NULL, censor_c2n = TRUE) {
  covs <- paste(c(COVS_ADNI, extra_covs), collapse = " + ")
  purrr::map_dfr(names(analytes), function(an) {
    ms <- analytes[[an]]
    if (length(ms) < 2) return(NULL)
    purrr::map_dfr(combn(ms, 2, simplify = FALSE), function(pr) {
      dp <- if (is.null(pair_rows)) data else pair_rows(data, pr)
      # The data for the slope and the corrected difference: the floor values
      # of C2N p-tau217 are set aside when that platform is in the pair.
      dk <- dp
      uses_c2n <- C2N_P217 %in% pr
      if (uses_c2n) dk[[C2N_P217]][c2n_below_limit(dk)] <- NA_real_
      slope_sample <- if (uses_c2n) "C2N p-tau217 above the detection limit"
                      else "all samples with both values"
      slope_n <- sum(stats::complete.cases(dk[[pr[1]]], dk[[pr[2]]]))
      ks <- vapply(LAMBDAS, function(l) deming_slope(dk[[pr[2]]], dk[[pr[1]]], lambda = l),
                   numeric(1))
      dp$diff <- dp[[pr[1]]] - dp[[pr[2]]]
      for (i in seq_along(LAMBDAS))
        dp[[paste0("diff_", lam_tag(LAMBDAS[i]))]] <- dk[[pr[1]]] - ks[i] * dk[[pr[2]]]
      purrr::map_dfr(seq_len(nrow(factors)), function(j) {
        zv <- paste0("z_", factors$var[j])
        g1 <- if (uses_c2n && censor_c2n) fit_diff_censored(dp, pr, zv, covs) else
          fit_cl(as.formula(sprintf("diff ~ %s + %s", zv, covs)), dp, zv) %>%
            mutate(n_censored = NA_integer_, method = "least squares")
        sc <- purrr::map(seq_along(LAMBDAS), function(i) {
          tg <- lam_tag(LAMBDAS[i])
          g2 <- fit_cl(as.formula(sprintf("diff_%s ~ %s + %s", tg, zv, covs)), dp, zv)
          out <- tibble(ks[i], pct(g2$b), pct(g2$lo), pct(g2$hi))
          names(out) <- paste0(c("slope_", "corrected_", "corrected_lo_", "corrected_hi_"), tg)
          out
        }) %>% bind_cols()
        bind_cols(
          tibble(sample = label, analyte = an,
                 pair = paste(short(pr[1]), "minus", short(pr[2])),
                 factor = factors$label[j], var = factors$var[j], role = factors$role[j]),
          g1,
          equivalence(g1$b, g1$se),
          sc,
          tibble(slope_sample = slope_sample, slope_n = slope_n))
      })
    })
  }) %>%
    mutate(pct_diff = pct(b), pct_lo = pct(lo), pct_hi = pct(hi))
}

# Main analysis. The family is all platform pairs by all factors of one role.
T2 <- difference_table(d, ANALYTES, FAC) %>% add_verdict(role)
write_csv(T2, file.path(TAB_DIR, "T2_platform_difference.csv"))

# The summary counts all three verdicts, so that the pairs with no difference
# detected are shown rather than left to be inferred. range_of_90pct_limits is
# the lowest lower and the highest upper 90% limit across the pairs; it is the
# outer envelope of the intervals, not an estimated difference.
T2b <- T2 %>%
  filter(role == "primary") %>%
  group_by(analyte, factor = sub(" \\(.*", "", factor)) %>%
  summarise(pairs = n(),
            equivalent = sum(verdict == VERDICT_AGREE),
            differ = sum(verdict == VERDICT_DIFFER),
            no_difference_detected = sum(verdict == VERDICT_NONE),
            range_of_90pct_limits = sprintf("%.2f to %.2f", min(lo90), max(hi90)),
            .groups = "drop") %>%
  arrange(factor, analyte)
write_csv(T2b, file.path(TAB_DIR, "T2b_equivalence_summary.csv"))
say("Verdicts by analyte, primary factors:")
print(as.data.frame(T2b), row.names = FALSE)

say("p-tau217, primary factors, every pair:")
print(as.data.frame(T2 %>% filter(analyte == "pTau217", role == "primary") %>%
  transmute(pair, factor = sub(" \\(.*", "", factor),
            estimate = sprintf("%.2f (%.2f to %.2f)", pct_diff, pct_lo, pct_hi),
            p = signif(p, 2), q = signif(p_fdr, 2), verdict)),
  row.names = FALSE)

# ---- Is the difference robust ------------------------------------------------
# Each of these repeats the creatine kinase result in a different sample. The
# question is whether the difference survives.
PAIRS_OF_INTEREST <- list(pTau181 = ANALYTES$pTau181,
                          pTau217 = ANALYTES$pTau217)

# Statin use can change between draws, so the subsample is defined at each draw.
# A person-level indicator fixed at the first draw would count a later draw taken
# after a statin was started as a non-user's draw. The per-draw indicator comes
# from step 01; with older data that lack it, the first-draw indicator is used
# and the run says so.
STATIN_VAR <- if ("statin_at_draw" %in% names(d)) "statin_at_draw" else "statin_by_first_plasma"
if (STATIN_VAR != "statin_at_draw" && STATIN_VAR %in% names(d)) {
  msg <- paste("adni_visit.rds has no statin_at_draw, so the statin subsamples fall back",
               "to statin status at the first draw (statin_by_first_plasma). Rerun step 01.")
  say("WARNING: %s", msg)
  warning(msg, call. = FALSE)
}

# The first visit is each person's first draw that has both platforms of the
# pair. visit_no == 1 would drop anyone whose first draw was on other platforms
# only, and would give a different sample for no reason connected to the pair.
first_draw_with_pair <- function(data, pr) {
  data %>%
    filter(!is.na(.data[[pr[1]]]), !is.na(.data[[pr[2]]])) %>%
    group_by(RID) %>%
    filter(visit_date == min(visit_date)) %>%
    slice(1) %>%
    ungroup()
}

sens <- list()
if (STATIN_VAR %in% names(d)) {
  sens[["statin non-users"]] <- list(data = d %>% filter(.data[[STATIN_VAR]] == 0),
                                     rule = sprintf("%s == 0", STATIN_VAR))
  sens[["statin users"]]     <- list(data = d %>% filter(.data[[STATIN_VAR]] == 1),
                                     rule = sprintf("%s == 1", STATIN_VAR))
}
# Statin use at each draw as a covariate, and draws on a kidney-acting drug set
# aside: the two planned checks that the creatine kinase differences are not
# produced by statin myopathy or by drugs acting on the kidney.
if (STATIN_VAR %in% names(d))
  sens[["statin adjusted"]] <- list(data = d %>% filter(!is.na(.data[[STATIN_VAR]])),
                                    extra_covs = STATIN_VAR,
                                    rule = sprintf("%s added as a covariate", STATIN_VAR))
if ("renal_med_at_draw" %in% names(d)) {
  sens[["no kidney-acting medication"]] <- list(data = d %>% filter(renal_med_at_draw == 0),
                                                rule = "renal_med_at_draw == 0")
  sens[["kidney-acting medication"]]    <- list(data = d %>% filter(renal_med_at_draw == 1),
                                                rule = "renal_med_at_draw == 1")
} else say("  renal_med_at_draw is not in the data; rerun step 01 for the kidney-medication check.")
sens[["men"]]         <- list(data = d %>% filter(female == 0), rule = "female == 0")
sens[["women"]]       <- list(data = d %>% filter(female == 1), rule = "female == 1")
sens[["first visit"]] <- list(data = d, pair_rows = first_draw_with_pair,
                              rule = "first draw with both platforms of the pair")
# Recency is the gap between each draw and the laboratory date, not the gap
# from the first draw: every draw enters the model, and a person whose
# screening laboratory was a month before the first draw is years from a later
# one (the same definition as T1f in step 03).
d$lab_gap_draw_years <- abs(as.numeric(d$lab_date - d$visit_date)) / YEAR_DAYS
sens[["laboratory within 1 year"]]  <- list(data = d %>% filter(!is.na(lab_gap_draw_years),
                                                                lab_gap_draw_years <= 1),
                                            rule = "laboratory within 1 year of the draw")
sens[["laboratory within 3 years"]] <- list(data = d %>% filter(!is.na(lab_gap_draw_years),
                                                                lab_gap_draw_years <= 3),
                                            rule = "laboratory within 3 years of the draw")

# The "all" rows restate the main analysis, so they carry its q and verdict
# from T2. Each subsample is its own sensitivity analysis and is adjusted within
# itself, by role (eMethods 9), so that no subsample can change another's q or
# the main analysis's.
T2c_all <- T2 %>%
  filter(analyte %in% names(PAIRS_OF_INTEREST), var == "l_ck") %>%
  mutate(sample = "all", subsample_rule = "all draws",
         fdr_family = "main analysis (T2), all pairs by all factors of one role")
T2c_sub <- purrr::map_dfr(names(sens), function(nm)
  difference_table(sens[[nm]]$data, PAIRS_OF_INTEREST,
                   FAC %>% filter(var == "l_ck"), label = nm,
                   pair_rows = sens[[nm]]$pair_rows,
                   extra_covs = sens[[nm]]$extra_covs) %>%
    mutate(subsample_rule = sens[[nm]]$rule) %>%
    add_verdict(role) %>%
    mutate(fdr_family = "within this subsample, pairs by factors of one role"))
T2c <- bind_rows(T2c_all, T2c_sub)
write_csv(T2c, file.path(TAB_DIR, "T2c_sensitivity.csv"))
say("Creatine kinase and p-tau181, in every subsample:")
print(as.data.frame(T2c %>% filter(analyte == "pTau181") %>%
  transmute(sample, n_people, estimate = sprintf("%.2f (%.2f to %.2f)", pct_diff, pct_lo, pct_hi),
            q = signif(p_fdr, 2), verdict)),
  row.names = FALSE)

# ---- The latent concentration model ------------------------------------------
# The factor acts on a latent concentration that all platforms measure, and an
# assay-specific effect is a direct path from the factor to one platform.
#
# Only the omnibus test is reported. With a direct path allowed on every
# platform but one, the effect on the latent concentration (the "common"
# effect) is identified only by assuming that the remaining platform has no
# assay-specific effect, and its size, even its sign, changes with the platform
# chosen: for creatine kinase and p-tau217 it ran from -0.6% to 12.8% on the C2N
# scale. The model with no direct paths gives a single common effect, but that
# is the model the omnibus test rejects, and a single direct path fitted in its
# own model carries a common effect of its own. Neither can be read as the
# split of an association into a blood part and an assay part. The omnibus
# test, whether there is any assay-specific effect at all, is the same whichever
# platform is the reference, so it is the identified quantity and it is what is
# reported, with the fit statistics it is built from.
#
# The omnibus test compares the model with direct paths on every platform but the
# first against the model with none, by the scaled (Satorra-Bentler 2001)
# difference test. q is Benjamini-Hochberg across all analytes by all factors of
# one role, the same family as the pairwise differences.
mimic_omnibus <- function(data, ms, zv) {
  if (!HAVE_LAVAAN || length(ms) < 3) return(NULL)
  # Every covariate is required here, so that the diagnosis dummies built below
  # line up row for row with the data they are bound to.
  dd <- data %>%
    select(RID, all_of(ms), all_of(zv), age_at_visit, female, e4, centiloid, dx, plasma_time) %>%
    filter(!is.na(.data[[zv]]), !is.na(centiloid), !is.na(age_at_visit), !is.na(e4),
           !is.na(female), !is.na(dx), !is.na(plasma_time))
  if (nrow(dd) < 50) return(NULL)
  ind <- paste0("y", seq_along(ms))
  names(dd)[match(ms, names(dd))] <- ind
  # The indicators stay on the log scale, with the first as the marker (loading
  # fixed at 1), as in eMethods 8.8.
  dx_d <- model.matrix(~ dx, dd)[, -1, drop = FALSE]
  colnames(dx_d) <- paste0("dx", seq_len(ncol(dx_d)))
  dd <- bind_cols(dd, as.data.frame(dx_d))
  dd$x <- dd[[zv]]
  covs <- paste(c("age_at_visit", "female", "e4", "centiloid", colnames(dx_d), "plasma_time"),
                collapse = " + ")
  base <- sprintf("eta =~ %s\n eta ~ x + %s", paste(ind, collapse = " + "), covs)
  # Warnings are recorded rather than hidden, so that a fit that did not
  # converge, or has a negative variance, is visible in the table.
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
         chisq_no_direct = fm(f0$fit, "chisq.scaled"), df_no_direct = fm(f0$fit, "df"),
         chisq_all_direct = fm(f_all$fit, "chisq.scaled"), df_all_direct = fm(f_all$fit, "df"),
         omnibus_chisq = if (is.null(lrt)) NA_real_ else lrt$`Chisq diff`[2],
         omnibus_df = if (is.null(lrt)) NA_real_ else lrt$`Df diff`[2],
         omnibus_p = if (is.null(lrt)) NA_real_ else lrt$`Pr(>Chisq)`[2],
         fit_warnings = paste(unique(c(f0$warnings, f_all$warnings)), collapse = " | "))
}
run_mimic <- function(data, analytes, factors, specification)
  purrr::map_dfr(names(analytes), function(an)
    purrr::map_dfr(seq_len(nrow(factors)), function(j) {
      r <- mimic_omnibus(data, analytes[[an]], paste0("z_", factors$var[j]))
      if (is.null(r)) return(NULL)
      bind_cols(tibble(specification = specification, analyte = an,
                       platforms = paste(short(analytes[[an]]), collapse = ", "),
                       factor = factors$label[j], var = factors$var[j],
                       role = factors$role[j]), r)
    }))

# C2N p-tau217 below its limit is reported as 0.65, which is not a measurement
# (eMethods 5), and lavaan cannot fit a censored indicator with this estimator.
# Main: p-tau217 without C2N, on the three platforms with no floor. Sensitivity:
# C2N kept, with its values below the limit set to missing, so that full
# information maximum likelihood uses its quantified values; that assumes the
# missingness is unrelated to the value, which it is not, and is why this is
# the sensitivity and not the main model. Each is its own false discovery rate
# family by role.
MIMIC_MAIN <- ANALYTES
MIMIC_MAIN$pTau217 <- setdiff(ANALYTES$pTau217, C2N_P217)
d_mimic_sens <- d
d_mimic_sens[[C2N_P217]][c2n_below_limit(d_mimic_sens)] <- NA_real_

if (HAVE_LAVAAN) {
  say("Fitting the latent concentration model")
  T3 <- bind_rows(
    run_mimic(d, MIMIC_MAIN, FAC, "main") %>%
      group_by(role) %>% mutate(omnibus_q = p.adjust(omnibus_p, "BH")) %>% ungroup(),
    run_mimic(d_mimic_sens, ANALYTES["pTau217"], FAC,
              "C2N included, values below the limit set to missing") %>%
      group_by(role) %>% mutate(omnibus_q = p.adjust(omnibus_p, "BH")) %>% ungroup()) %>%
    relocate(omnibus_q, .after = omnibus_p)
  write_csv(T3, file.path(TAB_DIR, "T3_latent_direct_effects.csv"))
  if (any(nzchar(T3$fit_warnings))) say("Some latent model fits gave warnings; see fit_warnings in T3.")
  say("Omnibus test for an assay-specific effect, primary factors (p, q):")
  print(as.data.frame(T3 %>% filter(role == "primary") %>%
    transmute(specification, analyte, factor = sub(" \\(.*", "", factor),
              chisq = round(omnibus_chisq, 2), df = omnibus_df,
              p = signif(omnibus_p, 2), q = signif(omnibus_q, 2))), row.names = FALSE)
}

# ---- The C2N detection limit -------------------------------------------------
# C2N p-tau217 is reported as half the limit of detection when it falls below
# it. Excluding those values selects on the outcome, so the censored regression
# is the main analysis and the exclusion is the sensitivity analysis.
d_excl <- d
d_excl$C2N_plasma_ptau217[c2n_below_limit(d_excl)] <- NA_real_

# The counts are taken from the rows the censored regression can use, those
# with the outcome, the factor and every covariate. Counting before the model
# drops incomplete rows would give the same 1,179 samples for every factor.
tobit_rows_used <- function(formula, data) {
  vars <- unique(c(all.vars(formula), "RID"))
  stats::complete.cases(as.data.frame(data)[vars])
}
T6 <- purrr::map_dfr(seq_len(nrow(FAC)), function(j) {
  zv <- paste0("z_", FAC$var[j])
  fml <- as.formula(sprintf("C2N_plasma_ptau217 ~ %s + %s", zv, COVS_ADNI))
  as_is <- fit_cl(fml, d, zv)
  excl  <- fit_cl(fml, d_excl, zv)
  cens  <- fit_tobit(fml, d, zv, censor_at = C2N_LOD)
  used  <- tobit_rows_used(fml, d)
  tibble(factor = FAC$label[j], var = FAC$var[j], role = FAC$role[j],
         n_obs = sum(used),
         n_censored = sum(used & d$C2N_plasma_ptau217 <= log(C2N_LOD) + 1e-8),
         pct_as_delivered = pct(as_is$b),
         pct_floor_excluded = pct(excl$b),
         pct_censored = pct(cens$b), lo_censored = pct(cens$lo), hi_censored = pct(cens$hi),
         p_censored = cens$p)
}) %>%
  # A sensitivity analysis, adjusted within itself by role.
  group_by(role) %>% mutate(p_fdr_censored = p.adjust(p_censored, "BH")) %>% ungroup()
write_csv(T6, file.path(TAB_DIR, "T6_floor_sensitivity.csv"))
say("C2N p-tau217, the three ways of handling the detection limit:")
print(as.data.frame(T6 %>% filter(role == "primary") %>%
  transmute(factor = sub(" \\(.*", "", factor), n_obs, n_censored,
            as_delivered = round(pct_as_delivered, 2),
            floor_excluded = round(pct_floor_excluded, 2),
            censored = sprintf("%.2f (%.2f to %.2f)", pct_censored, lo_censored, hi_censored))),
  row.names = FALSE)

# The floor also affects the pair differences involving C2N, so those are
# repeated with the floor values set aside. This is a sensitivity analysis, so
# it is adjusted within itself: the three C2N pairs by the factors of one role.
T6b <- difference_table(d_excl, ANALYTES["pTau217"], FAC, label = "C2N floor excluded",
                        censor_c2n = FALSE) %>%
  filter(grepl("C2N", pair)) %>%
  add_verdict(role)
write_csv(T6b, file.path(TAB_DIR, "T6b_difference_floor_excluded.csv"))

say("Written to %s", TAB_DIR)
say("T2 %d rows, T2c %d, T6 %d, T6b %d%s", nrow(T2), nrow(T2c), nrow(T6), nrow(T6b),
    if (HAVE_LAVAAN) sprintf(", T3 %d", nrow(T3)) else "")
})


# #############################################################################
# Section 4.2. Creatine kinase interactions, difference of coefficients, date of the draw (eTables 36 to 38)
# Developed as 14_review_checks.R version 1.1.0, sections 1, 1b and 2
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
T2   <- read_tab("T2_platform_difference")

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
# The between-platform difference, as in step 04
# =============================================================================
# A difference involving C2N p-tau217 is left-censored where the C2N value is
# below the limit, at log(1.30) minus the other platform's log value; any other
# difference is least squares. Standard errors are clustered by person in both.
pair_rows <- function(data, pr, rhs) {
  dd <- as.data.frame(data)
  # The same complete cases as step 04: the terms, the cluster and the two
  # platforms. The C2N flag is not required, and a missing flag counts as not
  # censored, as in step 04.
  vars <- unique(c(all.vars(stats::as.formula(paste("~", rhs))), "RID", pr))
  dd <- dd[stats::complete.cases(dd[, vars, drop = FALSE]), , drop = FALSE]
  rownames(dd) <- NULL
  dd
}
pair_fit <- function(data, pr, rhs, term) {
  dd <- pair_rows(data, pr, rhs)
  dd$diff <- dd[[pr[1]]] - dd[[pr[2]]]
  f <- stats::as.formula(paste("diff ~", rhs))
  if (pr[1] == C2N_P217) {
    cens <- dplyr::coalesce(as.logical(dd$c2n_at_floor), FALSE)
    r <- tryCatch(fit_tobit(f, dd, term, censor_at = C2N_LOD / exp(dd[[pr[2]]]), censored = cens),
                  error = function(e) tibble(b = NA_real_, se = NA_real_, lo = NA_real_,
                                             hi = NA_real_, p = NA_real_))
    method <- "censored regression"
  } else {
    r <- fit_cl(f, dd, term)
    method <- "least squares"
  }
  tibble(n_obs = nrow(dd), n_people = dplyr::n_distinct(dd$RID), b = r$b, se = r$se,
         lo = r$lo, hi = r$hi, p = r$p, method = method)
}
pair_name <- function(pr) paste(short(pr[1]), "minus", short(pr[2]))
check_pair <- function(pr, v, data = adni) {
  mine <- pair_fit(data, pr, sprintf("z_%s + %s", v, COVS_ADNI), paste0("z_", v))
  ref  <- T2 %>% filter(sample == "main", pair == pair_name(pr), var == v)
  if (nrow(ref) != 1 || is.na(mine$b) || abs(mine$b - ref$b) > 1e-6)
    stop(sprintf("Step 14 does not reproduce T2 for %s and %s: %.6f here, %s in T2",
                 pair_name(pr), v, mine$b, paste(signif(ref$b, 6), collapse = ", ")))
  invisible(mine)
}

# The C2N difference a second way (section 1b). The censored difference above
# takes its censoring point, log(1.30) minus the other platform's value, as
# fixed, although that value carries the other platform's measurement error.
# Here the C2N level is fitted by censored regression with the fixed limit of
# 1.30 pg/mL, the other platform's level by least squares, on the same rows and
# with the same terms, and the difference of the two coefficients is the
# estimate. Its interval comes from resampling participants.
c2n_level_coef <- function(dd, rhs, term) {
  dd$.obs <- !dplyr::coalesce(as.logical(dd$c2n_at_floor), FALSE)
  dd$.y   <- ifelse(dd$.obs, dd[[C2N_P217]], log(C2N_LOD))
  f <- stats::as.formula(sprintf("survival::Surv(.y, .obs, type = 'left') ~ %s", rhs))
  fit <- tryCatch(survival::survreg(f, data = dd, dist = "gaussian"), error = function(e) NULL)
  if (is.null(fit) || !term %in% names(stats::coef(fit))) return(NA_real_)
  unname(stats::coef(fit)[term])
}
ls_coef <- function(dd, y, rhs, term) {
  cf <- stats::coef(stats::lm(stats::as.formula(sprintf("`%s` ~ %s", y, rhs)), data = dd))
  if (!term %in% names(cf)) NA_real_ else unname(cf[term])
}
coef_difference <- function(data, pr, rhs, term) {
  dd <- pair_rows(data, pr, rhs)
  fun <- function(x) c(difference = c2n_level_coef(x, rhs, term) - ls_coef(x, pr[2], rhs, term))
  boot_pct(dd, fun) %>%
    transmute(n_obs = nrow(dd), n_people = n_distinct(dd$RID), pct = pct(est),
              lo = pct(lo), hi = pct(hi), usable)
}

# =============================================================================
# 1. Creatine kinase differences by sex, statin use and kidney-acting medication
# =============================================================================
# Each interaction model lets every covariate differ by the modifier, so the
# effect in each group is that of a model fitted in that group alone (exactly
# for least squares; for the censored model the residual spread is shared).
# The effect in each group comes from the same model with the modifier coded
# each way round. The interaction is reported as a percentage: the ratio of
# the per-SD fold-difference in the group coded 1 to that in the group coded
# 0, minus one, times 100. The modifier is dropped from the covariates where it
# is one of them (sex).
# Creatine kinase was measured once, at screening, so the main medication
# modifier is the status at each participant's first draw, the draw nearest
# the laboratory panel; status at each draw, as in eTable 12, is the
# alternative. The six main interactions form one false discovery rate family;
# the alternatives and the two sensitivity models are reported without
# adjustment.
say("1. Creatine kinase differences: interaction tests")
CK_PAIRS <- list(c(C2N_P217, "AlzPath_plasma_ptau217"), P181)
for (pr in CK_PAIRS) check_pair(pr, "l_ck")
say("  check: both creatine kinase differences reproduce T2")
adni <- adni %>% arrange(RID, visit_no) %>% group_by(RID) %>%
  mutate(statin_first = statin_at_draw[1], renal_first = renal_med_at_draw[1]) %>%
  ungroup()
MODIFIERS <- tibble::tribble(
  ~mod,                ~label,                                        ~level0,        ~level1,             ~set,
  "female",            "sex",                                         "men",          "women",             "main",
  "statin_first",      "statin use at the first draw",                "no statin",    "statin",            "main",
  "renal_first",       "kidney-acting medication at the first draw",  "no such drug", "kidney-acting drug", "main",
  "statin_at_draw",    "statin use at each draw",                     "no statin",    "statin",            "alternative",
  "renal_med_at_draw", "kidney-acting medication at each draw",       "no such drug", "kidney-acting drug", "alternative")
interaction_rows <- function(pr, k, extra = "") {
  m <- MODIFIERS$mod[k]
  d <- adni %>% filter(!is.na(.data[[m]]))
  d$.m1 <- as.integer(d[[m]] == 1)
  d$.m0 <- 1L - d$.m1
  covs <- covs_without(m)
  rhs1 <- sprintf("(z_l_ck + %s) * .m1%s", covs, extra)
  rhs0 <- sprintf("(z_l_ck + %s) * .m0%s", covs, extra)
  int <- pair_fit(d, pr, rhs1, "z_l_ck:.m1")
  g0  <- pair_fit(d, pr, rhs1, "z_l_ck")
  g1  <- pair_fit(d, pr, rhs0, "z_l_ck")
  tibble(pair = pair_name(pr), modifier = MODIFIERS$label[k], set = MODIFIERS$set[k],
         method = int$method, n_obs = int$n_obs, n_people = int$n_people,
         n_people_group1 = n_distinct(pair_rows(d, pr, rhs1)$RID[pair_rows(d, pr, rhs1)$.m1 == 1]),
         group0 = MODIFIERS$level0[k],
         group0_pct = pct(g0$b), group0_lo = pct(g0$lo), group0_hi = pct(g0$hi),
         group1 = MODIFIERS$level1[k],
         group1_pct = pct(g1$b), group1_lo = pct(g1$lo), group1_hi = pct(g1$hi),
         interaction_pct = pct(int$b), interaction_lo = pct(int$lo), interaction_hi = pct(int$hi),
         p = int$p, se_log = int$se)
}
R1 <- purrr::map_dfr(CK_PAIRS, function(pr) purrr::map_dfr(seq_len(nrow(MODIFIERS)),
                                                            function(k) interaction_rows(pr, k)))
# Sensitivity: a squared creatine kinase term, so that a curved creatine
# kinase effect cannot appear as a sex interaction (men have higher values).
R1_sq <- purrr::map_dfr(CK_PAIRS, function(pr)
  interaction_rows(pr, which(MODIFIERS$mod == "female"), extra = " + I(z_l_ck^2)")) %>%
  mutate(set = "sensitivity: squared creatine kinase term")
# Sensitivity: the three main modifiers in one model, each interacting with
# creatine kinase, on the participants with all three.
R1_mut <- purrr::map_dfr(CK_PAIRS, function(pr) {
  d <- adni %>% filter(!is.na(statin_first), !is.na(renal_first))
  rhs <- sprintf("z_l_ck * (female + statin_first + renal_first) + %s", covs_without("female"))
  purrr::map_dfr(c("female", "statin_first", "renal_first"), function(m) {
    r <- pair_fit(d, pr, rhs, paste0("z_l_ck:", m))
    tibble(pair = pair_name(pr), modifier = MODIFIERS$label[MODIFIERS$mod == m],
           set = "sensitivity: three modifiers in one model", method = r$method,
           n_obs = r$n_obs, n_people = r$n_people,
           interaction_pct = pct(r$b), interaction_lo = pct(r$lo), interaction_hi = pct(r$hi),
           p = r$p, se_log = r$se)
  })
})
R1 <- bind_rows(R1, R1_sq, R1_mut) %>%
  group_by(set) %>%
  mutate(q = if (first(set) == "main") p.adjust(p, "BH") else NA_real_) %>%
  ungroup() %>%
  # The smallest interaction detectable with 80% power at a two-sided 5% level,
  # from the standard error of the interaction on the log scale.
  mutate(detectable_80pct_power = pct(2.80 * se_log),
         fdr_family = if_else(set == "main", "the six main interactions", "not adjusted"))
put_rev(R1, "review_1_ck_interactions")
print(as.data.frame(R1 %>% transmute(pair, modifier, set,
  group0 = ifelse(is.na(group0), "", paste(group0, fmt(group0_pct, group0_lo, group0_hi))),
  group1 = ifelse(is.na(group1), "", paste(group1, fmt(group1_pct, group1_lo, group1_hi))),
  interaction = fmt(interaction_pct, interaction_lo, interaction_hi),
  p = signif(p, 2), q = signif(q, 2), detectable = round(detectable_80pct_power, 2))), row.names = FALSE)

# =============================================================================
# 1b. The C2N differences as a difference of level coefficients
# =============================================================================
say("1b. The C2N differences as a difference of level coefficients (participant bootstrap)")
C2N_PAIRS <- combn(P217, 2, simplify = FALSE)
C2N_PAIRS <- C2N_PAIRS[vapply(C2N_PAIRS, function(pr) pr[1] == C2N_P217, logical(1))]
PRIMARY <- FACTORS %>% filter(role == "primary", paste0("z_", var) %in% names(adni))
R1b_main <- purrr::map_dfr(C2N_PAIRS, function(pr) purrr::map_dfr(PRIMARY$var, function(v) {
  zv <- paste0("z_", v)
  ref <- T2 %>% filter(sample == "main", pair == pair_name(pr), var == v)
  coef_difference(adni, pr, sprintf("%s + %s", zv, COVS_ADNI), zv) %>%
    mutate(analysis = "main difference", pair = pair_name(pr), var = v,
           censored_difference_pct = ref$pct_diff, censored_lo = ref$pct_lo,
           censored_hi = ref$pct_hi, .before = 1)
}))
R1b_int <- purrr::map_dfr(which(MODIFIERS$set == "main"), function(k) {
  pr <- CK_PAIRS[[1]]
  m <- MODIFIERS$mod[k]
  d <- adni %>% filter(!is.na(.data[[m]]))
  d$.m1 <- as.integer(d[[m]] == 1)
  rhs <- sprintf("(z_l_ck + %s) * .m1", covs_without(m))
  coef_difference(d, pr, rhs, "z_l_ck:.m1") %>%
    mutate(analysis = paste("creatine kinase interaction with", MODIFIERS$label[k]),
           pair = pair_name(pr), var = "l_ck", .before = 1)
})
R1b <- bind_rows(R1b_main, R1b_int)
put_rev(R1b, "review_1b_c2n_difference_of_coefficients")
print(as.data.frame(R1b %>% transmute(analysis, pair, var, difference_of_coefficients = fmt(pct, lo, hi),
  censored_difference = fmt(censored_difference_pct, censored_lo, censored_hi), usable)), row.names = FALSE)

# =============================================================================
# 2. Calendar date of the draw
# =============================================================================
# Every sample was assayed in one campaign, so the date of the draw is also its
# storage time. With time since the first draw already in the model, the date
# is identified only between participants, by when they enrolled, so it indexes
# enrolment era and ADNI phase as much as storage, and its own coefficient is
# not read as degradation per year. The question is whether the factor
# estimates move. Every p-tau217 pair and the p-tau181 pair, for eGFR and
# creatine kinase, on the same rows as the main model.
say("2. The calendar date of the draw as a covariate")
adni$draw_year <- lubridate::decimal_date(adni$visit_date)
ALL_PAIRS <- c(combn(P217, 2, simplify = FALSE), list(P181))
R2 <- purrr::map_dfr(ALL_PAIRS, function(pr) purrr::map_dfr(c("egfr", "l_ck"), function(v) {
  zv <- paste0("z_", v)
  base <- check_pair(pr, v)
  adj  <- pair_fit(adni, pr, sprintf("%s + %s + draw_year", zv, COVS_ADNI), zv)
  yr   <- pair_fit(adni, pr, sprintf("%s + %s + draw_year", zv, COVS_ADNI), "draw_year")
  ok <- !is.na(adni[[zv]]) & !is.na(adni$draw_year)
  tibble(pair = pair_name(pr), var = v, method = base$method,
         n_obs = base$n_obs, n_obs_adjusted = adj$n_obs,
         main_pct = pct(base$b), main_lo = pct(base$lo), main_hi = pct(base$hi),
         adjusted_pct = pct(adj$b), adjusted_lo = pct(adj$lo), adjusted_hi = pct(adj$hi),
         per_year_earlier_draw_pct = pct(-yr$b), per_year_earlier_lo = pct(-yr$hi),
         per_year_earlier_hi = pct(-yr$lo),
         factor_draw_date_correlation = stats::cor(adni[[zv]][ok], adni$draw_year[ok]))
}))
if (any(R2$n_obs != R2$n_obs_adjusted)) stop("Adding the draw date changed the rows of a model.")
put_rev(R2, "review_2_draw_date")
print(as.data.frame(R2 %>% transmute(pair, var, main = fmt(main_pct, main_lo, main_hi),
  with_draw_date = fmt(adjusted_pct, adjusted_lo, adjusted_hi),
  per_year_earlier = fmt(per_year_earlier_draw_pct, per_year_earlier_lo, per_year_earlier_hi),
  r_factor_date = round(factor_draw_date_correlation, 2))), row.names = FALSE)
})
