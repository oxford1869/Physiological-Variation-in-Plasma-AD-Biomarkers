# =============================================================================
# 02_data_a4learn.R
# Stage 2. The A4 and LEARN analysis tables, built from the raw release.
# Writes paper1_out/a4_screen.rds, a4_long.rds, a4_lab.rds, the laboratory
# extract a4_lab_raw_v3.rds, and checks/a4_build.csv, a4_sources.csv.
# =============================================================================

# #############################################################################
# Section 2.1. The A4 and LEARN analysis tables
# Developed as 02_build_a4.R version 1.1.3, whole
# #############################################################################
local({
# =============================================================================
# 02_build_a4.R
# Paper 1: build the A4/LEARN analysis tables from the raw release.
#
# A4/LEARN carries the A4 half of research question 1 and the whole of question
# 4. It measures one platform per analyte, so it cannot speak to questions 2
# and 3, but it has two things ADNI does not: everyone screened, including the
# 1,394 who failed the amyloid PET screen, and repeated laboratory values.
#
# Builds three tables.
#   a4_screen  one row per person in the screening plasma panel (2,330), with
#              the Roche biomarkers on the log scale and the laboratory values
#              nearest that draw. This is the cross-section for question 1.
#   a4_long    one row per p-tau217 draw in A4 and LEARN, with the laboratory
#              values nearest that draw. This is the table for question 4.
#   a4_lab     every laboratory visit of the people in a4_long, which is what
#              the reliability calculations need.
#
# Reads    External Data/biomarker_Plasma_Roche_Results.csv  screening panel
#          External Data/biomarker_pTau217.csv               repeated p-tau217
#          External Data/clrm_lab.csv                        laboratory values
#          Derived Data/SUBJINFO.csv                         age, sex, APOE, BMI, Centiloid
#          All read only.
#
# Writes   paper1_out/a4_screen.rds, a4_long.rds, a4_lab.rds
#          paper1_out/a4_lab_raw_v2.rds   the laboratory extract, cached
#          paper1_out/checks/a4_build.csv, a4_sources.csv
#
# The laboratory file is 168 MB, so the first run reads it and saves the
# extract; later runs use the extract. Delete a4_lab_raw_v2.rds to force a reread.
# Dates in A4/LEARN are days since consent, not calendar dates.
# =============================================================================
# Version 1.1.3 renames the row "screen failures" to "screened, not enrolled in
# A4 or LEARN" in checks/a4_build.csv. No value changes.
SCRIPT_VERSION <- "1.1.3"
cat("02_build_a4.R version", SCRIPT_VERSION, "\n")

if (!exists("PAPER1_VERSION")) source("00_setup.R")

# A laboratory visit counts as concurrent with a blood draw within this many
# days. Screening draws and screening bloods are taken close together; the
# repeated draws are matched more tightly.
WINDOW_SCREEN_DAYS <- 90
WINDOW_DRAW_DAYS   <- 45

# ---- 1. People ---------------------------------------------------------------
say("Reading the A4/LEARN files")
si <- read_source(a4_file("subjinfo"), "SUBJINFO")

person <- si %>%
  transmute(BID,
            substudy  = SUBSTUDY,
            arm       = if_else(is.na(TX) | TX == "", "LEARN", TX),
            age0      = suppressWarnings(as.numeric(AGEYR)),
            female    = as.integer(SEX == "1"),
            # A missing genotype is unknown, not a non-carrier. Version 1.0.0
            # coalesced it to "" and so coded 23 screening people as e4 = 0.
            e4        = if_else(is.na(APOEGN) | str_trim(APOEGN) == "", NA_integer_,
                                as.integer(str_detect(APOEGN, "4"))),
            prev_bid  = if ("PREVBID" %in% names(si)) na_if(str_trim(PREVBID), "") else NA_character_,
            bmi       = suppressWarnings(as.numeric(BMIBL)),
            centiloid = suppressWarnings(as.numeric(AMYLCENT))) %>%
  filter(!is.na(BID))
say("SUBJINFO: %d people (%s)", nrow(person),
    paste(sprintf("%s %d", names(table(person$substudy)), table(person$substudy)), collapse = ", "))

# ---- 2. Laboratory values ----------------------------------------------------
# The raw file is large, so it is read once and the extract kept. The cache
# name carries a version: version 1.1.0 needs the text result and the result
# type (SIRESC, RPTRESC, RPTRTYP) for the below-limit CRP values, which the
# first cache, a4_lab_raw.rds, does not have, so that file is never read.
# Version 1.1.2 adds random glucose (RCT1669) to the panel, so it needs a new
# extract, a4_lab_raw_v3.rds; the v2 file lacks that test and is never read.
lab_cache <- file.path(OUT_DIR, "a4_lab_raw_v3.rds")
if (file.exists(lab_cache)) {
  lab_raw <- readRDS(lab_cache)
  say("  laboratory extract read from %s (%d rows). Delete that file to reread the raw file.",
      basename(lab_cache), nrow(lab_raw))
} else {
  say("  reading the laboratory file, which is large and takes a minute")
  lab_raw <- readr::read_csv(
    a4_file("lab"),
    col_select = any_of(c("SUBSTUDY", "BID", "VISCODE", "LBDTM_DAYS_CONSENT",
                          "LBTESTCD", "LBTEST", "SIRESN", "SIU",
                          "SIRESC", "RPTRESC", "RPTRTYP")),
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE, progress = FALSE)
  names(lab_raw) <- toupper(names(lab_raw))
  lab_raw <- lab_raw %>% filter(str_trim(LBTESTCD) %in% A4_LAB_CODES)
  saveRDS(lab_raw, lab_cache)
  say("  laboratory file read: %d rows kept for the %d tests in the panel",
      nrow(lab_raw), length(A4_LAB_CODES))
}

stopifnot(all(c("SIRESC", "RPTRESC", "RPTRTYP") %in% names(lab_raw)))
code_to_name <- setNames(names(A4_LAB_CODES), A4_LAB_CODES)
# CRP below the limit of the assay is reported as text, "<0.20" mg/L, with
# result type L and no number, and was dropped in version 1.0.0, which removed
# the low tail of CRP selectively (454 results). Such a result is now set to
# the limit it reports, parsed from the text result in SI units (SIRESC, or
# RPTRESC where SIRESC is blank), so that no draw loses its CRP for being low.
# Only CRP is treated this way; the few below-limit results of other tests are
# counted below but left out, as before.
lab_raw <- lab_raw %>%
  mutate(.value = suppressWarnings(as.numeric(SIRESN)),
         .below = is.na(.value) & coalesce(str_trim(RPTRTYP), "") == "L" &
                  str_detect(coalesce(SIRESC, RPTRESC, ""), "^\\s*<\\s*[0-9.]+"),
         .limit = if_else(.below, suppressWarnings(readr::parse_number(coalesce(SIRESC, RPTRESC))),
                          NA_real_))
below_limit <- lab_raw %>%
  filter(.below) %>%
  count(test = code_to_name[str_trim(LBTESTCD)], limit_text = coalesce(SIRESC, RPTRESC),
        name = "results")
say("  laboratory results reported below a limit:")
print(as.data.frame(below_limit), row.names = FALSE)
lab <- lab_raw %>%
  transmute(BID,
            day  = suppressWarnings(as.numeric(LBDTM_DAYS_CONSENT)),
            test = code_to_name[str_trim(LBTESTCD)],
            value = if_else(test == "crp" & .below, .limit, .value)) %>%
  filter(!is.na(BID), !is.na(day), !is.na(test), !is.na(value), value > 0) %>%
  group_by(BID, day, test) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = test, values_from = value) %>%
  left_join(person %>% select(BID, age0, female), by = "BID") %>%
  # Creatinine is in micromoles per litre here and the equation takes mg/dL.
  # The same plausibility limit on creatinine as in ADNI (eMethods 6), below
  # 20 mg/dL, applied after conversion from micromoles per litre.
  mutate(creatinine = if_else(!is.na(creatinine) & creatinine / UMOL_PER_MGDL < 20,
                              creatinine, NA_real_)) %>%
  mutate(age_at_lab = age0 + day / YEAR_DAYS,
         egfr = egfr_ckdepi_2021(creatinine / UMOL_PER_MGDL, age_at_lab,
                                 if_else(female == 1, "Female", "Male")),
         l_ck = safe_log(ck), l_ast = safe_log(ast), l_alt = safe_log(alt),
         l_crp = safe_log(crp)) %>%
  select(-age0, -female)
say("Laboratory: %d visits for %d people, median %d visits each",
    nrow(lab), n_distinct(lab$BID),
    median(lab %>% count(BID) %>% pull(n)))

# Laboratory values nearest a given day, one row per person and draw.
LAB_VARS <- c("egfr", "creatinine", "l_ck", "hb", "hba1c", "urea", "l_ast",
              "l_alt", "albumin", "l_crp", "glucose")
nearest_lab <- function(draws, window) {
  j <- draws %>%
    select(BID, draw_day) %>%
    inner_join(lab, by = "BID", relationship = "many-to-many") %>%
    mutate(gap = abs(day - draw_day)) %>%
    filter(gap <= window)
  # Each factor is taken from its own nearest visit, because a visit may have
  # some of the panel and not the rest.
  out <- lapply(intersect(LAB_VARS, names(j)), function(v) {
    j %>% filter(!is.na(.data[[v]])) %>%
      group_by(BID, draw_day) %>%
      slice_min(gap, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(BID, draw_day, var = v, val = .data[[v]])
  })
  bind_rows(out) %>% pivot_wider(names_from = var, values_from = val)
}

# ---- 3. The screening cross-section ------------------------------------------
ro <- read_source(a4_file("roche"), "Roche panel")

ROCHE <- c(pTau181 = "TPP181", NfL = "NF-L", GFAP = "GFAP",
           Abeta40 = "AMYLB40", Abeta42 = "AMYLB42")
# Results below the limit of quantification are reported as "BLQ" with no
# number. They are left out, as before, but counted per analyte so that the
# exclusion can be disclosed (checks/a4_build.csv).
roche_blq_rows <- ro %>%
  transmute(BID, code = str_trim(LBTESTCD), blq = str_trim(coalesce(LABRESC, "")) == "BLQ") %>%
  filter(code %in% ROCHE, blq) %>%
  mutate(analyte = names(ROCHE)[match(code, ROCHE)])
roche <- ro %>%
  transmute(BID,
            code  = str_trim(LBTESTCD),
            unit  = str_to_lower(str_trim(coalesce(LABORESU, ""))),
            value = suppressWarnings(as.numeric(LABRESN)),
            draw_day = suppressWarnings(as.numeric(LABD_DAYS_CONSENT))) %>%
  filter(code %in% ROCHE, !is.na(value), value > 0) %>%
  # ng/mL to pg/mL, so that the amyloid ratio is a ratio of like with like
  mutate(value = if_else(unit == "ng/ml", value * 1000, value),
         analyte = names(ROCHE)[match(code, ROCHE)]) %>%
  group_by(BID, analyte) %>%
  summarise(value = mean(value), draw_day = first(draw_day), .groups = "drop")

screen_day <- roche %>% group_by(BID) %>% summarise(draw_day = first(draw_day), .groups = "drop")
screen_wide <- roche %>%
  mutate(lv = log(value)) %>%
  select(BID, analyte, lv) %>%
  pivot_wider(names_from = analyte, values_from = lv) %>%
  mutate(Ab42_40 = Abeta42 - Abeta40)

a4_screen <- screen_day %>%
  left_join(screen_wide, by = "BID") %>%
  left_join(nearest_lab(screen_day, WINDOW_SCREEN_DAYS), by = c("BID", "draw_day")) %>%
  inner_join(person, by = "BID") %>%
  mutate(age_at_visit = age0 + draw_day / YEAR_DAYS,
         cohort = "A4 screening") %>%
  arrange(BID)

# A person who failed screening and was screened again has a new BID, linked
# to the earlier one through PREVBID in SUBJINFO. Where both BIDs have a
# screening panel the person would be counted twice (15 people in this
# release), so the earlier record is dropped and the most recent screening
# record, the one that determined the person's final status, is kept. No chain
# of more than two BIDs occurs, and the check stops if one ever does.
stopifnot(!any(person$prev_bid[!is.na(person$prev_bid)] %in%
               person$BID[!is.na(person$prev_bid)]))
rescreen_earlier <- a4_screen$prev_bid[!is.na(a4_screen$prev_bid) &
                                         a4_screen$prev_bid %in% a4_screen$BID]
a4_screen <- a4_screen %>% filter(!BID %in% rescreen_earlier) %>% select(-prev_bid)
say("Screening panel: %d earlier records of re-screened people dropped, %d people kept",
    length(rescreen_earlier), nrow(a4_screen))
stopifnot(!anyDuplicated(a4_screen$BID))
roche_blq <- roche_blq_rows %>%
  filter(BID %in% a4_screen$BID) %>%
  count(analyte, name = "blq")
say("Roche results below the limit of quantification, left out, in the screening sample: %s",
    paste(sprintf("%s %d", roche_blq$analyte, roche_blq$blq), collapse = ", "))

a4_screen_factors <- FACTORS %>% filter(var %in% names(a4_screen))
for (v in a4_screen_factors$var)
  a4_screen[[paste0("z_", v)]] <- FACTOR_SIGN[[v]] * zf(a4_screen[[v]])

# ---- 4. The repeated p-tau217 draws ------------------------------------------
pt <- read_source(a4_file("ptau217"), "p-tau217")

# The assay reports "<LLOQ" or ">ULOQ" rather than a number for draws outside
# its quantification range, but the reading itself is in ORRESRAW, and for
# every quantified draw ORRESRAW equals the reported value exactly. The reading
# is therefore used throughout and the out-of-range draws are flagged instead
# of dropped: dropping them would select on the outcome, and they are 15% of
# the series.
# Three draws above the ULOQ have no reading ("> 2x ULOQ not reported"). They
# cannot be given a value and are left out, but they are counted and printed,
# because they are the highest values in the series.
uloq_no_reading <- pt %>%
  filter(SUBSTUDY %in% c("A4", "LEARN"),
         str_detect(coalesce(ORRES, ""), "ULOQ"),
         is.na(suppressWarnings(as.numeric(ORRESRAW))))
say("p-tau217: %d draws above the ULOQ have no reading and are left out (%d people)",
    nrow(uloq_no_reading), n_distinct(uloq_no_reading$BID))
draws <- pt %>%
  filter(SUBSTUDY %in% c("A4", "LEARN")) %>%
  transmute(BID,
            draw_day   = suppressWarnings(as.numeric(COLLECTION_DATE_DAYS_CONSENT)),
            reported   = suppressWarnings(as.numeric(ORRES)),
            reading    = suppressWarnings(as.numeric(ORRESRAW)),
            below_lloq = str_detect(coalesce(ORRES, ""), "LLOQ"),
            above_uloq = str_detect(coalesce(ORRES, ""), "ULOQ")) %>%
  mutate(value = coalesce(reading, reported)) %>%
  filter(!is.na(BID), !is.na(draw_day), !is.na(value), value > 0) %>%
  group_by(BID, draw_day) %>%
  summarise(value = mean(value), quantified = all(!below_lloq & !above_uloq),
            below_lloq = any(below_lloq), above_uloq = any(above_uloq),
            .groups = "drop")

a4_long <- draws %>%
  left_join(nearest_lab(draws, WINDOW_DRAW_DAYS), by = c("BID", "draw_day")) %>%
  inner_join(person, by = "BID") %>%
  group_by(BID) %>%
  mutate(first_draw = min(draw_day)) %>%
  ungroup() %>%
  mutate(pTau217      = log(value),
         plasma_time  = (draw_day - first_draw) / YEAR_DAYS,
         age_at_visit = age0 + draw_day / YEAR_DAYS,
         cohort       = "A4 longitudinal") %>%
  arrange(BID, draw_day)

# The same re-screen rule for the repeated draws. In this release no earlier
# (screen-failure) BID has a p-tau217 draw, so nothing is dropped, and the check
# stops rather than letting one person enter twice if a later release differs.
long_rescreen_pairs <- sum(!is.na(a4_long$prev_bid) & a4_long$prev_bid %in% a4_long$BID)
stopifnot(long_rescreen_pairs == 0)
a4_long <- a4_long %>% select(-prev_bid)

a4_long_factors <- FACTORS %>% filter(var %in% names(a4_long))
for (v in a4_long_factors$var)
  a4_long[[paste0("z_", v)]] <- FACTOR_SIGN[[v]] * zf(a4_long[[v]])

# Every laboratory visit of these people, in years since their first draw, for
# the reliability calculations in question 4.
a4_lab <- lab %>%
  semi_join(a4_long %>% distinct(BID), by = "BID") %>%
  left_join(a4_long %>% distinct(BID, first_draw), by = "BID") %>%
  mutate(lab_time = (day - first_draw) / YEAR_DAYS) %>%
  arrange(BID, day)

save_p1(a4_screen, "a4_screen")
save_p1(a4_long,   "a4_long")
save_p1(a4_lab,    "a4_lab")

# ---- 5. What was built -------------------------------------------------------
n_with <- function(d, x) sum(!is.na(d[[x]]))
chk <- tibble::tribble(
  ~item, ~value,
  "Screening panel: people", nrow(a4_screen),
  "  earlier records of re-screened people dropped", length(rescreen_earlier),
  "  A4 randomised", sum(a4_screen$substudy == "A4"),
  "  LEARN", sum(a4_screen$substudy == "LEARN"),
  "  screened, not enrolled in A4 or LEARN", sum(a4_screen$substudy == "SF"),
  "  with p-tau181", n_with(a4_screen, "pTau181"),
  "  with eGFR", n_with(a4_screen, "egfr"),
  "  with creatine kinase", n_with(a4_screen, "l_ck"),
  "  with Centiloid", n_with(a4_screen, "centiloid"),
  "  with APOE genotype", n_with(a4_screen, "e4"),
  "  with CRP", n_with(a4_screen, "l_crp"),
  "  Roche Abeta40 below the limit, left out", sum(roche_blq$blq[roche_blq$analyte == "Abeta40"]),
  "  Roche Abeta42 below the limit, left out", sum(roche_blq$blq[roche_blq$analyte == "Abeta42"]),
  "  Roche GFAP below the limit, left out", sum(roche_blq$blq[roche_blq$analyte == "GFAP"]),
  "  Roche NfL below the limit, left out", sum(roche_blq$blq[roche_blq$analyte == "NfL"]),
  "  Roche p-tau181 below the limit, left out", sum(roche_blq$blq[roche_blq$analyte == "pTau181"]),
  "  Abeta42/40 ratio missing", sum(is.na(a4_screen$Ab42_40)),
  "Repeated p-tau217: draws", nrow(a4_long),
  "  people", n_distinct(a4_long$BID),
  "  draws below the limit of quantification", sum(a4_long$below_lloq),
  "  draws above it", sum(a4_long$above_uloq),
  "  draws above it with no reading, left out", nrow(uloq_no_reading),
  "  people with APOE genotype", n_distinct(a4_long$BID[!is.na(a4_long$e4)]),
  "  people with at least one draw below the limit", n_distinct(a4_long$BID[a4_long$below_lloq]),
  "  draws with laboratory values within 45 days", n_with(a4_long, "egfr"),
  "  draws with CRP within 45 days", n_with(a4_long, "l_crp"),
  "CRP results below the limit, set to the limit", sum(below_limit$results[below_limit$test == "crp"]),
  "Laboratory visits of those people", nrow(a4_lab),
  "  people with three or more creatinine values",
    sum((a4_lab %>% filter(!is.na(egfr)) %>% count(BID) %>% pull(n)) >= 3)
)
print(as.data.frame(chk), row.names = FALSE)
write_csv(chk, file.path(CHECK_DIR, "a4_build.csv"))
write_csv(source_table(), file.path(CHECK_DIR, "a4_sources.csv"))

say("Follow-up of the repeated draws: median %.2f years, longest %.2f",
    median(a4_long %>% group_by(BID) %>% summarise(s = max(plasma_time)) %>% pull(s)),
    max(a4_long$plasma_time))
say("a4_screen: %d rows. a4_long: %d rows. a4_lab: %d rows.",
    nrow(a4_screen), nrow(a4_long), nrow(a4_lab))
})
