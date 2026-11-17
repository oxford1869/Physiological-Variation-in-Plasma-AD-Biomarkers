# =============================================================================
# 01_data_adni.R
# Stage 1. The ADNI analysis table, built from the raw release.
# Writes paper1_out/adni_visit.rds and adni_long.rds, and checks/adni_*.csv.
# =============================================================================

# #############################################################################
# Section 1.1. The ADNI analysis table
# Developed as 01_build_adni.R version 1.2.1, whole
# #############################################################################
local({
# =============================================================================
# 01_build_adni.R
# Paper 1: build the ADNI analysis table from the raw release.
#
# ADNI measured the same aliquot on up to six platforms, so it carries research
# questions 2 and 3 and the ADNI half of question 1.
#
# Sample   the 393 people with at least one plasma sample measured by a platform
#          other than Quanterix. Samples that only Quanterix measured are dropped,
#          because nothing can be compared across platforms at those.
#
# Reads    FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.csv  plasma, every platform
#          ADNIMERGE_DATA.csv                       age, sex, APOE4, diagnosis
#          UCBERKELEY_AMY_6MM.csv                   amyloid PET, Centiloid
#          LABDATA.csv                              creatinine and the screening panel
#          VITALS.csv                               weight and height for BMI
#          RECCMEDS.csv                             statin use
#          All read only.
#
# Writes   paper1_out/adni_visit.rds  one row per person and plasma sample, with
#                                     one column per measure on the log scale
#          paper1_out/adni_long.rds   the same in the shared long shape, one row
#                                     per person, sample and measure
#          paper1_out/checks/adni_build.csv
#
# Takes about half a minute.
# =============================================================================
SCRIPT_VERSION <- "1.2.1"
cat("01_build_adni.R version", SCRIPT_VERSION, "\n")

if (!exists("PAPER1_VERSION")) source("00_setup.R")

# Amyloid PET counts as concurrent within this many days of the blood draw, and
# weight within a year. Both windows follow the papers 2 and 3 code.
WINDOW_AMYLOID_DAYS <- 180
WINDOW_BMI_DAYS     <- 365

# Which platform and analyte each measure in the raw file belongs to. role is
# "measure" for anything analysed in its own right and "component" for the
# parts of a ratio.
MEASURES <- tibble::tribble(
  ~measure,                      ~analyte,   ~platform,   ~role,
  "C2N_plasma_ptau217",          "pTau217",  "C2N",       "measure",
  "Fuji_plasma_ptau217",         "pTau217",  "Fujirebio", "measure",
  "AlzPath_plasma_ptau217",      "pTau217",  "ALZpath",   "measure",
  "Janssen_plasma_ptau217",      "pTau217",  "Janssen",   "measure",
  "C2N_plasma_nptau217",         "npTau217", "C2N",       "measure",
  "C2N_plasma_ptau217_ratio",    "pTau217 ratio", "C2N",  "measure",
  "Roche_plasma_ptau181",        "pTau181",  "Roche",     "measure",
  "QX_plasma_ptau181",           "pTau181",  "Quanterix", "measure",
  "Roche_plasma_NfL",            "NfL",      "Roche",     "measure",
  "QX_plasma_NfL",               "NfL",      "Quanterix", "measure",
  "Roche_plasma_GFAP",           "GFAP",     "Roche",     "measure",
  "QX_plasma_GFAP",              "GFAP",     "Quanterix", "measure",
  "C2N_plasma_Abeta40",          "Abeta40",  "C2N",       "measure",
  "Fuji_plasma_Ab40",            "Abeta40",  "Fujirebio", "measure",
  "Roche_plasma_Ab40",           "Abeta40",  "Roche",     "measure",
  "QX_plasma_Ab40",              "Abeta40",  "Quanterix", "measure",
  "C2N_plasma_Abeta42",          "Abeta42",  "C2N",       "measure",
  "Fuji_plasma_Ab42",            "Abeta42",  "Fujirebio", "measure",
  "Roche_plasma_Ab42",           "Abeta42",  "Roche",     "measure",
  "QX_plasma_Ab42",              "Abeta42",  "Quanterix", "measure",
  "C2N_plasma_Abeta42_Abeta40",  "Ab42_40",  "C2N",       "measure"
)
# Ratios the laboratory does not supply, derived below as the difference of the
# two logged components.
RATIOS <- tibble::tribble(
  ~measure,                  ~numerator,          ~denominator,        ~platform,
  "Fuji_plasma_Ab42_Ab40",   "Fuji_plasma_Ab42",  "Fuji_plasma_Ab40",  "Fujirebio",
  "Roche_plasma_Ab42_Ab40",  "Roche_plasma_Ab42", "Roche_plasma_Ab40", "Roche",
  "QX_plasma_Ab42_Ab40",     "QX_plasma_Ab42",    "QX_plasma_Ab40",    "Quanterix"
)

# ---- 1. Plasma: one row per person, sample and measure -----------------------
say("Reading the raw files")
raw <- read_source(adni_file("plasma"), "plasma")

pl <- raw %>%
  transmute(
    RID       = as.integer(RID),
    measure   = PLASMA_BIOMARKER,
    exam_date = parse_date_dmy(EXAMDATE),
    value_raw = TESTVALUE,
    units     = UNITS,
    comments  = COMMENTS
  ) %>%
  filter(!is.na(RID), !is.na(exam_date)) %>%
  mutate(
    # A value may carry a comparator such as "<1.3". Keep that, and anything the
    # comments say about the detection limit, before coercing to a number.
    comparator = str_trim(str_extract(value_raw, "^\\s*[<>]")),
    value      = if_else(str_detect(coalesce(value_raw, ""), "[0-9]"),
                         suppressWarnings(readr::parse_number(value_raw)), NA_real_),
    at_floor   = (!is.na(comparator) & comparator == "<") |
                 str_detect(str_to_lower(coalesce(comments, "")), "lod|lloq|below|detection"),
    # ng/mL to pg/mL so that every concentration is in the same unit
    value      = if_else(str_to_lower(str_trim(coalesce(units, ""))) == "ng/ml",
                         value * 1000, value)
  ) %>%
  filter(!is.na(value), value > 0)

# A sample is one blood draw. Dates within a week of each other for the same
# person are the same draw.
visits <- pl %>%
  distinct(RID, exam_date) %>%
  arrange(RID, exam_date) %>%
  group_by(RID) %>%
  mutate(new = is.na(lag(exam_date)) | as.integer(exam_date - lag(exam_date)) > VISIT_GAP_DAYS,
         visit_no = cumsum(new)) %>%
  group_by(RID, visit_no) %>%
  mutate(visit_date = min(exam_date)) %>%
  ungroup() %>%
  select(RID, exam_date, visit_no, visit_date)

pl <- pl %>%
  left_join(visits, by = c("RID", "exam_date")) %>%
  left_join(MEASURES, by = "measure") %>%
  mutate(platform = coalesce(platform, "unmapped"))

unmapped <- pl %>% filter(platform == "unmapped") %>% distinct(measure)
if (nrow(unmapped) > 0) {
  say("WARNING: %d measures in the raw file are not in the MEASURES table:", nrow(unmapped))
  print(as.data.frame(unmapped))
}

# ---- 2. The sample -----------------------------------------------------------
# Quanterix measured a partly different set of samples from the other five
# laboratories. A person is in the sample if at least one of their draws was
# measured by another platform, and only those draws are kept.
draw_platforms <- pl %>%
  filter(platform != "unmapped") %>%
  group_by(RID, visit_no) %>%
  summarise(n_platforms = n_distinct(platform),
            quanterix_only = all(platform == "Quanterix"), .groups = "drop")
keep_people <- draw_platforms %>% filter(!quanterix_only) %>% distinct(RID)

# Time is counted from a person's first blood draw of any kind, including one
# that only Quanterix measured, because that is when their sampling began.
person_first <- pl %>%
  semi_join(keep_people, by = "RID") %>%
  group_by(RID) %>%
  summarise(first_plasma = min(visit_date), .groups = "drop")

pl <- pl %>%
  semi_join(keep_people, by = "RID") %>%
  left_join(draw_platforms, by = c("RID", "visit_no")) %>%
  filter(!quanterix_only)

sample_visits <- pl %>%
  distinct(RID, visit_no, visit_date, n_platforms) %>%
  left_join(person_first, by = "RID") %>%
  mutate(plasma_time = as.numeric(visit_date - first_plasma) / YEAR_DAYS)
say("Sample: %d people, %d blood draws", n_distinct(sample_visits$RID), nrow(sample_visits))

# ---- 3. Wide on the log scale, and the derived ratios ------------------------
wide <- pl %>%
  filter(platform != "unmapped") %>%
  group_by(RID, visit_no, measure) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  mutate(lv = log(value)) %>%
  select(-value) %>%
  pivot_wider(names_from = measure, values_from = lv)

for (i in seq_len(nrow(RATIOS))) {
  num <- RATIOS$numerator[i]; den <- RATIOS$denominator[i]
  if (all(c(num, den) %in% names(wide))) wide[[RATIOS$measure[i]]] <- wide[[num]] - wide[[den]]
}
MEASURES <- bind_rows(MEASURES,
                      RATIOS %>% transmute(measure, analyte = "Ab42_40", platform, role = "measure"))

# Which C2N p-tau217 values sit at the detection limit. The laboratory reports
# half the limit for these, so they are flagged rather than trusted.
floor_flag <- pl %>%
  filter(measure == "C2N_plasma_ptau217") %>%
  group_by(RID, visit_no) %>%
  summarise(c2n_at_floor = any(at_floor) | any(abs(value - C2N_REPORTED) < 1e-8),
            .groups = "drop")

# ---- 4. Age, sex, APOE4 and diagnosis from ADNIMERGE -------------------------
mg_raw <- read_source(adni_file("merge"), "ADNIMERGE")
m_vis  <- pick_col(mg_raw, c("VISCODE", "VISCODE2"))
m_date <- pick_col(mg_raw, c("EXAMDATE"))
m_dx   <- pick_col(mg_raw, c("DX"), required = FALSE)

mg <- mg_raw %>%
  mutate(RID = as.integer(RID),
         exam_date = parse_date_ymd(.data[[m_date]]),
         viscode   = recode_viscode(.data[[m_vis]])) %>%
  semi_join(keep_people, by = "RID")

# Baseline is the row coded bl, or the earliest dated row if there is none.
person <- mg %>%
  arrange(RID, viscode != "bl", exam_date) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  transmute(RID,
            adni_bl_date = exam_date,
            age_bl       = suppressWarnings(as.numeric(AGE)),
            sex          = PTGENDER,
            educ         = suppressWarnings(as.numeric(PTEDUCAT)),
            apoe4_n      = suppressWarnings(as.integer(APOE4))) %>%
  mutate(female = as.integer(sex == "Female"),
         e4     = as.integer(apoe4_n >= 1))

# Diagnosis at the ADNIMERGE visit nearest the first blood draw.
dx_first <- NULL
if (!is.na(m_dx)) {
  dx_first <- mg %>%
    filter(!is.na(.data[[m_dx]]), .data[[m_dx]] != "") %>%
    transmute(RID, dx_date = exam_date, dx_raw = .data[[m_dx]]) %>%
    inner_join(sample_visits %>% distinct(RID, first_plasma), by = "RID") %>%
    mutate(gap = as.integer(dx_date - first_plasma)) %>%
    group_by(RID) %>% slice_min(abs(gap), n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(RID, dx = dx_raw, dx_gap_days = gap)
  person <- person %>% left_join(dx_first, by = "RID")
}

# ---- 5. Amyloid PET ----------------------------------------------------------
amy_raw <- read_source(adni_file("amyloid"), "amyloid")
a_date  <- pick_col(amy_raw, c("SCANDATE", "EXAMDATE"))
a_cl    <- pick_col(amy_raw, c("CENTILOIDS", "CENTILOID"))

amy <- amy_raw %>%
  transmute(RID = as.integer(RID),
            scan_date = parse_date_ymd(.data[[a_date]]),
            centiloid = suppressWarnings(as.numeric(.data[[a_cl]]))) %>%
  filter(!is.na(RID), !is.na(scan_date), !is.na(centiloid)) %>%
  semi_join(keep_people, by = "RID") %>%
  # Each draw takes its own nearest scan within the window (see nearest_visit).
  nearest_visit("scan_date", sample_visits %>% select(RID, visit_no, visit_date),
                window = WINDOW_AMYLOID_DAYS) %>%
  select(RID, visit_no, centiloid, centiloid_gap_days = gap_days)

# ---- 6. Weight and height for BMI --------------------------------------------
vt_raw <- read_source(adni_file("vitals"), "vitals")
v_date <- pick_col(vt_raw, c("EXAMDATE", "VISDATE"))
num <- function(x) { x <- suppressWarnings(as.numeric(x)); if_else(!is.na(x) & x > 0, x, NA_real_) }

vt <- vt_raw %>%
  transmute(RID = as.integer(RID),
            exam_date = parse_date_ymd(.data[[v_date]]),
            wt  = num(VSWEIGHT),
            wtu = if ("VSWTUNIT" %in% names(vt_raw)) VSWTUNIT else NA_character_,
            ht  = if ("VSHEIGHT" %in% names(vt_raw)) num(VSHEIGHT) else NA_real_,
            htu = if ("VSHTUNIT" %in% names(vt_raw)) VSHTUNIT else NA_character_) %>%
  filter(!is.na(RID), !is.na(exam_date)) %>%
  semi_join(keep_people, by = "RID") %>%
  # ADNI unit codes: weight 1 pounds, 2 kilograms; height 1 inches, 2 centimetres
  mutate(weight_kg = case_when(wtu == "1" ~ wt * 0.45359237,
                               wtu == "2" ~ wt,
                               is.na(wtu) & wt > 130 ~ wt * 0.45359237,
                               TRUE ~ wt),
         height_cm = case_when(htu == "1" ~ ht * 2.54,
                               htu == "2" ~ ht,
                               is.na(htu) & ht < 100 ~ ht * 2.54,
                               TRUE ~ ht),
         weight_kg = if_else(weight_kg >= 30  & weight_kg <= 250, weight_kg, NA_real_),
         height_cm = if_else(height_cm >= 120 & height_cm <= 220, height_cm, NA_real_)) %>%
  filter(!is.na(weight_kg) | !is.na(height_cm)) %>%
  distinct(RID, exam_date, .keep_all = TRUE) %>%
  # Height is measured at screening and often missing later, so a person's own
  # median measured height fills the visits that have none.
  group_by(RID) %>%
  mutate(height_used = coalesce(height_cm, median(height_cm, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(bmi = weight_kg / (height_used / 100)^2) %>%
  filter(!is.na(bmi)) %>%
  # Each draw takes its own nearest vitals visit within the window.
  nearest_visit("exam_date", sample_visits %>% select(RID, visit_no, visit_date),
                window = WINDOW_BMI_DAYS) %>%
  select(RID, visit_no, bmi, weight_kg, bmi_gap_days = gap_days)

# ---- 7. Creatinine, eGFR and the screening laboratory panel ------------------
lb_raw <- read_source(adni_file("lab"), "laboratory")
l_date <- pick_col(lb_raw, c("EXAMDATE", "VISDATE"))
l_vis  <- pick_col(lb_raw, c("VISCODE2", "VISCODE"))
l_vis1 <- pick_col(lb_raw, c("VISCODE"), required = FALSE)

# Creatinine is measured at many visits. The value used is the one closest to
# the first blood draw, and its gap is kept so that the timing can be checked.
creat <- lb_raw %>%
  transmute(RID = as.integer(RID),
            exam_date = parse_date_ymd(.data[[l_date]]),
            cr_raw = RCT392) %>%
  filter(!is.na(RID), !is.na(exam_date)) %>%
  semi_join(keep_people, by = "RID") %>%
  mutate(creatinine = if_else(str_detect(coalesce(cr_raw, ""), "^\\s*[0-9]*\\.?[0-9]+\\s*$"),
                              suppressWarnings(as.numeric(cr_raw)), NA_real_),
         creatinine = if_else(creatinine > 0 & creatinine < 20, creatinine, NA_real_)) %>%
  filter(!is.na(creatinine)) %>%
  distinct(RID, exam_date, .keep_all = TRUE) %>%
  left_join(person %>% select(RID, age_bl, adni_bl_date, sex), by = "RID") %>%
  left_join(sample_visits %>% distinct(RID, first_plasma), by = "RID") %>%
  mutate(age_at_lab = age_bl + as.numeric(exam_date - adni_bl_date) / YEAR_DAYS,
         egfr       = egfr_ckdepi_2021(creatinine, age_at_lab, sex),
         egfr_2009  = egfr_ckdepi_2009(creatinine, age_at_lab, sex),
         egfr_gap_years = as.numeric(exam_date - first_plasma) / YEAR_DAYS) %>%
  group_by(RID) %>%
  slice_min(abs(egfr_gap_years), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(RID, creatinine, egfr, egfr_2009, egfr_gap_years, egfr_date = exam_date)

# The rest of the panel is measured once, at screening.
# The visit is read from VISCODE2 and, where VISCODE2 is blank, from VISCODE:
# some screening rows (RIDs 4548 and 4835) carry VISCODE "sc" with no VISCODE2,
# and testing VISCODE2 alone lost their whole panel while their creatinine,
# from the same row, was used.
# A screening result can be unusable, for example a haemolysed sample coded
# SCC368 in place of the haemoglobin, and ADNI then re-draws at an unscheduled
# visit coded "nv" (RIDs 123 and 4386, 6 and 11 days after screening). An "nv"
# row is therefore used, but only to fill a test that has no valid screening
# value and only if it is dated on or after the screening row, so that it
# replaces an unusable result and never a valid one. lab_date stays the
# screening date.
lab_vis <- if (!is.na(l_vis1) && l_vis1 != l_vis)
  coalesce(na_if(str_trim(lb_raw[[l_vis]]), ""), str_trim(lb_raw[[l_vis1]])) else
  str_trim(lb_raw[[l_vis]])
panel_rows <- lb_raw %>%
  mutate(.vis = lab_vis, .file_order = row_number()) %>%
  filter(.vis %in% c("sc", "scmri", "f", "nv")) %>%
  transmute(RID = as.integer(RID),
            lab_date = parse_date_ymd(.data[[l_date]]),
            redraw = .vis == "nv",
            .file_order,
            across(all_of(unname(ADNI_LAB_CODES)), ~ suppressWarnings(as.numeric(.x)))) %>%
  rename(!!!ADNI_LAB_CODES) %>%
  semi_join(keep_people, by = "RID") %>%
  mutate(across(all_of(names(ADNI_LAB_CODES)), ~ if_else(!is.na(.x) & .x > 0, .x, NA_real_))) %>%
  # ADNI reports haemoglobin in g/dL and A4/LEARN in g/L. The analyses are on
  # z-scores and so do not care, but a descriptive table does, and a raw
  # coefficient would. Everything is put into g/L here.
  mutate(hb = hb * 10) %>%
  group_by(RID) %>%
  mutate(sc_date = suppressWarnings(min(lab_date[!redraw], na.rm = TRUE))) %>%
  ungroup() %>%
  # A re-draw is kept only after a screening row of the same person.
  filter(!redraw | (is.finite(sc_date) & !is.na(lab_date) & lab_date >= sc_date))
# Which values came from a re-draw, for the checks file.
panel_redraw <- panel_rows %>%
  group_by(RID) %>%
  summarise(across(all_of(names(ADNI_LAB_CODES)),
                   ~ all(is.na(.x[!redraw])) & any(!is.na(.x[redraw]))), .groups = "drop") %>%
  tidyr::pivot_longer(-RID, names_to = "test", values_to = "from_redraw") %>%
  filter(from_redraw)
panel <- panel_rows %>%
  # Screening rows first, in file order, then the re-draws by date.
  arrange(RID, redraw, if_else(redraw, as.numeric(lab_date), as.numeric(.file_order))) %>%
  mutate(lab_date = if_else(redraw, as.Date(NA), lab_date)) %>%
  select(-redraw, -.file_order, -sc_date) %>%
  group_by(RID) %>%
  summarise(across(everything(), ~ first(na.omit(.x))), .groups = "drop") %>%
  left_join(sample_visits %>% distinct(RID, first_plasma), by = "RID") %>%
  mutate(lab_gap_years = as.numeric(lab_date - first_plasma) / YEAR_DAYS) %>%
  select(-first_plasma)

# ---- 8. Statin use -----------------------------------------------------------
# "statin" on its own also matches the antifungal nystatin and its brand and
# misspelt forms (mycostatin, nistatin), so those are excluded explicitly. The
# misspellings "rosuvastaton", "simvaststin" and "rasuvastine" occur in
# RECCMEDS for people in the sample and are added. Every string the
# pattern matches is written to checks/adni_statin_strings.csv so that the
# list can be inspected.
STATIN_PATTERN <- paste0("statin|atorvastatin|simvastatin|rosuvastatin|rosuvastaton|simvaststin|rasuvastine|pravastatin|",
                         "lovastatin|fluvastatin|pitavastatin|cerivastatin|",
                         "lipitor|zocor|crestor|pravachol|mevacor|altoprev|lescol|",
                         "livalo|vytorin|caduet|advicor|simcor")
NOT_STATIN_PATTERN <- "nystatin|mycostatin|nistatin"
med_raw <- read_source(adni_file("meds"), "medication")
m_mdate <- pick_col(med_raw, c("VISDATE", "EXAMDATE", "USERDATE"))
med <- med_raw %>%
  transmute(RID = as.integer(RID),
            med_date = parse_date_ymd(.data[[m_mdate]]),
            med = str_to_lower(coalesce(CMMED, ""))) %>%
  filter(!is.na(RID)) %>%
  semi_join(keep_people, by = "RID") %>%
  mutate(is_statin = str_detect(med, STATIN_PATTERN) & !str_detect(med, NOT_STATIN_PATTERN))
statin_strings <- med %>%
  filter(is_statin) %>%
  count(med, name = "records") %>%
  arrange(desc(records), med)
write_csv(statin_strings, file.path(CHECK_DIR, "adni_statin_strings.csv"))
say("Medication strings counted as a statin: %d distinct (checks/adni_statin_strings.csv)",
    nrow(statin_strings))
print(as.data.frame(statin_strings), row.names = FALSE)
# Statin use at each draw, on the same rule as the first-draw indicator: a
# statin record dated on or before that draw's date (or undated) counts.
statin_draw <- med %>%
  filter(is_statin) %>%
  inner_join(sample_visits %>% select(RID, visit_no, visit_date), by = "RID",
             relationship = "many-to-many") %>%
  group_by(RID, visit_no) %>%
  summarise(statin_at_draw = as.integer(any(is.na(med_date) | med_date <= visit_date)),
            .groups = "drop")
# ---- 8b. Medication acting on the kidney -----------------------------------
# Drugs that change glomerular filtration or the handling of creatinine, and so
# could move eGFR, creatine kinase clearance or a platform difference through
# the kidney: loop, thiazide and potassium-sparing diuretics, ACE inhibitors,
# angiotensin receptor blockers, angiotensin receptor-neprilysin inhibitors and
# SGLT2 inhibitors. Generic and common brand
# names are listed explicitly rather than matched by suffix, so that nothing
# unintended is caught; every matched string is written to
# checks/adni_renal_medication_strings.csv for inspection.
RENAL_PATTERN <- paste0(
  "furosemide|frusemide|lasix|bumetanide|bumex|torsemide|torasemide|demadex|",
  "hydrochlorothiazide|hctz|chlorthalidone|chlortalidone|indapamide|metolazone|chlorothiazide|",
  "spironolactone|aldactone|eplerenone|inspra|amiloride|triamterene|dyazide|maxzide|",
  "lisinopril|enalapril|ramipril|benazepril|quinapril|fosinopril|captopril|perindopril|",
  "trandolapril|moexipril|zestril|prinivil|zestoretic|vasotec|altace|lotensin|lotrel|accupril|",
  "losartan|valsartan|irbesartan|candesartan|olmesartan|telmisartan|azilsartan|eprosartan|",
  "cozaar|hyzaar|diovan|exforge|avapro|avalide|benicar|azor|micardis|atacand|edarbi|",
  "empagliflozin|dapagliflozin|canagliflozin|ertugliflozin|jardiance|farxiga|invokana|invokamet|",
  # angiotensin receptor-neprilysin inhibitor
  "entresto|sacubitril|",
  # further brand names and the misspellings found in RECCMEDS
  "monopril|aceon|microzide|temisartan|quinopril|linsinopril|",
  "thiazide|thizide|thazide|triazide|athyzide|hydrochlorathiazide|hydrochlorothyazide|hydrocholoro")
med <- med %>% mutate(is_renal = str_detect(med, RENAL_PATTERN))
renal_strings <- med %>% filter(is_renal) %>% count(med, name = "records") %>%
  arrange(desc(records), med)
write_csv(renal_strings, file.path(CHECK_DIR, "adni_renal_medication_strings.csv"))
say("Medication strings counted as acting on the kidney: %d distinct (checks/adni_renal_medication_strings.csv)",
    nrow(renal_strings))
# Use at each draw on the same rule as statins: a record dated on or before the
# draw, or undated, counts.
renal_draw <- med %>%
  filter(is_renal) %>%
  inner_join(sample_visits %>% select(RID, visit_no, visit_date), by = "RID",
             relationship = "many-to-many") %>%
  group_by(RID, visit_no) %>%
  summarise(renal_med_at_draw = as.integer(any(is.na(med_date) | med_date <= visit_date)),
            .groups = "drop")

statin <- med %>%
  left_join(sample_visits %>% distinct(RID, first_plasma), by = "RID") %>%
  group_by(RID) %>%
  summarise(statin_ever = as.integer(any(is_statin)),
            statin_by_first_plasma = as.integer(any(is_statin &
              (is.na(med_date) | med_date <= first_plasma))),
            n_statin_records = sum(is_statin),
            .groups = "drop")
say("Medication file: date column %s, %d records undated; statin records %d for %d people",
    m_mdate, sum(is.na(parse_date_ymd(med_raw[[m_mdate]]))),
    sum(statin$n_statin_records), sum(statin$statin_ever))

# ---- 9. Assemble -------------------------------------------------------------
adni_visit <- sample_visits %>%
  left_join(person, by = "RID") %>%
  left_join(amy,    by = c("RID", "visit_no")) %>%
  left_join(vt,     by = c("RID", "visit_no")) %>%
  left_join(creat,  by = "RID") %>%
  left_join(panel,  by = "RID") %>%
  left_join(statin, by = "RID") %>%
  left_join(statin_draw, by = c("RID", "visit_no")) %>%
  left_join(renal_draw, by = c("RID", "visit_no")) %>%
  left_join(floor_flag, by = c("RID", "visit_no")) %>%
  left_join(wide,   by = c("RID", "visit_no")) %>%
  mutate(age_at_visit = age_bl + as.numeric(visit_date - adni_bl_date) / YEAR_DAYS,
         dx    = factor(dx),
         l_ck  = safe_log(ck),
         l_ast = safe_log(ast),
         l_alt = safe_log(alt),
         # People with medication records but no statin, and draws before a
         # statin was recorded, are non-users.
         statin_at_draw = if_else(!is.na(statin_ever), coalesce(statin_at_draw, 0L),
                                  NA_integer_),
         # People with medication records (statin_ever is not missing) but no
         # kidney-acting drug by a draw are non-users at that draw.
         renal_med_at_draw = if_else(!is.na(statin_ever), coalesce(renal_med_at_draw, 0L),
                                     NA_integer_),
         cohort = "ADNI") %>%
  arrange(RID, visit_no)

# The factors this cohort has, z-scored and signed so that a positive estimate
# always means the direction in the label. Scaling uses the analysis sample, so
# it must happen here and not inside each model.
adni_factors <- FACTORS %>% filter(var %in% names(adni_visit))
for (v in adni_factors$var)
  adni_visit[[paste0("z_", v)]] <- FACTOR_SIGN[[v]] * zf(adni_visit[[v]])

# The same thing in the long shape that the A4 table also uses, so that
# research question 1 can run the same code over both cohorts.
measure_cols <- intersect(MEASURES$measure, names(adni_visit))
adni_long <- adni_visit %>%
  select(cohort, RID, visit_no, plasma_time, age_at_visit, female, e4, dx, centiloid,
         all_of(measure_cols), starts_with("z_"), c2n_at_floor) %>%
  pivot_longer(all_of(measure_cols), names_to = "measure", values_to = "log_value") %>%
  filter(!is.na(log_value)) %>%
  left_join(MEASURES %>% select(measure, analyte, platform), by = "measure") %>%
  mutate(at_floor = measure == "C2N_plasma_ptau217" & coalesce(c2n_at_floor, FALSE)) %>%
  select(-c2n_at_floor) %>%
  arrange(RID, visit_no, analyte, platform)

save_p1(adni_visit, "adni_visit")
save_p1(adni_long,  "adni_long")

# ---- 10. What was built ------------------------------------------------------
n_people <- function(x) n_distinct(adni_visit$RID[!is.na(x)])
chk <- tibble::tribble(
  ~item, ~value,
  "People", n_distinct(adni_visit$RID),
  "Blood draws", nrow(adni_visit),
  "  draws measured by all six platforms", sum(adni_visit$n_platforms == 6, na.rm = TRUE),
  "People with eGFR", n_people(adni_visit$egfr),
  "People with creatine kinase", n_people(adni_visit$ck),
  "People with haemoglobin", n_people(adni_visit$hb),
  "People with BMI at some draw", n_people(adni_visit$bmi),
  "Draws with concurrent amyloid PET", sum(!is.na(adni_visit$centiloid)),
  "Draws with C2N p-tau217 at the detection limit", sum(coalesce(adni_visit$c2n_at_floor, FALSE)),
  "People taking a statin by the first draw", sum(adni_visit$statin_by_first_plasma[!duplicated(adni_visit$RID)], na.rm = TRUE),
  "Draws with a statin recorded by that draw", sum(adni_visit$statin_at_draw, na.rm = TRUE),
  "Statin records (nystatin excluded)", sum(statin$n_statin_records),
  "Distinct medication strings counted as a statin", nrow(statin_strings),
  "Draws with a kidney-acting medication recorded by that draw", sum(adni_visit$renal_med_at_draw, na.rm = TRUE),
  "Distinct medication strings counted as kidney-acting", nrow(renal_strings),
  "People with a screening value taken from an nv re-draw", n_distinct(panel_redraw$RID)
)
print(as.data.frame(chk), row.names = FALSE)
write_csv(chk, file.path(CHECK_DIR, "adni_build.csv"))
write_csv(source_table(), file.path(CHECK_DIR, "adni_sources.csv"))

say("Screening laboratory values, years before the first blood draw: median %.2f, 75th %.2f, 90th %.2f, largest %.2f",
    -median(adni_visit$lab_gap_years, na.rm = TRUE),
    -quantile(adni_visit$lab_gap_years, 0.25, na.rm = TRUE),
    -quantile(adni_visit$lab_gap_years, 0.10, na.rm = TRUE),
    -min(adni_visit$lab_gap_years, na.rm = TRUE))
say("adni_visit: %d rows, %d columns. adni_long: %d rows, %d measures",
    nrow(adni_visit), ncol(adni_visit), nrow(adni_long), n_distinct(adni_long$measure))
})
