# =============================================================================
# run_all.R
# The driver. Runs the Paper 1 pipeline in order, from the raw releases to the
# tables and figures, and writes a log of the run. Source it by its full path,
# from any working directory, or open it in RStudio and press Source.
#
# Before sourcing, set where the data are and where the output should go:
#
#   Sys.setenv(PAPER1_ADNI_DIR = "<folder with the ADNI files>",
#              PAPER1_A4_DIR   = "<folder with the A4/LEARN release>",
#              PAPER1_OUT_DIR  = "<output folder>")
#
# Without them, 00_setup.R looks for ADNI/ and A4&learn/ beside the folder
# that holds this code and writes to paper1_out/ there.
#
# STEPS   which stages to run, by their two-digit prefix; "all" runs every
#         stage. Each stage reads what the earlier ones wrote, so a later
#         stage can be run on its own once the earlier ones have run.
#           01 ADNI data            06 change on change
#           02 A4/LEARN data        07 sensitivity and robustness
#           03 levels               08 intervals added in the review
#           04 platform differences 09 formatted tables
#           05 ratios               10 figures
# QUICK   TRUE turns off the bootstraps of the joint growth model and skips
#         stage 08, for checking the point estimates in about an hour (the
#         1,000-resample bootstraps of sections 4.2, 5.4 and 7.3, which take
#         minutes, still run). FALSE, the default, runs the paper's bootstraps:
#         1,000 resamples of participants each, seed 1. The full run takes
#         about a day, almost all of it in the joint growth model bootstraps of
#         stages 06 to 08; the Mac is kept awake while R runs.
#
# Write each run to an empty output folder: a rerun into a used one reuses
# some earlier results (the joint model bootstrap rows of section 6.1, the
# comorbidity interval of section 7.2 when its bootstrap is off, the saved
# tasks of section 8.2 and the laboratory extract of section 2.1).
#
# Every bootstrap setting is set here, for the paper's run or for QUICK, so an
# earlier QUICK run in the same R session cannot carry its zeros into a full
# run. Each stage writes its own tables; the log of the run, with every
# message, warning and error, and the session information (R and package
# versions) are written to <output folder>/logs/.
# =============================================================================
if (!exists("STEPS")) STEPS <- "all"
if (!exists("QUICK")) QUICK <- FALSE

# ---- find this file's folder and move there ---------------------------------
here <- NULL
for (i in seq_len(sys.nframe())) {
  f <- sys.frame(i)$ofile
  if (!is.null(f)) here <- f
}
if (is.null(here)) {
  a <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(a)) here <- sub("^--file=", "", a[1])
}
if (is.null(here) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  here <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL)
if (is.null(here) || !nzchar(here) || !file.exists(here)) stop(
  "\n  This file cannot work out where it is. Either press Source in RStudio, or run",
  "\n  source() with the full path of run_all.R.")
setwd(dirname(normalizePath(here)))

# ---- the stages -------------------------------------------------------------
ALL_STEPS <- c(
  "01" = "01_data_adni.R",
  "02" = "02_data_a4learn.R",
  "03" = "03_levels.R",
  "04" = "04_platform_differences.R",
  "05" = "05_ratios.R",
  "06" = "06_change.R",
  "07" = "07_sensitivity_and_robustness.R",
  "08" = "08_intervals.R",
  "09" = "09_tables.R",
  "10" = "10_figures.R")
want <- if (identical(STEPS, "all")) names(ALL_STEPS) else as.character(STEPS)
want <- intersect(names(ALL_STEPS), want)
if (isTRUE(QUICK)) want <- setdiff(want, "08")
if (length(want) == 0) stop("STEPS does not name any stage. Use \"all\" or some of: ",
                            paste(names(ALL_STEPS), collapse = ", "))
missing_files <- ALL_STEPS[want][!file.exists(ALL_STEPS[want])]
if (length(missing_files)) stop("Not in this folder: ", paste(missing_files, collapse = ", "))

# ---- bootstrap settings -------------------------------------------------------
# The paper's settings: 1,000 resamples for every interval, the joint growth
# model bootstrapped for eGFR, creatine kinase and haemoglobin (Table 5,
# eTable 19), seed 1.
n_boot <- if (isTRUE(QUICK)) "0" else "1000"
Sys.setenv(PAPER1_BOOT = n_boot, PAPER1_BOOT_VARS = "egfr,l_ck,hb", PAPER1_BOOT_SEED = "1",
           PAPER1_N_BOOT = n_boot, PAPER1_COMORB_BOOT = n_boot,
           PAPER1_REVIEW_BOOT = n_boot, PAPER1_REVIEW_ARM_BOOT = n_boot,
           PAPER1_ARM_R_BOOT = n_boot, PAPER1_MI_BOOT = "1000")

source("00_setup.R")

# Keep the Mac awake while this R session runs (macOS only; nothing is changed
# when R exits).
if (Sys.info()[["sysname"]] == "Darwin")
  try(system2("caffeinate", c("-i", "-w", Sys.getpid()), wait = FALSE), silent = TRUE)

# ---- log file ----------------------------------------------------------------
# Everything printed, and every message, warning and error, is also written to
# a log file named after the stages and the start time. Messages and warnings
# still appear in the console; an error is copied to the log and then stops
# the run. The log is closed whatever happens.
LOG_DIR <- file.path(OUT_DIR, "logs")
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
LOG_FILE <- file.path(LOG_DIR, sprintf("run_all_%s%s_%s.log", paste(want, collapse = "-"),
                                       if (isTRUE(QUICK)) "_quick" else "", stamp))
log_con    <- file(LOG_FILE, open = "wt")
sink_level <- sink.number()
sink(log_con, split = TRUE)
log_line <- function(...) cat(..., file = log_con, sep = "")
tryCatch(
  withCallingHandlers({
    cat("Log of run_all.R, written to ", LOG_FILE, "\n", sep = "")
    cat("Paper 1, version ", PAPER1_VERSION, "; ", R.version.string, "\n", sep = "")
    cat("  code    ", CODE_DIR, "\n  ADNI    ", paste(ADNI_DIRS, collapse = "; "),
        "\n  A4      ", A4_DIR, "\n  output  ", OUT_DIR, "\n", sep = "")
    cat("Bootstraps: ", if (isTRUE(QUICK)) "off (QUICK = TRUE)" else "the paper's settings", "\n", sep = "")
    env <- Sys.getenv()
    env <- env[grepl("^PAPER1_", names(env))]
    if (length(env)) cat("Settings: ", paste(names(env), env, sep = " = ", collapse = "; "), "\n", sep = "")

    needed <- c(survival = "the censored regressions (stages 03 to 05, 07)",
                lme4     = "the reliability models (stages 03 and 06)",
                nlme     = "the joint growth models (stages 06 to 08)",
                lavaan   = "the MIMIC model (stages 04 and 05; skipped with a message if absent)",
                splines  = "the eGFR splines (stages 03 and 10)",
                parallel = "the interval tasks of stage 08")
    miss <- names(needed)[!vapply(names(needed), requireNamespace, logical(1), quietly = TRUE)]
    if (length(miss) > 0) {
      cat("\nThese packages are not installed. The stages that need them will stop:\n")
      for (m in miss) cat("  ", m, " for ", needed[[m]], "\n", sep = "")
      cat("  install.packages(c(", paste(sprintf("\"%s\"", miss), collapse = ", "), "))\n\n", sep = "")
    }

    cat("\nStages: ", paste(want, collapse = ", "), "\nStarted at ", format(Sys.time()), "\n", sep = "")
    for (s in want) {
      f <- ALL_STEPS[[s]]
      cat("\n############################################################\n## ", f,
          "   ", format(Sys.time()), "\n############################################################\n", sep = "")
      flush(log_con)
      source(f, echo = FALSE)
    }
    cat("\nFinished at ", format(Sys.time()), "\nOutput is in ", OUT_DIR, "\n", sep = "")
    si <- utils::capture.output(print(utils::sessionInfo()))
    writeLines(si, file.path(LOG_DIR, sprintf("session_info_%s.txt", stamp)))
    cat("\n", paste(si, collapse = "\n"), "\n", sep = "")
  },
  message   = function(m) { x <- conditionMessage(m); log_line(x, if (!grepl("\n$", x)) "\n") },
  warning   = function(w) log_line("Warning: ", conditionMessage(w), "\n"),
  error     = function(e) log_line("\nError at ", format(Sys.time()), ": ", conditionMessage(e), "\n"),
  interrupt = function(i) log_line("\nStopped by the user at ", format(Sys.time()), "\n")),
  finally = {
    while (sink.number() > sink_level) sink()
    close(log_con)
    cat("The log of this run is in ", LOG_FILE, "\n", sep = "")
  })
