# =============================================================================
# 10_figures.R
# Stage 10. Figures 1 and 2 and eFigure 2.
# Writes paper1_out/figures/Figure1, Figure2 and eFigure1 (eFigure 2 of the
# supplement), each as .pdf and .png.
# =============================================================================

# #############################################################################
# Section 10.1. Figures
# Developed as 11_figures.R version 1.9.0, whole
# #############################################################################
local({
# =============================================================================
# 11_figures.R
# The manuscript figures and one supplementary figure.
#
#   Figure 1   Between-platform differences on the same aliquot, by
#              physiological factor, with the equivalence margin drawn on:
#              the six p-tau217 pairs and, below them, Roche against Quanterix
#              for p-tau181, NfL and GFAP.
#
#   Figure 2   The eGFR association modelled flexibly, as the percentage
#              difference from an eGFR of 90 across the observed range: one
#              panel per biomarker (p-tau217, p-tau181, NfL, GFAP), with every
#              platform and cohort overlaid, each with its 95% confidence band.
#
#   eFigure 1  Each of those 13 curves in a panel of its own, with its 95%
#              confidence interval, five panels to a row on a landscape page.
#
# Figure 1 is drawn from the estimates already computed in step 04. The curves
# cannot be, because the tables hold three summary points and no interval
# across the range, so the spline model of step 03 is refitted here on the
# same rows, and the curve and its pointwise interval are computed from the
# fitted coefficients and the same cluster-robust variance used everywhere
# else. Every curve drawn has its row in eTable 8.
#
# All figures are written as PDF, which is vector and is what a journal
# should receive, and as PNG at 300 dpi for looking at. Base graphics are used
# throughout, so no plotting package has to be installed.
#
# Reads    paper1_out/tables/T2_platform_difference.csv
#          paper1_out/adni_visit.rds, a4_screen.rds
# Writes   paper1_out/figures/Figure1, Figure2 and eFigure1, each .pdf and .png
#
# Takes under a minute.
# =============================================================================
SCRIPT_VERSION <- "1.9.0"
cat("11_figures.R version", SCRIPT_VERSION, "\n")

if (!exists("PAPER1_VERSION")) source("00_setup.R")
library(splines)

FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Draw the same plot to both devices, so the two can never drift apart.
# expr is a quoted expression, evaluated once per device. A promise would be
# evaluated only on the first pass and the second file would come out empty.
# Neurology prints a figure at most 6.8 inches wide, so both are drawn at that
# width with 9 point text: what is drawn is what is printed, with no scaling.
FIG_WIDTH <- 6.8
POINTSIZE <- 9
draw_both <- function(name, width, height, expr) {
  for (dev in c("pdf", "png")) {
    f <- file.path(FIG_DIR, paste0(name, ".", dev))
    if (dev == "pdf") grDevices::pdf(f, width = width, height = height, useDingbats = FALSE,
                                     pointsize = POINTSIZE)
    else grDevices::png(f, width = width, height = height, units = "in", res = 300,
                        pointsize = POINTSIZE)
    eval(expr, envir = parent.frame())
    grDevices::dev.off()
    say("  %-14s written", paste0(name, ".", dev))
  }
}

# =============================================================================
# Figure 1. Between-platform differences for p-tau217
# =============================================================================
say("Figure 1, differences between p-tau217 platforms")

T2 <- read_csv(file.path(TAB_DIR, "T2_platform_difference.csv"), show_col_types = FALSE)

nice_platform <- c(C2N_ptau217 = "C2N", Fuji_ptau217 = "Fujirebio",
                   AlzPath_ptau217 = "ALZpath", Janssen_ptau217 = "Janssen")
nice_one <- function(z) unname(ifelse(is.na(nice_platform[z]), z, nice_platform[z]))
nice_pair <- function(x) {
  p <- strsplit(x, " minus ", fixed = TRUE)
  vapply(p, function(z) sprintf("%s vs %s", nice_one(z[1]), nice_one(z[2])), character(1))
}

# Two groups of rows. The six p-tau217 pairs come first; below them, the one
# Roche against Quanterix pair available for each of p-tau181, NfL and GFAP.
# The second group shows the same test reaching equivalence where the platforms
# agree, and the creatine kinase difference recurring in a second analyte.
f1 <- T2 %>%
  filter(sample == "main", role == "primary",
         analyte %in% c("pTau217", "pTau181", "NfL", "GFAP")) %>%
  mutate(factor_short = sub(" \\(.*", "", factor),
         pair_label = ifelse(analyte == "pTau217", nice_pair(pair),
                             c(pTau181 = "p-tau181", NfL = "NfL", GFAP = "GFAP")[analyte]))

FACTOR_ORDER <- c("eGFR", "Body mass index", "Creatine kinase", "Haemoglobin")
PAIR_ORDER <- unique(f1$pair_label[f1$analyte == "pTau217"])
OTHER_ORDER <- c("p-tau181", "NfL", "GFAP")
ROW_ORDER <- c(PAIR_ORDER, OTHER_ORDER)
# Vertical positions: the six pairs, a gap with a heading, then the three
# single pairs. Each group has a heading row of its own above it.
Y_HEAD_1 <- 11
Y_ROWS   <- c(10:5, 3:1)
Y_HEAD_2 <- 4
names(Y_ROWS) <- ROW_ORDER
f1$factor_short <- factor(f1$factor_short, levels = FACTOR_ORDER)
stopifnot(all(ROW_ORDER %in% f1$pair_label), nrow(f1) == 4 * length(ROW_ORDER))

# Verdicts are carried twice, by colour and by fill, so that neither is needed
# alone and the figure reads the same in greyscale or to a colour-blind reader.
# Equivalent and differ are opposite findings and take opposed hues; no
# difference detected is the absence of a finding and stays neutral grey.
#   differ      orange, solid point
#   equivalent  blue,   open point
#   no finding  grey,   solid point
# Checked with the palette validator: blue against orange separates by 24.70
# under colour blindness and 33.60 in normal vision, grey against orange by 9.80
# and 17.60, all above their floors, and each clears 3:1 against white.
COL_DIFFER <- "#eb6834"
COL_EQUIV  <- "#2a78d6"
COL_NONE   <- "#898781"
COL_BAND   <- "#f0efec"   # the margin: a neutral reference region, not data
COL_AXIS   <- "#52514e"
COL_FRAME  <- "#c3c2b7"
point_for <- function(v) {
  ifelse(grepl("^platforms differ", v), 16,
         ifelse(grepl("agree to within", v), 21, 16))
}
colour_for <- function(v) {
  ifelse(grepl("^platforms differ", v), COL_DIFFER,
         ifelse(grepl("agree to within", v), COL_EQUIV, COL_NONE))
}

fig1 <- quote({
  # Two by two. The rows are the same nine comparisons in the same order in
  # every panel, so their names are written once, against the left column.
  op <- par(mfrow = c(2, 2), mar = c(2.6, 0.5, 2.0, 0.6), oma = c(4.6, 9.6, 0.4, 0.4),
            mgp = c(2.0, 0.55, 0), tcl = -0.35, cex = 1)
  # One scale for all four panels, so the factors can be compared by eye.
  xlim <- range(c(f1$pct_lo, f1$pct_hi, -EQUIV_MARGIN_PCT, EQUIV_MARGIN_PCT), na.rm = TRUE)
  xlim <- xlim + c(-1, 1) * 0.04 * diff(xlim)
  for (k in seq_along(FACTOR_ORDER)) {
    fk <- FACTOR_ORDER[k]
    d <- f1[f1$factor_short == fk, , drop = FALSE]
    d <- d[match(ROW_ORDER, d$pair_label), , drop = FALSE]
    y <- unname(Y_ROWS[d$pair_label])
    plot(NA, xlim = xlim, ylim = c(0.4, Y_HEAD_1 + 0.6), axes = FALSE, xlab = "", ylab = "")
    # the equivalence margin
    rect(-EQUIV_MARGIN_PCT, par("usr")[3], EQUIV_MARGIN_PCT, par("usr")[4],
         col = COL_BAND, border = NA)
    abline(v = 0, col = COL_AXIS, lwd = 0.8)
    # a light rule between the two groups of rows
    abline(h = Y_HEAD_2 + 0.5, col = COL_FRAME, lwd = 0.8)
    # Two intervals. The thin line is the 95% interval, which decides whether
    # the platforms differ. The thick line is the 90% interval, which decides
    # equivalence, so a reader can see an equivalent pair's thick line sitting
    # wholly inside the band even where the thin line reaches beyond it.
    segments(d$pct_lo, y, d$pct_hi, y, lwd = 1.1, col = colour_for(d$verdict))
    segments(d$lo90, y, d$hi90, y, lwd = 3.6, col = colour_for(d$verdict), lend = 1)
    points(d$pct_diff, y, pch = point_for(d$verdict), cex = 1.15,
           col = colour_for(d$verdict), bg = "white", lwd = 1.4)
    axis(1, at = pretty(xlim, n = 5), col = COL_FRAME, col.axis = COL_AXIS, cex.axis = 0.95)
    # Row names and group headings sit to the left of the frame, right-aligned
    # against it, and may run into the outer margin so that none is clipped.
    if (k %% 2 == 1) {
      text(par("usr")[1], Y_ROWS, labels = ROW_ORDER, pos = 2,
           offset = 0.45, xpd = NA, cex = 1.0)
      text(par("usr")[1], c(Y_HEAD_1, Y_HEAD_2), labels = c("p-tau217", "Roche vs Quanterix"),
           pos = 2, offset = 0.45, xpd = NA, cex = 1.0, font = 2)
    }
    title(main = fk, line = 0.6, cex.main = 1.08, font.main = 2)
    box(col = COL_FRAME)
  }
  mtext("Difference between platforms, % per standard deviation of the factor",
        side = 1, outer = TRUE, line = 1.3, cex = 1.0)
  # A legend across the foot of the whole figure rather than inside one panel.
  # Centred under the panels, like the axis title, not under the whole page.
  mid_ndc <- grconvertX(0.5, from = "nic", to = "ndc")
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  leg_txt <- c("platforms differ", "equivalent within the margin",
               "no difference detected")
  # Drawn item by item so the gaps between the three items are equal.
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  sym_w <- strwidth("M", cex = 0.95) * 1.2   # symbol plus the space after it
  gap   <- strwidth("M", cex = 0.95) * 2.5   # space between items
  item_w <- sym_w + strwidth(leg_txt, cex = 0.95)
  x0 <- grconvertX(mid_ndc, from = "ndc", to = "user") -
    (sum(item_w) + gap * (length(leg_txt) - 1)) / 2
  y0 <- grconvertY(0.018, from = "ndc", to = "user")
  for (m in seq_along(leg_txt)) {
    points(x0 + sym_w * 0.3, y0, pch = c(16, 21, 16)[m], cex = 1.0,
           col = c(COL_DIFFER, COL_EQUIV, COL_NONE)[m], bg = "white",
           lwd = c(1, 1.4, 1)[m], xpd = NA)
    text(x0 + sym_w, y0, leg_txt[m], adj = c(0, 0.5), cex = 0.95,
         col = COL_AXIS, xpd = NA)
    x0 <- x0 + item_w[m] + gap
  }
  # The meaning of the thick and thin lines is given in the figure legend of the
  # manuscript, where there is room for it.
  par(op)
})
draw_both("Figure1", FIG_WIDTH, 7.6, fig1)

# =============================================================================
# Figure 2. The shape of the eGFR association
# =============================================================================
say("Figure 2, the eGFR association modelled flexibly")

adni <- readRDS(file.path(OUT_DIR, "adni_visit.rds"))
a4   <- readRDS(file.path(OUT_DIR, "a4_screen.rds"))

KNOTS <- c(45, 60, 75)
REFERENCE <- 90          # the eGFR the curve is drawn as a difference from

# The curve and its pointwise interval. The spline basis is built once and put
# into the data, so that the basis at any eGFR on the grid is the same object
# the model was fitted with. The quantity plotted is a contrast between two
# points on the curve, so every covariate cancels and only the spline columns
# contribute.
# C2N p-tau217 values at the detection limit are not measurements. From
# version 1.9.0 its curve, like its row in eTable 8, is fitted by left-censored
# regression on every sample (spline_censored in 00_setup.R) rather than by
# least squares on the samples above the limit, which selected on the outcome.
# drop_c2n_floor is kept, but now returns the data unchanged, so the grid below
# is set from every sample the censored curve uses.
drop_c2n_floor <- function(d, measure) d
spline_curve <- function(d, measure, covs, id, grid) {
  if (measure == "C2N_plasma_ptau217") {
    r <- spline_censored(d, measure, covs, id, at = grid, knots = KNOTS, reference = REFERENCE)
    return(tibble(egfr = grid, pct = r$points$pct, lo = r$points$lo, hi = r$points$hi,
                  n_obs = r$n_obs, n_people = r$n_people))
  }
  d <- drop_c2n_floor(as.data.frame(d), measure)
  # The same rows step 03 uses for eTable 8, so the curve and the table are the
  # same model: the measure, eGFR and amyloid PET must all be present before the
  # spline basis is built, because the basis boundaries are taken from the eGFR
  # values of the rows it is built on.
  d <- d[!is.na(d$egfr) & !is.na(d[[measure]]) & !is.na(d$centiloid), , drop = FALSE]
  # Row names are reset so that positions and names agree. Without this the
  # rows kept after filtering carry their original names, and indexing by those
  # names selects the wrong participants.
  rownames(d) <- NULL
  if (nrow(d) < 100) return(NULL)
  B <- ns(d$egfr, knots = KNOTS, Boundary.knots = range(d$egfr))
  sp <- paste0("sp", seq_len(ncol(B)))
  for (i in seq_along(sp)) d[[sp[i]]] <- B[, i]
  f <- stats::as.formula(sprintf("`%s` ~ %s + %s", measure,
                                 paste(sp, collapse = " + "), covs))
  mf  <- stats::model.frame(f, data = d, na.action = stats::na.omit)
  idx <- as.integer(rownames(mf))
  dd  <- d[idx, , drop = FALSE]
  m   <- stats::lm(f, data = dd)
  cf  <- stats::coef(m); keep <- !is.na(cf)
  X <- stats::model.matrix(m)[, keep, drop = FALSE]
  e <- stats::residuals(m); cl <- dd[[id]]
  G <- length(unique(cl)); n <- nrow(X); k <- ncol(X)
  bread <- solve(crossprod(X))
  meat  <- Reduce(`+`, lapply(split(seq_len(n), cl), function(i) {
    s <- crossprod(X[i, , drop = FALSE], e[i]); s %*% t(s) }))
  V <- bread %*% meat %*% bread * (G / (G - 1)) * ((n - 1) / (n - k))
  b <- cf[keep]
  nm <- names(b)
  pos <- match(sp, nm)
  if (anyNA(pos)) return(NULL)
  b0 <- as.numeric(predict(B, newx = REFERENCE))
  out <- lapply(grid, function(g) {
    cvec <- rep(0, length(b))
    cvec[pos] <- as.numeric(predict(B, newx = g)) - b0
    est <- sum(cvec * b)
    se  <- sqrt(as.numeric(t(cvec) %*% V %*% cvec))
    c(est = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
  })
  out <- do.call(rbind, out)
  tibble(egfr = grid, pct = pct(out[, "est"]),
         lo = pct(out[, "lo"]), hi = pct(out[, "hi"]),
         n_obs = nrow(dd), n_people = G)
}

# Every curve is its own model, one biomarker on one platform in one cohort,
# drawn in a panel of its own. Curves are not overlaid: their intervals come
# from separate models and overlap heavily, so an overlaid figure invites a
# comparison between platforms that only the same-aliquot test in Figure 1 can
# make.
#   Figure 2   p-tau217 on the four ADNI platforms, two by two.
#   eFigure 1  p-tau181, NfL and GFAP, one row each, with one column for each
#              of Roche in A4 and LEARN, Roche in ADNI and Quanterix in ADNI.
# Within a figure the x axis is shared, and the y axis is shared along the
# comparison a reader makes: across all four panels of Figure 2, and along
# each row of eFigure 1.
SPECS <- tibble::tribble(
  ~key,             ~biomarker, ~column,                  ~data,  ~measure,
  "C2N",            "p-tau217", "C2N",                    "adni", "C2N_plasma_ptau217",
  "Fujirebio",      "p-tau217", "Fujirebio",              "adni", "Fuji_plasma_ptau217",
  "ALZpath",        "p-tau217", "ALZpath",                "adni", "AlzPath_plasma_ptau217",
  "Janssen",        "p-tau217", "Janssen",                "adni", "Janssen_plasma_ptau217",
  "ptau181_a4",     "p-tau181", "A4 and LEARN, Roche",    "a4",   "pTau181",
  "ptau181_roche",  "p-tau181", "ADNI, Roche",            "adni", "Roche_plasma_ptau181",
  "ptau181_qx",     "p-tau181", "ADNI, Quanterix",        "adni", "QX_plasma_ptau181",
  "nfl_a4",         "NfL",      "A4 and LEARN, Roche",    "a4",   "NfL",
  "nfl_roche",      "NfL",      "ADNI, Roche",            "adni", "Roche_plasma_NfL",
  "nfl_qx",         "NfL",      "ADNI, Quanterix",        "adni", "QX_plasma_NfL",
  "gfap_a4",        "GFAP",     "A4 and LEARN, Roche",    "a4",   "GFAP",
  "gfap_roche",     "GFAP",     "ADNI, Roche",            "adni", "Roche_plasma_GFAP",
  "gfap_qx",        "GFAP",     "ADNI, Quanterix",        "adni", "QX_plasma_GFAP")

curves <- lapply(seq_len(nrow(SPECS)), function(k) {
  p <- SPECS[k, ]
  d <- if (p$data == "a4") a4 else adni
  if (!p$measure %in% names(d)) return(NULL)
  covs <- if (p$data == "a4") COVS_A4 else COVS_ADNI
  id   <- if (p$data == "a4") "BID" else "RID"
  # The grid stops at the 1st and 99th centiles of eGFR in the rows the curve
  # is fitted on, because a spline says little beyond the data.
  d <- drop_c2n_floor(as.data.frame(d), p$measure)
  fit_rows <- d[!is.na(d$egfr) & !is.na(d[[p$measure]]) & !is.na(d$centiloid), , drop = FALSE]
  lim <- stats::quantile(fit_rows$egfr, c(0.01, 0.99), na.rm = TRUE)
  g <- seq(max(25, lim[1]), min(120, lim[2]), length.out = 120)
  r <- spline_curve(d, p$measure, covs, id, g)
  if (is.null(r)) return(NULL)
  r$key <- p$key
  r
})
names(curves) <- SPECS$key
missing_curves <- SPECS$key[vapply(curves, is.null, logical(1))]
if (length(missing_curves) > 0)
  stop("These curves could not be fitted: ", paste(missing_curves, collapse = ", "))
for (nm in names(curves))
  say("  %-14s %d observations, %d people, eGFR %.2f to %.2f", nm,
      curves[[nm]]$n_obs[1], curves[[nm]]$n_people[1],
      min(curves[[nm]]$egfr), max(curves[[nm]]$egfr))

# Each platform keeps one colour wherever it appears, so a reader can follow
# Roche and Quanterix across eFigure 1 and tell the four p-tau217 platforms
# apart in Figure 2. The curve is drawn in the platform colour and its interval
# in a pale tint of the same colour. The six hues are the reference
# categorical palette, checked with the palette validator: Roche against
# Quanterix, and the four p-tau217 platforms against one another, all clear the
# colour-blind and normal-vision floors. ALZpath, Janssen and Quanterix fall
# below 3:1 against white, so no panel relies on colour alone: every panel is
# named by its platform in its title or column heading.
COL_PLATFORM <- c(Roche = "#2a78d6", Quanterix = "#1baf7a", C2N = "#4a3aa7",
                  Fujirebio = "#008300", ALZpath = "#e87ba4", Janssen = "#eda100")
platform_of <- function(key) {
  col <- SPECS$column[match(key, SPECS$key)]
  ifelse(key %in% names(COL_PLATFORM), key, sub(".*, ", "", col))
}
XLIM <- range(unlist(lapply(curves, `[[`, "egfr")))

y_range <- function(keys) {
  yl <- range(unlist(lapply(curves[keys], function(r) c(r$lo, r$hi, 0))), na.rm = TRUE)
  yl + c(-1, 1) * 0.04 * diff(yl)
}
one_panel <- function(r, ylim, col, main = NULL, xaxis = TRUE, yaxis = TRUE, xat = NULL) {
  plot(NA, xlim = XLIM, ylim = ylim, axes = FALSE, xlab = "", ylab = "")
  polygon(c(r$egfr, rev(r$egfr)), c(r$lo, rev(r$hi)), col = adjustcolor(col, alpha.f = 0.20), border = NA)
  abline(h = 0, col = COL_AXIS, lwd = 0.8)
  abline(v = REFERENCE, col = COL_AXIS, lty = 3, lwd = 0.8)
  lines(r$egfr, r$pct, lwd = 2, col = col)
  axis(1, at = if (is.null(xat)) axTicks(1) else xat, col = COL_FRAME, col.axis = COL_AXIS,
       cex.axis = 0.95, labels = xaxis)
  axis(2, las = 1, col = COL_FRAME, col.axis = COL_AXIS, cex.axis = 0.95, labels = yaxis)
  if (!is.null(main)) title(main = main, line = 0.6, cex.main = 1.08, font.main = 2)
  box(col = COL_FRAME)
}

# ---- Figure 2 ----------------------------------------------------------------
# One panel per biomarker, every platform and cohort overlaid, each curve with
# its 95% confidence interval as a translucent band. Colour marks the platform;
# ADNI curves are solid and A4 and LEARN curves dashed. The bands come from
# separate models, so their overlap is not a test of agreement between
# platforms; that test is the same-aliquot comparison in Figure 1, and the
# figure legend says so.
F2_PANELS <- c("p-tau217", "p-tau181", "NfL", "GFAP")
cohort_of <- function(key) SPECS$data[match(key, SPECS$key)]
fig2 <- quote({
  op <- par(mfrow = c(2, 2), mar = c(2.6, 3.0, 2.0, 0.6), oma = c(4.2, 1.8, 0.4, 0.4),
            mgp = c(2.0, 0.55, 0), tcl = -0.35, cex = 1)
  for (b in F2_PANELS) {
    keys <- SPECS$key[SPECS$biomarker == b]
    yl <- range(unlist(lapply(curves[keys], function(r) c(r$lo, r$hi))))
    plot(NA, xlim = c(38, 104), ylim = yl, axes = FALSE, xlab = "", ylab = "")
    abline(h = 0, col = COL_AXIS, lwd = 0.8)
    abline(v = REFERENCE, col = COL_AXIS, lty = 3, lwd = 0.8)
    # bands first, then every curve on top of every band
    for (k in keys) {
      r <- curves[[k]]
      polygon(c(r$egfr, rev(r$egfr)), c(r$lo, rev(r$hi)),
              col = adjustcolor(COL_PLATFORM[[platform_of(k)]], alpha.f = 0.10), border = NA)
    }
    for (k in keys) lines(curves[[k]]$egfr, curves[[k]]$pct, lwd = 2,
                          col = COL_PLATFORM[[platform_of(k)]],
                          lty = if (cohort_of(k) == "a4") 2 else 1)
    axis(1, col = COL_FRAME, col.axis = COL_AXIS)
    axis(2, las = 1, col = COL_FRAME, col.axis = COL_AXIS)
    box(col = COL_FRAME)
    title(main = b, line = 0.6, cex.main = 1.08)
  }
  mtext(expression("Estimated glomerular filtration rate, mL/min/1.73 m"^2),
        side = 1, outer = TRUE, line = 0.8)
  mtext("Difference from an eGFR of 90, %", side = 2, outer = TRUE, line = 0.4)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  leg_txt <- c(names(COL_PLATFORM), "ADNI (solid)", "A4 and LEARN (dashed)")
  # Column widths set from the longest label, so the legend stays well inside
  # the page whatever the device.
  legend("bottom", horiz = FALSE, ncol = 4, bty = "n", cex = 0.9, inset = 0.005,
         legend = leg_txt, text.width = max(strwidth(leg_txt, cex = 0.9)) * 0.62,
         col = c(COL_PLATFORM, COL_AXIS, COL_AXIS), lwd = 2,
         lty = c(rep(1, length(COL_PLATFORM)), 1, 2), seg.len = 2.2)
  par(op)
})
draw_both("Figure2", FIG_WIDTH, 5.8, fig2)

# ---- eFigure 1 ---------------------------------------------------------------
# Every curve with its interval, one panel each: the 13 panels in reading
# order, five to a row, on a landscape page, since the supplement is not bound
# by the journal's figure limits. Each panel is headed by its biomarker and,
# beneath, its cohort and platform. The y axis is shared by all panels of one
# biomarker, so its curves can be compared across cohorts and platforms.
E1_ORDER <- c("C2N", "Fujirebio", "ALZpath", "Janssen",
              "ptau181_a4", "ptau181_roche", "ptau181_qx",
              "nfl_a4", "nfl_roche", "nfl_qx",
              "gfap_a4", "gfap_roche", "gfap_qx")
E1_PER_ROW <- 5
E1_HEAD <- function(k) {
  b <- SPECS$biomarker[match(k, SPECS$key)]
  sub_ <- if (k %in% c("C2N", "Fujirebio", "ALZpath", "Janssen")) paste0("ADNI, ", k)
          else SPECS$column[match(k, SPECS$key)]
  c(b, sub_)
}
E1_YLIM <- lapply(split(SPECS$key, SPECS$biomarker), y_range)
efig1 <- quote({
  n_rows <- ceiling(length(E1_ORDER) / E1_PER_ROW)
  op <- par(mfrow = c(n_rows, E1_PER_ROW), mar = c(2.2, 2.6, 2.6, 0.5),
            oma = c(2.6, 1.8, 0.2, 0.4), mgp = c(2.0, 0.55, 0), tcl = -0.35, cex = 1)
  for (k in E1_ORDER) {
    b <- SPECS$biomarker[match(k, SPECS$key)]
    # Ticks every 20 in these narrow panels, so every tick keeps its label.
    one_panel(curves[[k]], E1_YLIM[[b]], col = COL_PLATFORM[[platform_of(k)]],
              xat = seq(40, 100, 20))
    h <- E1_HEAD(k)
    mtext(h[1], side = 3, line = 1.3, cex = 0.9, font = 2)
    mtext(h[2], side = 3, line = 0.35, cex = 0.8)
  }
  mtext(expression("Estimated glomerular filtration rate, mL/min/1.73 m"^2),
        side = 1, outer = TRUE, line = 1.2, cex = 1.0)
  mtext("Difference from an eGFR of 90, %", side = 2, outer = TRUE, line = 0.4, cex = 1.0)
  par(op)
})
draw_both("eFigure1", 10.0, 6.6, efig1)

say("Written to %s", FIG_DIR)
say("Figure 1 is drawn from the estimates in T2_platform_difference.csv, so it")
say("cannot disagree with Table 3. Figure 2 and eFigure 1 refit the spline")
say("of step 03 on the same rows, so each curve passes through the values in eTable 8.")
})
