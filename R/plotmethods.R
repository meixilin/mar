# a collection of plot methods
.anncol <- "#A9A9A9"
.theorycol <- "#377EB8"
# martheory data frames only carry N, M, thetaw, thetapi (no E, A, or Asq)
.Mtype_theory <- c("M", "thetaw", "thetapi")

# ---- plot methods -----------------------------------------------------------

# define the plotting method for marmaps
# methods(class = "marmaps") > [1] plot
#' Plot Method for marmaps Class
#'
#' Creates a spatial plot of sample locations and density
#'
#' @param x An object of class "marmaps" containing sample mapping information
#' @param ... Additional arguments passed to plot
#'
#' @return Invisibly returns NULL
#' @export
#'
#' @examples
#' plot(gm1001g$maps)
plot.marmaps <- function(x, ...) {
    sm = .get_samplemap(x)
    terra::plot(sm, plg = list(title = "# of samples"), ...)
    points(x$lonlat, pch = 20, col = "#D3D3D380")
    return(invisible())
}

#' Plot Method for Site Frequency Spectrum
#'
#' Creates a line/point plot of the site frequency spectrum
#'
#' @param x An object of class "sfs" containing site frequency spectrum data
#' @param expected What expected (neutral) SFS to overlay: `NULL` (default, no
#'   overlay), an `sfs` object (output of \link{expsfs}, used as-is), or a
#'   `genomaps` object (an `sfs` is generated from it via \link{expsfs}). When
#'   supplied, a legend distinguishing observed vs. expected is added in the
#'   top right corner.
#' @param ... Additional arguments passed to \link[graphics]{plot}, other than
#'   `type`, `xlab` and `ylab`, which the method sets itself. To examine rare
#'   variants more closely, pass `log = "x"` to spread out the low allele-count
#'   bins.
#'
#' @return Invisibly returns NULL
#' @export
#'
#' @examples
#' sfs_object <- sfs(AC = c(1, 1, 0, 2, 0, 1, 1, 0, 0, 2, 30), N = 50, ploidy = 2)
#' plot(sfs_object)
#' \dontrun{
#' # overlay the neutral expectation, generated automatically from gm
#' plot(sfs_object, expected = gm1001g)
#' # spread out rare variants (low allele counts) on a log-x scale
#' plot(sfs_object, log = "x")
#' }
plot.sfs <- function(x, expected = NULL, ...) {
    data <- as.vector(x)
    bins <- as.integer(names(x))
    # plot
    graphics::plot(x = bins, y = data, type = "b", xlab = "Allele Count", ylab = "Number of Alleles", ...)
    expected <- .resolve_expsfs(x, expected)
    if (!is.null(expected)) {
        .overlay_expsfs(bins, x, expected)
    }
    # dataline = 1 because the observed spectrum is drawn as points AND line
    .ann_legend("topright", list(...), "observed",
        overlaylab = if (!is.null(expected)) "expected", dataline = 1
    )
    return(invisible())
}

# this also works for classes belonging to marextinct (also a sampling process)
#' Plot Method for Mutations-Area Relationship sample data
#'
#' Creates a plot of genetic diversity against area, optionally with fitted power law
#'
#' @param x A data frame containing sampling data
#' @param fit Whether to draw the fitted power-law curve. `FALSE` (default) draws
#'   nothing; `TRUE` fits both `c` and `z` via \link{MARcalc}; a named
#'   `c(c = , z = )` vector supplies one or both directly, and leaves the rest
#'   to be fitted via \link{MARcalc}.
#' @param Mtype Type of genetic diversity metric to plot
#' @param Atype Type of area metric to plot
#' @param theory What theoretical MAR prediction to overlay: `NULL` (default,
#'   no overlay), a `martheory` object (output of \link{MARtheory}, used
#'   as-is), or a `genomaps` object (a `martheory` is generated from it via
#'   \link{MARtheory}). Only supported when `Atype = "N"`. When supplied
#'   (together with `fit`), a legend distinguishing observed data / fitted
#'   equation / theory is added in the bottom right corner.
#' @param ... Additional arguments passed to \link[graphics]{plot}. Pass
#'   `log = "xy"` (or `"x"` / `"y"`) for a log-log (or semi-log) plot; the
#'   fitted curve follows the axes.
#'
#' @return Invisibly returns NULL
#' @export
#'
#' @examples
#' \dontrun{
#' plot(marsamp_object, fit = c(c = 0.5, z = 0.25))
#' # fit c/z automatically, and overlay the MARtheory prediction generated from gm
#' plot(marsamp_object, Atype = "N", fit = TRUE, theory = gm1001g)
#' # the power law as a straight line on log-log axes
#' plot(marsamp_object, fit = TRUE, log = "xy")
#' }
plot.marsamp <- function(x, fit = FALSE, Mtype = .Mtype, Atype = .Atype, theory = NULL, ...) {
    Mtype <- match.arg(Mtype)
    Atype <- match.arg(Atype)
    tmpdf <- x[, c(Atype, Mtype)]
    tmpdf <- tmpdf[(tmpdf[, Mtype] > 0 & !is.na(tmpdf[, Mtype])), ]
    # previously had a check for length(unique(tmpdf[,Mtype])) < 4
    if (nrow(tmpdf) == 0) {
        warning(paste0("MAR for Mtype = ", Mtype, ", Atype = ", Atype, " cannot be plotted"))
        graphics::plot.new()
        return(invisible())
    }
    cz <- .resolve_fit(tmpdf, fit, Mtype = Mtype, Atype = Atype)
    c <- cz$c
    z <- cz$z
    # plot
    graphics::plot(x = tmpdf[, Atype], y = tmpdf[, Mtype], xlab = Atype, ylab = Mtype, ...)
    # M = c*A^z. curve() follows the axes it is drawn on, so this stays correct
    # when the caller asks for log = "x" / "y" / "xy"
    if (!is.null(c) & !is.null(z)) {
        curve(c * x^z, add = TRUE, col = .anncol)
    }
    theory <- .resolve_theory(theory)
    if (!is.null(theory)) {
        .overlay_theory_marsamp(theory, Mtype = Mtype, Atype = Atype)
    }
    .ann_legend("bottomright", list(...), "observed",
        eq = if (!is.null(c) & !is.null(z)) .eq_marsamp(c, z, Mtype = Mtype, Atype = Atype),
        overlaylab = if (!is.null(theory)) "theory"
    )
    return(invisible())
}

#' Plot Method for Mutations Extinction Curves
#'
#' Creates a plot showing relationship between area loss and genetic diversity loss
#'
#' @param x A data frame containing extinction data
#' @param fit Whether to draw the fitted extinction curve. `FALSE` (default)
#'   draws nothing; `TRUE` fits `z` via \link{MARcalc}; a named `c(z = )`
#'   vector supplies it directly.
#' @param Mtype Type of genetic diversity metric to plot
#' @param Atype Type of area metric to plot
#' @param theory What theoretical MAR prediction to overlay (rescaled to the
#'   same % lost / % remained axes): `NULL` (default, no overlay), a
#'   `martheory` object (output of \link{MARtheory}, used as-is), or a
#'   `genomaps` object (a `martheory` is generated from it via
#'   \link{MARtheory}). Only supported when `Atype = "N"`.
#' @param ... Additional arguments passed to plot
#'
#' @return Invisibly returns NULL
#' @export
#'
#' @examples
#' \dontrun{
#' plot(marextinct_object, fit = c(z = 0.25))
#' # fit z automatically, and overlay the MARtheory prediction generated from gm
#' plot(marextinct_object, Atype = "N", fit = TRUE, theory = gm1001g)
#' }
plot.marextinct <- function(x, fit = FALSE, Mtype = .Mtype, Atype = .Atype, theory = NULL, ...) {
    Mtype <- match.arg(Mtype)
    Atype <- match.arg(Atype)
    # remove NA or zero data
    tmpdf <- x[, c(Atype, Mtype)]
    tmpdf <- tmpdf[(tmpdf[, Mtype] > 0 & !is.na(tmpdf[, Mtype])), ]
    # previously had a check for length(unique(tmpdf[,Mtype])) < 4
    if (nrow(tmpdf) == 0) {
        stop(paste0("MAR for Mtype = ", Mtype, ", Atype = ", Atype, " cannot be plotted"))
    }
    # generate area loss data (scale by the first value not the max value)
    a_per <- 1 - tmpdf[, Atype] / (tmpdf[1, Atype])
    m_per <- tmpdf[, Mtype] / (tmpdf[1, Mtype])
    z <- .resolve_fit(tmpdf, fit, Mtype = Mtype, Atype = Atype)$z
    # plot
    # m_per = 1 - (1-a_per)^z
    graphics::plot(x = a_per, y = m_per, xlab = paste0("% of ", Atype, " lost"), ylab = paste0("% of ", Mtype, " remained"), ...)
    if (!is.null(z)) {
        curve((1 - x)^z, add = TRUE, col = .anncol)
    }
    theory <- .resolve_theory(theory)
    if (!is.null(theory)) {
        .overlay_theory_marextinct(theory, Mtype = Mtype, Atype = Atype)
    }
    .ann_legend("bottomleft", list(...), "observed",
        eq = if (!is.null(z)) bquote((1 - m) == (1 - a)^.(round(z, 2))),
        overlaylab = if (!is.null(theory)) "theory"
    )
    return(invisible())
}

#' Plot Method for Theoretical Mutations-Area Relationship Predictions
#'
#' Creates a plot of the theoretical MAR prediction (from \link{MARtheory})
#' against sample size, optionally with a fitted power-law curve.
#'
#' @param x An object of class "martheory" (output of \link{MARtheory})
#' @param fit Whether to draw the fitted power-law curve. `FALSE` (default) draws
#'   nothing; `TRUE` fits both `c` and `z` via \link{MARcalc}; a named
#'   `c(c = , z = )` vector supplies one or both directly, and leaves the rest
#'   to be fitted via \link{MARcalc}.
#' @param Mtype Type of genetic diversity metric to plot. Allowed values are `r toString(.Mtype_theory)`.
#' @param ... Additional arguments passed to \link[graphics]{plot}. Pass
#'   `log = "xy"` (or `"x"` / `"y"`) for a log-log (or semi-log) plot; the
#'   fitted curve follows the axes.
#'
#' @return Invisibly returns NULL
#' @export
#'
#' @examples
#' \dontrun{
#' theory <- MARtheory(gm1001g)
#' plot(theory, fit = TRUE)
#' plot(theory, fit = TRUE, log = "xy")
#' }
plot.martheory <- function(x, fit = FALSE, Mtype = .Mtype_theory, ...) {
    Mtype <- match.arg(Mtype)
    Atype <- "N" # martheory is always predicted against N (number of individuals)
    tmpdf <- x[, c(Atype, Mtype)]
    tmpdf <- tmpdf[(tmpdf[, Mtype] > 0 & !is.na(tmpdf[, Mtype])), ]
    if (nrow(tmpdf) == 0) {
        warning(paste0("MARtheory for Mtype = ", Mtype, " cannot be plotted"))
        graphics::plot.new()
        return(invisible())
    }
    cz <- .resolve_fit(tmpdf, fit, Mtype = Mtype, Atype = Atype)
    c <- cz$c
    z <- cz$z
    # plot
    graphics::plot(x = tmpdf[, Atype], y = tmpdf[, Mtype], xlab = Atype, ylab = Mtype, ...)
    # M = c*N^z. curve() follows the axes it is drawn on, so this stays correct
    # when the caller asks for log = "x" / "y" / "xy"
    if (!is.null(c) & !is.null(z)) {
        curve(c * x^z, add = TRUE, col = .anncol)
    }
    # the plotted points here ARE the theory rather than observed data
    .ann_legend("bottomright", list(...), "theory",
        eq = if (!is.null(c) & !is.null(z)) .eq_marsamp(c, z, Mtype = Mtype, Atype = Atype)
    )
    return(invisible())
}

# ---- internal helpers --------------------------------------------------------

# renders Mtype as its natural plotmath symbol: thetaw / thetapi as
# theta[w] / theta[pi] (rendered as the Greek letter with a subscript),
# anything else (M, E, ...) as its own name
.mtype_label <- function(Mtype) {
    switch(Mtype,
        thetaw = quote(theta[w]),
        thetapi = quote(theta[pi]),
        as.name(Mtype)
    )
}

# the fitted power-law equation as a plotmath expression, shared by the combined
# data/fit/theory legends of plot.marsamp and plot.martheory
.eq_marsamp <- function(c, z, Mtype = "M", Atype = "A") {
    bquote(.(.mtype_label(Mtype)) == .(round(c, 2)) * .(as.name(Atype))^.(round(z, 2)))
}

# the legend shared by every plot method: one key for the plotted data, plus
# optional keys for the fitted equation (`eq`) and for an overlaid curve
# (`overlaylab`, "theory" or "expected"). the data key's colour and symbol are
# read back out of the plot call's arguments so the key matches the points
# actually drawn -- with [[ ]], since $ partial-matches and would let col.axis
# stand in for col. `dataline` is the data key's line type: NA for a scatter,
# 1 for a plot whose points are joined (plot.sfs, type = "b"). the eq and
# overlay keys are always line-only. a lone data key is not worth a legend, so
# nothing is drawn in that case.
.ann_legend <- function(location, plotargs, datalab, eq = NULL, overlaylab = NULL, dataline = NA) {
    leglab <- list(datalab)
    legcol <- if (!is.null(plotargs[["col"]])) plotargs[["col"]] else graphics::par("col")
    legpch <- if (!is.null(plotargs[["pch"]])) plotargs[["pch"]] else graphics::par("pch")
    leglty <- dataline
    if (!is.null(eq)) {
        leglab <- c(leglab, eq)
        legcol <- c(legcol, .anncol)
        legpch <- c(legpch, NA)
        leglty <- c(leglty, 1)
    }
    if (!is.null(overlaylab)) {
        leglab <- c(leglab, overlaylab)
        legcol <- c(legcol, .theorycol)
        legpch <- c(legpch, NA)
        leglty <- c(leglty, 1)
    }
    if (length(leglab) > 1) {
        legend(location, legend = as.expression(leglab), col = legcol, pch = legpch, lty = leglty, bty = "n")
    }
    return(invisible())
}

# resolve the expected (neutral) sfs to overlay from a single `expected`
# argument: NULL (no overlay), an `sfs` object (used as-is), or a `genomaps`
# object (auto-generates one via expsfs(), matching `x`'s folding/zero-handling
# so the bins line up)
.resolve_expsfs <- function(x, expected) {
    if (is.null(expected)) {
        return(NULL)
    }
    if (inherits(expected, "sfs")) {
        return(expected)
    }
    if (inherits(expected, "genomaps")) {
        return(expsfs(expected, folded = attr(x, "folded"), nozero = attr(x, "nozero")))
    }
    stop("`expected` must be an `sfs` object, a `genomaps` object, or NULL")
}

# overlay the expected sfs on top of the observed sfs line/point plot of `x`,
# aligning bins by name in case `expected` uses a different bin set.
.overlay_expsfs <- function(bins, x, expected, col = .theorycol) {
    aligned <- as.vector(expected)[match(names(x), names(expected))]
    lines(bins, aligned, col = col, type = 'b')
    return(invisible())
}

# resolve the SAR power-law parameters (c, z) to plot, from a single `fit`
# argument: FALSE/NULL (no curve), TRUE (fit both via MARcalc), or a named
# c(c = , z = ) vector supplying one or both directly -- whichever of the two
# is missing or NA is still filled in via MARcalc.
.resolve_fit <- function(mardf, fit, Mtype, Atype) {
    if (isFALSE(fit) || is.null(fit)) {
        return(list(c = NULL, z = NULL))
    }
    if (isTRUE(fit)) {
        c <- NA_real_
        z <- NA_real_
    } else {
        stopifnot(
            "`fit` must be FALSE, TRUE, or a named c(c = , z = ) vector" =
                is.numeric(fit) && !is.null(names(fit)) && all(names(fit) %in% c("c", "z"))
        )
        c <- if ("c" %in% names(fit)) fit[["c"]] else NA_real_
        z <- if ("z" %in% names(fit)) fit[["z"]] else NA_real_
    }
    if (is.na(c) || is.na(z)) {
        mar <- MARcalc(mardf, Mtype = Mtype, Atype = Atype)
        marsum <- .marsummary(mar)
        if (is.na(c)) c <- marsum$c
        if (is.na(z)) z <- marsum$z
    }
    if (is.na(c)) c <- NULL
    if (is.na(z)) z <- NULL
    return(list(c = c, z = z))
}

# resolve the martheory object to overlay from a single `theory` argument:
# NULL (no overlay), a `martheory` object (used as-is), or a `genomaps` object
# (auto-generates one via MARtheory(), with the default maxind)
.resolve_theory <- function(theory) {
    if (is.null(theory)) {
        return(NULL)
    }
    if (inherits(theory, "martheory")) {
        return(theory)
    }
    if (inherits(theory, "genomaps")) {
        return(MARtheory(theory))
    }
    stop("`theory` must be a `martheory` object, a `genomaps` object, or NULL")
}

# overlay predicted M / thetaw / thetapi ~ N from a martheory object on top of
# a plot.marsamp plot. only Atype = "N" is comparable, since MARtheory() only
# predicts diversity as a function of the number of individuals
.overlay_theory_marsamp <- function(theory, Mtype, Atype, col = .theorycol) {
    if (Atype != "N" || !all(c(Atype, Mtype) %in% colnames(theory))) {
        warning("MARtheory overlay requires Atype = \"N\" and a matching Mtype column in `theory`; skipping overlay")
        return(invisible())
    }
    tdf <- theory[, c(Atype, Mtype)]
    tdf <- tdf[(tdf[, Mtype] > 0 & !is.na(tdf[, Mtype])), ]
    lines(x = tdf[, Atype], y = tdf[, Mtype], col = col)
    return(invisible())
}

# same, but rescaled to the % lost / % remained axes used by plot.marextinct.
# unlike marsamp/marextinct output (sorted with row 1 = full sample), martheory
# is always sorted with N ascending from 1, so the full-sample reference row is
# picked out by which.max() rather than assumed to be the first row
.overlay_theory_marextinct <- function(theory, Mtype, Atype, col = .theorycol) {
    if (Atype != "N" || !all(c(Atype, Mtype) %in% colnames(theory))) {
        warning("MARtheory overlay requires Atype = \"N\" and a matching Mtype column in `theory`; skipping overlay")
        return(invisible())
    }
    tdf <- theory[, c(Atype, Mtype)]
    tdf <- tdf[(tdf[, Mtype] > 0 & !is.na(tdf[, Mtype])), ]
    full <- tdf[which.max(tdf[, Atype]), ]
    a_per <- 1 - tdf[, Atype] / full[, Atype]
    m_per <- tdf[, Mtype] / full[, Mtype]
    lines(x = a_per, y = m_per, col = col)
    return(invisible())
}
