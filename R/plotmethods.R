# a collection of plot methods
.anncol = "darkgray"
# RColorBrewer::brewer.pal(9, "Set1")
.catcol = c("#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

.ann_marsamp <- function(c, z, location) {
    equation <- bquote(M == .(round(c,2)) * A^.(round(z,2)))
    legend(location, legend = as.expression(equation), bty = "n", text.col = .anncol)
}

.ann_marextinct <- function(z, location) {
    equation <- bquote(m == (1-a)^.(round(z,2)))
    legend(location, legend = as.expression(equation), bty = "n", text.col = .anncol)
}

.ann_marsadsfs <- function(aa, ll, location) {
    legend(location, legend = paste0("AIC = ", round(aa, 2), "\nLL = ", round(ll, 2)), bty = "n")
}

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
    # old_par <- par(no.readonly = T)
    # par(mar = c(5.1, 4.1, 4.1, 4.1))
    raster::plot(raster::extent(x$samplemap), xlab = 'lon', ylab = 'lat')
    raster::plot(x$samplemap, add = T, legend = F)
    points(x$lonlat)
    raster::plot(x$samplemap, add = T, legend.only = T, legend.mar = 3, legend.args = list(text = 'n'))
    # par(old_par)
    return(invisible())
}

#' Plot Method for Site Frequency Spectrum
#'
#' Creates a barplot of the site frequency spectrum
#'
#' @param x An object of class "sfs" containing site frequency spectrum data
#' @param ... Additional arguments passed to barplot
#'
#' @return Invisibly returns NULL
#' @export
#'
#' @examples
#' sfs_object <- sfs(AC=c(1,1,0,2,0,1,1,0,0,2,30), N=50, ploidy=2)
#' plot(sfs_object)
plot.sfs <- function(x, ...) {
    data <- as.vector(x)
    bins <- as.integer(names(x))
    graphics::barplot(data ~ bins, xlab = "Allele Count", ylab = "Number of Alleles", ...)
    return(invisible())
}

# this also works for classes belonging to marextinct (also a sampling process)
#' Plot Method for Mutations-Area Relationship sample data
#'
#' Creates a plot of genetic diversity against area, optionally with fitted power law
#'
#' @param x A data frame containing sampling data
#' @param c Intercept parameter of fitted power law curve
#' @param z Slope parameter of fitted power law curve
#' @param Mtype Type of genetic diversity metric to plot
#' @param Atype Type of area metric to plot
#' @param logscale Logical, whether to plot on log-log scale
#' @param ... Additional arguments passed to plot
#'
#' @return Invisibly returns NULL
#' @export
#'
#' @examples
#' plot(marsamp_object, c=0.5, z=0.25)
plot.marsamp <- function(x, c = NULL, z = NULL, Mtype = .Mtype, Atype = .Atype, logscale = FALSE, ...) {
    Mtype = match.arg(Mtype)
    Atype = match.arg(Atype)
    tmpdf = x[, c(Atype, Mtype)]
    tmpdf = tmpdf[(tmpdf[,Mtype] > 0 & !is.na(tmpdf[,Mtype])), ]
    # previously had a check for length(unique(tmpdf[,Mtype])) < 4
    if (nrow(tmpdf) == 0) {
        warning(paste0("MAR for Mtype = ", Mtype, ", Atype = ", Atype, " cannot be plotted"))
        graphics::plot.new()
    } else {
        # plot
        if (logscale) {
            graphics::plot(x = tmpdf[, Atype], y = tmpdf[, Mtype], log = 'xy', xlab = Atype, ylab = Mtype, ...)
            # log(M) = log(c) + A*z
            if(!is.null(c) & !is.null(z)) {
                abline(a = c, b = z, col = .anncol)
                .ann_marsamp(c, z, location = "topright")
            }
        } else{
            graphics::plot(x = tmpdf[, Atype], y = tmpdf[, Mtype], xlab = Atype, ylab = Mtype, ...)
            # M = c*A^z
            if (!is.null(c) & !is.null(z)) {
                curve(c * x^z, add = TRUE, col = .anncol)
                .ann_marsamp(c, z, location = "topright")
            }
        }
    }
    return(invisible())
}

#' Plot Method for Mutations Extinction Curves
#'
#' Creates a plot showing relationship between area loss and genetic diversity loss
#'
#' @param x A data frame containing extinction data
#' @param z Fitted extinction curve parameter
#' @param Mtype Type of genetic diversity metric to plot
#' @param Atype Type of area metric to plot
#'
#' @return Invisibly returns NULL
#' @export
#'
#' @examples
#' plot(marextinct_object, z=0.25)
plot.marextinct <- function(x, z = NULL, Mtype = .Mtype, Atype = .Atype) {
    Mtype = match.arg(Mtype)
    Atype = match.arg(Atype)
    # remove NA or zero data
    tmpdf = x[, c(Atype, Mtype)]
    tmpdf = tmpdf[(tmpdf[,Mtype] > 0 & !is.na(tmpdf[,Mtype])), ]
    # generate area loss data
    a_per = 1 - tmpdf[, Atype]/max(tmpdf[, Atype])
    m_per = tmpdf[, Mtype]/max(tmpdf[, Mtype])
    # previously had a check for length(unique(tmpdf[,Mtype])) < 4
    if (nrow(tmpdf) == 0) {
        stop(paste0("MAR for Mtype = ", Mtype, ", Atype = ", Atype, " cannot be plotted"))
    }
    # plot
    # m_per = 1 - (1-a_per)^z
    graphics::plot(x = a_per, y = m_per, xlab = paste0("% of ", Atype, " lost"), ylab = paste0("% of ", Mtype, " remained"))
    if (!is.null(z)) {
        curve((1-x)^z, add = TRUE, col = .anncol)
        .ann_marextinct(z, location = "topright")
    }
    return(invisible())
}


.pipe_plot.marsadsfs <- function(obj, AICtabs) {
    old_par <- par(no.readonly = T)
    par(mfrow = c(ceiling(length(obj$sfs) / 2), 2), mar = c(5.1, 4.1, 2.1, 2.1))
    for (ii in seq_along(obj$sfs)) {
        mname = names(obj$sfs)[ii]
        plot.sfs(obj$sfs[[ii]], col = .catcol[ii], border = NA, main = mname)
        .ann_marsadsfs(aa = AICtabs$AIC[attr(AICtabs, "row.names") == mname],
                       ll = obj$statdf[obj$statdf$model == mname, 'logLik'], location = "topright")
    }
    par(old_par)
}


.pipe_plot.marsamp <- function(mardf, mar) {
    old_par <- par(no.readonly = T)
    par(mfcol = c(length(unique(mar$M)), length(unique(mar$A))), mar = c(5.1, 4.1, 2.1, 2.1))
    for (ii in seq_along(mar)) {
        plot.marsamp(mardf, c = mar$c[ii], z = mar$z[ii], Mtype = mar$M[ii], Atype = mar$A[ii])
    }
    par(old_par)
    return(invisible())
}

.pipe_plot.marextinct <- function(extdf, ext) {
    old_par <- par(no.readonly = T)
    par(mfcol = c(length(unique(ext$M)), length(unique(ext$A))), mar = c(5.1, 4.1, 2.1, 2.1))
    for (ii in seq_along(ext)) {
        plot.marextinct(extdf, z = ext$z[ii], Mtype = ext$M[ii], Atype = ext$A[ii])
    }
    par(old_par)
    return(invisible())
}

