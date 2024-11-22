# a collection of plot methods
.anncol = "darkgray"
# RColorBrewer::brewer.pal(9, "Set1")
.catcol = c("#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

.ann_marsamp <- function(c, z, location) {
    equation <- bquote(M == .(round(c,2)) * A^.(round(z,2)))
    legend(location, legend = as.expression(equation), bty = "n", text.col = .anncol)
}

.ann_marextinct <- function(z, location) {
    equation <- bquote(m == 1 - (1-a)^.(round(z,2)))
    legend(location, legend = as.expression(equation), bty = "n", text.col = .anncol)
}

.ann_marsadsfs <- function(aa, ll, location) {
    legend(location, legend = paste0("AIC = ", round(aa, 2), "\nLL = ", round(ll, 2)), bty = "n")
}

# define the plotting method for marmaps
# methods(class = "marmaps") > [1] plot
#' Title
#'
#' @param marmaps
#'
#' @return
#' @export
#'
#' @examples
plot.marmaps <- function(obj) {
    old_par <- par(no.readonly = T)
    par(mar = c(5.1, 4.1, 4.1, 4.1))
    raster::plot(raster::extent(obj$samplemap), xlab = 'lon', ylab = 'lat')
    raster::plot(obj$samplemap, add = T, legend = F)
    points(obj$lonlat)
    raster::plot(obj$samplemap, add = T, legend.only = T, legend.mar = 3, legend.args = list(text = 'n'))
    par(old_par)
    return(invisible())
}

#' Title
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
plot.sfs <- function(obj, ...) {
    data <- as.vector(obj)
    bins <- as.integer(names(obj))
    graphics::barplot(data ~ bins, xlab = "Allele Count", ylab = "Number of Alleles", ...)
    return(invisible())
}

# this also works for classes belonging to marextinct (also a sampling process)
#' Title
#'
#' @param obj
#' @param c
#' @param z
#' @param Mtype
#' @param Atype
#' @param logscale
#'
#' @return
#' @export
#'
#' @examples
plot.marsamp <- function(obj, c = NULL, z = NULL, Mtype = .Mtype, Atype = .Atype, logscale = FALSE) {
    Mtype = match.arg(Mtype)
    Atype = match.arg(Atype)
    # remove NA or zero data
    tmpdf = obj[, c(Atype, Mtype)]
    tmpdf = tmpdf[(tmpdf[,Mtype] > 0 & !is.na(tmpdf[,Mtype])), ]
    # previously had a check for length(unique(tmpdf[,Mtype])) < 4
    if (nrow(tmpdf) == 0) {
        warning(paste0("MAR for Mtype = ", Mtype, ", Atype = ", Atype, " cannot be plotted"))
        graphics::plot.new()
    } else {
        # plot
        if (logscale) {
            graphics::plot(x = tmpdf[, Atype], y = tmpdf[, Mtype], log = 'xy', xlab = Atype, ylab = Mtype)
            # log(M) = log(c) + A*z
            if(!is.null(c) & !is.null(z)) {
                abline(a = c, b = z, col = .anncol)
                .ann_marsamp(c, z, location = "topright")
            }
        } else{
            graphics::plot(x = tmpdf[, Atype], y = tmpdf[, Mtype], xlab = Atype, ylab = Mtype)
            # M = c*A^z
            if (!is.null(c) & !is.null(z)) {
                curve(c * x^z, add = TRUE, col = .anncol)
                .ann_marsamp(c, z, location = "topright")
            }
        }
    }
    return(invisible())
}

#' Title
#'
#' @param obj
#' @param z
#' @param Mtype
#' @param Atype
#'
#' @return
#' @export
#'
#' @examples
plot.marextinct <- function(obj, z = NULL, Mtype = .Mtype, Atype = .Atype) {
    Mtype = match.arg(Mtype)
    Atype = match.arg(Atype)
    # remove NA or zero data
    tmpdf = obj[, c(Atype, Mtype)]
    tmpdf = tmpdf[(tmpdf[,Mtype] > 0 & !is.na(tmpdf[,Mtype])), ]
    # generate area loss data
    a_per = 1 - tmpdf[, Atype]/max(tmpdf[, Atype])
    m_per = 1- tmpdf[, Mtype]/max(tmpdf[, Mtype])
    # previously had a check for length(unique(tmpdf[,Mtype])) < 4
    if (nrow(tmpdf) == 0) {
        stop(paste0("MAR for Mtype = ", Mtype, ", Atype = ", Atype, " cannot be plotted"))
    }
    # plot
    # m_per = 1 - (1-a_per)^z
    graphics::plot(x = a_per, y = m_per, xlab = paste0("% of ", Atype, " lost"), ylab = paste0("% of ", Mtype, " lost"))
    if (!is.null(z)) {
        curve(1-(1-x)^z, add = TRUE, col = .anncol)
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

