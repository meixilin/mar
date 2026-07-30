.MARsampling_schemes <- c("random", "inwards", "outwards", "southnorth", "northsouth")

#' MAR sampling wrapper function
#'
#' @param gm a [genomaps] object created by [genomaps()]
#' @param scheme sampling schemes for spatial data. allowed are `r toString(.MARsampling_schemes)`.
#' @param nrep number of replicates.
#' @param xfrac fraction of the range to use for the step size in the sampling scheme.
#' @param quorum require all sampling grid to have samples. default is FALSE.
#' @param animate play an animation of the sampling boxes. default is FALSE.
#' @param myseed set seed for reproducibility. default is NULL.
#'
#' @return a `marsamp` object. consist of a data frame.
#' @export
#'
MARsampling <-
    function(gm,
             scheme = .MARsampling_schemes,
             nrep = 10,
             xfrac = 0.01,
             quorum = TRUE,
             animate = FALSE,
             myseed = NULL) {
        # set seed if specified
        if (!is.null(myseed)) {
            set.seed(myseed)
        }
        # match schemes (default to random)
        scheme <- match.arg(scheme)
        # if scheme is inwards, reverse the bounding box
        revbbox <- ifelse(scheme == "inwards", TRUE, FALSE)
        # unwrap the samplemap raster once for the terra operations below
        sm <- .get_samplemap(gm$maps)
        # calculate and store raster area in the given samplemap
        gmarea <- .areaofraster(sm)
        # the x and y number of cells in samplemap
        # y is Row is Lat. Selected by r1, r2.
        # x is Col is Lon. Selected by c1, c2.
        latrange <- dim(sm)[1]
        lonrange <- dim(sm)[2]
        # the maximum size box can become
        minrange <- min(latrange, lonrange)
        # the point where most samples are available (for inwards / outwards sampling)
        # TODO: user specified central point
        sm <- .as_spatraster(sm)
        maxrc <- terra::where.max(sm) # can have multiple rows when tied
        r0c0 <- terra::rowColFromCell(sm, maxrc[1, 'cell'])
        # calculate and store non-empty cell once, as a summed-area table for O(1) box occupancy queries
        pres <- !is.na(terra::as.matrix(sm, wide = TRUE))
        cumocc <- apply(pres, 2, cumsum)
        cumocc <- t(apply(cumocc, 1, cumsum))
        cumocc <- rbind(0L, cbind(0L, cumocc))
        # find right stepsize
        mystep <- ifelse(minrange > 100, ceiling(minrange * xfrac), 1)
        sidesize <- seq(1, minrange, by = mystep)
        # differences btw different schemes are in the bounding boxes selected for diversity calculations
        # differences in bounding boxes are in the sample probability settings
        bboxlist <- lapply(
            sidesize,
            .bblist_sample,
            scheme = scheme,
            nrep = nrep,
            quorum = quorum,
            latrange = latrange,
            lonrange = lonrange,
            r0c0 = r0c0,
            revbbox = revbbox,
            cumocc = cumocc
        )

        # if need to plot
        if (animate) {
            lapply(bboxlist, .animate_MARsampling, gm = gm)
        }
        bboxlist <- unlist(bboxlist, recursive = FALSE)
        # calculate area and genetic diversity in each bounding boxes
        outlist <- lapply(
            bboxlist,
            mutdiv.gridded,
            gm = gm,
            gmarea = gmarea,
            sm = sm,
            revbbox = revbbox
        )
        # use rbind to avoid importing dplyr
        outdf <- do.call(rbind, lapply(outlist, as.data.frame, stringsAsFactors = FALSE))
        # return bounding boxes as well
        outdf$extent <- unlist(lapply(bboxlist, paste0, collapse = ";"))
        if (revbbox) {
            outdf$extent <- paste0("-", outdf$extent)
        } # mark reverse selections
        # set outdf as a marsamp class
        class(outdf) <- c("marsamp", class(outdf)) # marsampling output class
        attr(outdf, "scheme") <- scheme
        return(outdf)
    }

# scheme based probfunc
.prob_sample <- function(xx) {
    if (length(xx) == 0) {
        return(xx)
    }
    pp <- stats::dgeom(xx * 2, prob = 0.5)
    if (all(pp == 0) || sum(pp, na.rm = TRUE) == 0) {
        # fallback to uniform sampling
        return(rep(1 / length(xx), length(xx)))
    }

    pp <- pp / sum(pp)

    # handle any NaN/Inf just in case
    if (any(!is.finite(pp))) {
        return(rep(1 / length(xx), length(xx)))
    }
    return(pp / sum(pp))
}

# southnorth / northsouth sample
.pole_prob <- function(rvars, from = c("N", "S")) {
    rxx <- rvars - 1
    rprob <- .prob_sample(rxx)
    if (from == "S") {
        rprob <- rev(rprob)
    }
    return(list(rprob, NULL))
}

# inwards / outwards sample
.point_prob <- function(rvars, cvars, r0c0, ss) {
    stopifnot(all(dim(r0c0) == c(1, 2)))
    rxx <- abs(rvars - (r0c0[1, 1] + 0.5 - 0.5 * ss))
    cxx <- abs(cvars - (r0c0[1, 2] + 0.5 - 0.5 * ss))
    rprob <- .prob_sample(rxx)
    cprob <- .prob_sample(cxx)
    return(list(rprob, cprob))
}

# sample bounding boxes
# when quorum, restrict sampling jointly to boxes with valid == TRUE (matrix lookup from the summed-area table)
.bbsample <- function(ss, nrep, rvars, cvars, rprob, cprob, quorum, valid) {
    if (length(rvars) == 0 || length(cvars) == 0) {
        return(list())
    }

    if (!is.null(rprob) && length(rprob) != length(rvars)) {
        rprob <- NULL
    }
    if (!is.null(cprob) && length(cprob) != length(cvars)) {
        cprob <- NULL
    }

    if (quorum) {
        if (!any(valid)) {
            warning("Cannot fulfill quorum, try set quorum = FALSE")
            return(list())
        }
        rp <- if (is.null(rprob)) rep(1, length(rvars)) else rprob
        cp <- if (is.null(cprob)) rep(1, length(cvars)) else cprob
        probmat <- outer(rp, cp) * valid
        idx <- sample.int(length(rvars) * length(cvars), size = nrep, prob = as.vector(probmat), replace = TRUE)
        rc <- arrayInd(idx, dim(probmat))
        r1 <- rvars[rc[, 1]]
        c1 <- cvars[rc[, 2]]
    } else {
        # allow replacement, so that the probability is respected
        r1 <- sample(rvars, size = nrep, prob = rprob, replace = TRUE)
        c1 <- sample(cvars, size = nrep, prob = cprob, replace = TRUE)
    }
    r2 <- r1 + ss - 1
    c2 <- c1 + ss - 1

    bblist <-
        lapply(1:nrep, function(ii) {
            c(r1[ii], r2[ii], c1[ii], c2[ii])
        })
    return(bblist)
}

# core sampling function
.bblist_sample <-
    function(ss,
             scheme,
             nrep,
             quorum,
             latrange,
             lonrange,
             r0c0,
             revbbox,
             cumocc) {
        # at this sidesize, the available row and column numbers
        # TODO: assumes row = 1, col = 1 of raster is northwest corner.
        # Add row, go south. Add col, go east.
        if ((latrange - ss + 1) <= 0 || (lonrange - ss + 1) <= 0) {
            return(list()) # skip invalid box sizes
        }
        rvars <- 1:(latrange - ss + 1)
        cvars <- 1:(lonrange - ss + 1)
        # generate rprob and cprob given the scheme specified and current rvars / cvars
        rcprob <- switch(scheme,
            random = list(NULL, NULL),
            # no prob
            inwards = .point_prob(rvars, cvars, r0c0, ss),
            outwards = .point_prob(rvars, cvars, r0c0, ss),
            southnorth = .pole_prob(rvars, from = "S"),
            northsouth = .pole_prob(rvars, from = "N")
        )
        # matrix lookup: occupancy count for every (r1, c1) box of this size, via the summed-area table cumocc
        valid <- NULL
        if (quorum) {
            countmat <- cumocc[rvars + ss, cvars + ss] - cumocc[rvars, cvars + ss] - cumocc[rvars + ss, cvars] + cumocc[rvars, cvars]
            if (revbbox) {
                countmat <- cumocc[nrow(cumocc), ncol(cumocc)] - countmat
            }
            valid <- countmat > 0
        }
        # sample bounding boxes; when quorum, sampling is restricted to valid boxes jointly
        bblist <- .bbsample(ss, nrep, rvars, cvars, rcprob[[1]], rcprob[[2]], quorum, valid)
        # # check for duplicates (some sampling method generate duplicates)
        # if (any(duplicated(bblist)))
        #     warning('Sampling generated duplicated bounding boxes')
        return(bblist)
    }

# animate the sampling results
.animate_MARsampling <- function(gm, bblist, pause = 0.2) {
    grDevices::dev.flush()
    plot.marmaps(gm$maps)
    sm <- .get_samplemap(gm$maps)
    for (ii in seq_along(bblist)) {
        terra::plot(terra::ext(sm[bblist[[ii]][1:2], bblist[[ii]][3:4], drop = FALSE]), add = T, legend = FALSE)
        Sys.sleep(pause)
    }
}
