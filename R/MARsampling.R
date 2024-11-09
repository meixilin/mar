.MARsampling_schemes = c("random", "inwards", "outwards", "southnorth", "northsouth")

# match with old MARsampling
# TODO: Add not gridded sampling
#' MAR sampling wrapper function
#'
#' @param gm a [genomaps] object created by [genomaps()]
#' @param scheme sampling schemes for spatial data. allowed are \code{\link{.MARsampling_schemes}}.
#' @param nrep number of replicates.
#' @param quorum require all sampling grid to have samples. default is FALSE.
#' @param animate play an animation of the sampling boxes. default is FALSE.
#' @param myseed set seed for reproducibility. default is NULL.
#'
#' @return a [marsamp] object. consist of a data frame.
#' @export
#'
#' @examples
MARsampling <-
    function(gm,
             scheme = .MARsampling_schemes,
             nrep = 10,
             xfrac = 0.01,
             quorum = FALSE,
             animate = FALSE,
             myseed = NULL) {
        # set seed if specified
        if (!is.null(myseed)) {
            set.seed(myseed)
        }
        # match schemes (default to random)
        scheme = match.arg(scheme)
        # if scheme is inwards, reverse the bounding box
        revbbox = ifelse(scheme == "inwards", TRUE, FALSE)
        # calculate and store raster area in the given gm$maps$samplemap
        gmarea = areaofraster(gm$maps$samplemap)
        # the x and y number of cells in gm$maps$samplemap
        # y is Row is Lat. Selected by r1, r2.
        # x is Col is Lon. Selected by c1, c2.
        latrange <- dim(gm$maps$samplemap)[1]
        lonrange <- dim(gm$maps$samplemap)[2]
        # the maximum size box can become
        minrange <- min(latrange, lonrange)
        # the point where most samples are available (for inwards / outwards sampling)
        maxids <- raster::which.max(gm$maps$samplemap)
        if (length(maxids) > 1) {
            warning('More than one cell with maximum samples')
        }
        r0c0 <- raster::rowColFromCell(gm$maps$samplemap, maxids[1])
        # find right stepsize
        mystep = ifelse(minrange > 100, ceiling(minrange * xfrac), 1)
        sidesize = seq(1, minrange, by = mystep)
        # differences btw different schemes are in the bounding boxes selected for diversity calculations
        # differences in bounding boxes are in the sample probability settings
        bboxlist = lapply(
            sidesize,
            .bblist_sample,
            gm = gm,
            scheme = scheme,
            nrep = nrep,
            quorum = quorum,
            latrange = latrange,
            lonrange = lonrange,
            r0c0 = r0c0,
            revbbox = revbbox
        )

        # if need to plot
        if (animate) {
            lapply(bboxlist, .animate_MARsampling, gm = gm)
        }
        bboxlist = unlist(bboxlist, recursive = FALSE)
        # calculate area and genetic diversity in each bounding boxes
        outlist = lapply(
            bboxlist,
            mutdiv.gridded,
            gm = gm,
            gmarea = gmarea,
            revbbox = revbbox
        )
        outdf = do.call(rbind, lapply(outlist, as.data.frame))
        # return bounding boxes as well
        outdf$extent = unlist(lapply(bboxlist, paste0, collapse = ';'))
        if (revbbox) {
            outdf$extent = paste0('-', outdf$extent)
        } # mark reverse selections
        # set outdf as a marsamp class
        class(outdf) <- c("data.frame", "marsamp") # marsampling output class
        attr(outdf, 'scheme') <- scheme
        return(outdf)
    }

# scheme based probfunc
.prob_sample <- function(xx) {
    pp = dgeom(xx * 2, prob = 0.5)
    return(pp / sum(pp))
}

# southnorth / northsouth sample
.pole_prob <- function(rvars, from = c('N', 'S')) {
    rxx = rvars - 1
    rprob = .prob_sample(rxx)
    if (from == 'S') {
        rprob = rev(rprob)
    }
    return(list(rprob, NULL))
}

# inwards / outwards sample
.point_prob <- function(rvars, cvars, r0c0, ss) {
    stopifnot(all(dim(r0c0) == c(1, 2)))
    rxx = abs(rvars - (r0c0[1, 1] + 0.5 - 0.5 * ss))
    cxx = abs(cvars - (r0c0[1, 2] + 0.5 - 0.5 * ss))
    rprob = .prob_sample(rxx)
    cprob = .prob_sample(cxx)
    return(list(rprob, cprob))
}

# sample bounding boxes
.bbsample <- function(ss, nrep, rvars, cvars, rprob, cprob) {
    # allow replacement, so that the probability is respected
    r1 <- sample(rvars, size = nrep, prob = rprob, replace = TRUE)
    c1 <- sample(cvars, size = nrep, prob = cprob, replace = TRUE)
    r2 <- r1 + ss - 1
    c2 <- c1 + ss - 1

    bblist <-
        lapply(1:nrep, function(ii)
            c(r1[ii], r2[ii], c1[ii], c2[ii]))
    return(bblist)
}

# core sampling function
.bblist_sample <-
    function(ss,
             gm,
             scheme,
             nrep,
             quorum,
             latrange,
             lonrange,
             r0c0,
             revbbox) {
        # at this sidesize, the available row and column numbers
        # TODO: assumes row = 1, col = 1 of raster is northwest corner.
        # Add row, go south. Add col, go east.
        rvars = 1:(latrange - ss + 1)
        cvars = 1:(lonrange - ss + 1)
        # generate rprob and cprob given the scheme specified and current rvars / cvars
        rcprob <- switch(
            scheme,
            random = list(NULL, NULL),
            # no prob
            inwards = .point_prob(rvars, cvars, r0c0, ss),
            outwards = .point_prob(rvars, cvars, r0c0, ss),
            southnorth = .pole_prob(rvars, from = 'S'),
            northsouth = .pole_prob(rvars, from = 'N')
        )
        # sample bounding boxes
        # run sampling
        bblist <- .bbsample(ss, nrep, rvars, cvars, rcprob[[1]], rcprob[[2]])
        # if need to have samples in all bounding boxes, rejection sampling
        if (quorum) {
            ncells <-
                sapply(lapply(
                    bblist,
                    rowcol_cellid,
                    mm = gm$maps,
                    revbbox = revbbox
                ),
                length)
            ntry <- 0
            nnrep = min(nrep * 5, min(latrange, lonrange) - 1)
            while (any(ncells == 0) & ntry < 50) {
                tbblist <-
                    .bbsample(ss, nrep = nnrep, rvars, cvars, rcprob[[1]], rcprob[[2]])
                tncells <-
                    sapply(lapply(
                        tbblist,
                        rowcol_cellid,
                        mm = gm$maps,
                        revbbox = revbbox
                    ),
                    length)
                nnewbb <- min(sum(tncells > 0), sum(ncells == 0))
                bbids <- which(ncells == 0)[nnewbb]
                bblist[bbids] <- tbblist[tncells > 0][nnewbb]
                ncells <-
                    sapply(lapply(
                        bblist,
                        rowcol_cellid,
                        mm = gm$maps,
                        revbbox = revbbox
                    ),
                    length)
                ntry = ntry + 1
            }
            if (ntry >= 50)
                warning('Cannot fulfill quorum, try set quorum = FALSE')
        }
        # # check for duplicates (some sampling method generate duplicates)
        # if (any(duplicated(bblist)))
        #     warning('Sampling generated duplicated bounding boxes')
        return(bblist)
    }

# animate the sampling results
.animate_MARsampling <- function(gm, bblist, pause = 0.2) {
    grDevices::dev.flush()
    plot(gm$maps)
    for (ii in seq_along(bblist)) {
        plot(rowcol_extent(gm$maps, bblist[[ii]]),
             add = T,
             col = 'black')
        Sys.sleep(pause)
    }
}
