#' Simulate the mutations-area relationship (MAR) extinction process
#'
#' This function performs a simulated extinction process on a map of genomic samples, similar to the
#' MARsampling function, but with the extinction happening at the cell level rather than
#' circling the grid with boxes. The function returns a data frame containing the area and
#' diversity metrics for each step of the extinction process.
#'
#' @param gm A [genomaps()] object created by [genomaps()].
#' @param scheme The order in which cells go extinct. Default is "random", allowed values are `r toString(.MARsampling_schemes)`.
#'   `inwards` / `outwards` remove cells starting from the cell holding the most samples, and `southnorth` / `northsouth`
#'   remove them starting from the southern or northern edge.
#' @param nrep The number of extinction replicates to perform. Default is 10.
#' @param xfrac The fraction of occupied cells removed at each extinction step.
#'   Rounded up, so at least one raster cell is removed per step. Default is 0.01.
#' @param animate If TRUE, the function will animate the extinction process and output to default plotting device. Default is FALSE.
#' @param myseed An optional seed value to ensure reproducibility. Implemented as `set.seed(myseed)`. Default is NULL.
#'
#' @return A `marextinct` object: a data frame with one row per extinction step
#'   and columns `N` (number of remaining samples), `M` (number of mutations),
#'   `E` (endemic mutations), `thetaw`, `thetapi`, `A` (remaining area), `extl`
#'   (the surviving cell ids) and `repid` (replicate id). Each replicate ends
#'   with a row of zeros, representing the loss of the whole range. The `scheme`
#'   used is stored as an attribute.
#' @seealso [MARsampling()] for the sampling counterpart and [plot.marextinct()]
#'   for plotting the extinction curve.
#' @export
#'
#' @examples
#' extdf <- MARextinction(gm1001g, nrep = 2, xfrac = 0.1, myseed = 42)
#' head(extdf[, c("N", "M", "thetaw", "thetapi", "A", "repid")])
#'
#' # diversity remaining against area lost, with the fitted extinction curve
#' plot(extdf, fit = TRUE)
#'
MARextinction <- function(gm, scheme = .MARsampling_schemes, nrep = 10, xfrac = 0.01, animate = FALSE, myseed = NULL) {
    # same as MARsampling ------------------------------------------------------
    # set seed if specified
    if (!is.null(myseed)) {
        set.seed(myseed)
    }
    # match schemes (default to random)
    scheme <- match.arg(scheme)
    # unwrap the samplemap raster once for the terra operations below
    sm <- .get_samplemap(gm$maps)
    # calculate and store raster area in the given samplemap
    gmarea <- .areaofraster(sm)
    # the point where most samples are available (for inwards / outwards sampling)
    # TODO: user specified central point
    maxrc <- terra::where.max(sm) # can have multiple rows when tied
    r0c0 <- terra::rowColFromCell(sm, maxrc[1, 'cell'])
    # End same as MARsampling --------------------------------------------------
    extlist <- .extlist_sample(gm, xfrac, scheme, nrep, r0c0)

    # if need to plot
    if (animate) {
        lapply(extlist, .animate_MARextinction, gm = gm)
    }

    # calculate area and genetic diversity in each extant cell list
    outlist <- lapply(seq_along(extlist), function(ii) {
        outl <- lapply(extlist[[ii]], mutdiv.cells, gm = gm, gmarea = gmarea)
        out <- do.call(rbind, lapply(outl, as.data.frame, stringsAsFactors = FALSE))
        # append end theta (zero in all)
        out[nrow(out) + 1, ] <- rep(0, ncol(out))
        # append extlist as well
        out$extl <- c(unlist(lapply(extlist[[ii]], paste0, collapse = ";")), "")
        out$repid <- ii # replicate id
        return(out)
    })
    outdf <- do.call(rbind, outlist)

    # set outdf as a marsamp class
    class(outdf) <- c("marextinct", class(outdf)) # marextinction output class
    attr(outdf, "scheme") <- scheme
    return(outdf)
}

.rcprob2myprob <- function(rcprob, gridpresent) {
    if (is.null(rcprob[[1]])) {
        myprob <- rcprob[[2]]
    } else {
        if (is.null(rcprob[[2]])) {
            myprob <- rcprob[[1]]
        } else {
            myprob <- rcprob[[1]] * rcprob[[2]]
            if (length(myprob) == 0 || sum(myprob) == 0) {
                return(NULL)
            }

            myprob <- myprob / sum(myprob)
        }
    }
    # add names if myprob is not NULL
    if (!is.null(myprob) && length(myprob) == length(gridpresent)) {
        names(myprob) <- gridpresent
    } else if (!is.null(myprob)) {
        # fallback: invalid probability vector
        return(NULL)
    }
    return(myprob)
}

.rescale_prob <- function(myprob) {
    if (!is.null(myprob)) {
        return(myprob / sum(myprob))
    } else {
        return(myprob)
    }
}

# core sampling function
.extlist_sample <- function(gm, xfrac, scheme, nrep, r0c0) {
    gridpresent <- sort(unique(gm$maps$cellid))
    gridrowcol <- terra::rowColFromCell(.get_samplemap(gm$maps), gridpresent)
    # find right stepsize
    mystep <- ceiling(length(gridpresent) * xfrac)
    rvars <- gridrowcol[, 1]
    cvars <- gridrowcol[, 2]

    # Calculate probability of all grids (rescale at each step). synonymous to the rcprob
    rcprob <- switch(scheme,
        random = list(NULL, NULL),
        # no prob
        inwards = lapply(.point_prob(rvars, cvars, r0c0, ss = 1), function(x) 1 - x),
        outwards = .point_prob(rvars, cvars, r0c0, ss = 1),
        southnorth = .pole_prob(rvars, from = "S"),
        northsouth = .pole_prob(rvars, from = "N")
    )
    myprob <- .rcprob2myprob(rcprob, gridpresent)
    extlist <- lapply(1:nrep, function(ii) .extsample(gridpresent, myprob, mystep))
    return(extlist)
}

.extsample <- function(gridpresent, myprob, mystep) {
    # Create a list to store grids that remain after each extinction step
    extl <- vector("list", length = ceiling(length(gridpresent) / mystep))
    extl[[1]] <- gridpresent

    # Simulate extinction process
    for (ii in 2:length(extl)) {
        if (length(gridpresent) <= mystep) {
            break
        }
        # Sample grids to become extinct
        toextinct <- base::sample(gridpresent, size = mystep, prob = myprob, replace = FALSE)

        # Update remaining grids
        gridpresent <- setdiff(gridpresent, toextinct)
        extl[[ii]] <- gridpresent
        myprob <- .rescale_prob(myprob[names(myprob) %in% gridpresent])
        # print(which.min(myprob))
        stopifnot(all(names(myprob) == gridpresent) | is.null(myprob))
    }
    # sanity check that the last on the list should be less than mystep
    stopifnot(length(extl[[length(extl)]]) <= mystep)

    return(extl)
}

.animate_MARextinction <- function(gm, extl, pause = 0.2) {
    grDevices::dev.flush()
    plot.marmaps(gm$maps)
    sm <- .get_samplemap(gm$maps)
    terra::values(sm) <- NA
    for (ii in seq_along(extl)) {
        sm[setdiff(gm$maps$cellid, extl[[ii]])] <- 1
        terra::plot(sm, add = T, col = "black", legend = FALSE)
        Sys.sleep(pause)
    }
}
