# not a reverse function of MARsampling. Here the operation is on a cell level but MARsampling circles the grid with boxes

MARextinction <- function(gm, scheme = .MARsampling_schemes, nrep = 10, xfrac = 0.01, animate = FALSE, myseed = NULL) {
    # same as MARsampling ------------------------------------------------------
    # set seed if specified
    if (!is.null(myseed)) {
        set.seed(myseed)
    }
    # match schemes (default to random)
    scheme = match.arg(scheme)
    # calculate and store raster area in the given gm$maps$samplemap
    gmarea = areaofraster(gm$maps$samplemap)
    # the point where most samples are available (for inwards / outwards sampling)
    maxids <- raster::which.max(gm$maps$samplemap)
    if (length(maxids) > 1) {
        warning('More than one cell with maximum samples')
    }
    r0c0 <- raster::rowColFromCell(gm$maps$samplemap, maxids[1])
    # End same as MARsampling --------------------------------------------------
    extlist <- .extlist_sample(gm, xfrac, scheme, nrep, r0c0)

    # if need to plot
    if (animate) {
        lapply(extlist, .animate_MARextinction, gm = gm)
    }

    # calculate area and genetic diversity in each extant cell list
    outlist <- lapply(seq_along(extlist), function(ii) {
        outl <- lapply(extlist[[ii]], mutdiv.cells, gm = gm, gmarea = gmarea)
        out <- do.call(rbind, lapply(outl, as.data.frame))
        # append end theta (zero in all)
        out[nrow(out)+1, ] <- rep(0, ncol(out))
        out$repid <- ii # replicate id
        return(out)
    })
    outdf <- do.call(rbind, outlist)

    # set outdf as a marsamp class
    class(outdf) <- c("data.frame", "marextinct") # marextinction output class
    attr(outdf, 'scheme') <- scheme
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
            myprob <- myprob / sum(myprob)
        }
    }
    # add names if myprob is not NULL
    if (!is.null(myprob)) {
        names(myprob) <- gridpresent
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
    gridrowcol <- raster::rowColFromCell(gm$maps$samplemap, gridpresent)
    # find right stepsize
    mystep <- ifelse(length(gridpresent) > 100, ceiling(length(gridpresent) * xfrac), 1)
    rvars <- gridrowcol[,'row']
    cvars <- gridrowcol[,'col']

    # Calculate probability of all grids (rescale at each step). synonymous to the rcprob
    rcprob <- switch(
        scheme,
        random = list(NULL, NULL),
        # no prob
        inwards = lapply(.point_prob(rvars, cvars, r0c0, ss=1), function(x) 1-x),
        outwards = .point_prob(rvars, cvars, r0c0, ss=1),
        southnorth = .pole_prob(rvars, from = 'S'),
        northsouth = .pole_prob(rvars, from = 'N')
    )
    myprob <- .rcprob2myprob(rcprob, gridpresent)
    extlist <- lapply(1:nrep, function(ii) .extsample(gridpresent, myprob, mystep))
    return(extlist)
}

.extsample <- function(gridpresent, myprob, mystep) {
    # Create a list to store grids that remain after each extinction step
    extl <- vector('list', length = ceiling(length(gridpresent) / mystep))
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
    tempmaps <- gm$maps
    tempids <- 1:raster::ncell(tempmaps$samplemap)
    for (ii in seq_along(extl)) {
        tempmaps$samplemap[setdiff(tempids, extl[[ii]])] <- NA
        tempmaps$lonlat = gm$maps$lonlat[gm$maps$cellid %in% extl[[ii]], ]
        if (!is.matrix(tempmaps$lonlat)) {
            tempmaps$lonlat = matrix(tempmaps$lonlat, ncol = 2)
        }
        plot(tempmaps)
        Sys.sleep(pause)
    }
    return(invisible())
}
