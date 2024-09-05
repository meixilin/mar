allowed_schemes = c("random", "inwards", "outwards", "southnorth", "northsouth")

# match with old MARsampling
# TODO: Add not gridded sampling
MARsampling <- function(gm, scheme = allowed_schemes, nrep = 10, gridded = TRUE, animate = FALSE, quorum = FALSE, myseed = NULL, debug = FALSE) {
    # set seed if specified
    if (!is.null(myseed)) {
        set.seed(myseed)
    }
    # calculate and store raster area in the given gm$samplemap
    gmarea = areaofraster(gm$samplemap)
    # the x and y number of cells in gm$samplemap
    latrange <- dim(gm$samplemap)[1] # y is Row is Lat. Selected by r1, r2.
    lonrange <- dim(gm$samplemap)[2] # x is Col is Lon. Selected by c1, c2.
    minrange <- min(latrange, lonrange) # the maximum size box can become
    # the point where most samples are available (for inwards / outwards sampling)
    maxids <- raster::which.max(gm$samplemap)
    if(length(maxids) > 1) {warning('More than one cell with maximum samples')}
    r0c0 <- rowColFromCell(gm$samplemap, maxids[1])
    # match schemes (default to random)
    scheme = match.arg(scheme)
    # find right stepsize
    mystep = ifelse(minrange > 100, ceiling(minrange/100), 1)
    sidesize = seq(1, minrange, by = mystep)
    # differences btw different schemes are in the bounding boxes selected for diversity calculations
    bboxlist = vector('list', length = nrep)
    for (ii in 1:nrep) {
        bboxlist[[ii]] = switch(scheme,
                                random = .random_MARsampling(latrange, lonrange, sidesize),
                                inwards = .point_MARsampling(latrange, lonrange, sidesize, r0c0),
                                outwards = .point_MARsampling(latrange, lonrange, sidesize, r0c0),
                                southnorth = .pole_MARsampling(latrange, lonrange, sidesize, from = 'S'),
                                northsouth = .pole_MARsampling(latrange, lonrange, sidesize, from = 'N')
        )
    }
    # if need to plot
    if (animate) {
        lapply(bboxlist, .animate_MARsampling, gm = gm)
    }
    bboxlist = unlist(bboxlist, recursive = FALSE)
    # if scheme is inwards, reverse the bounding box
    revbbox = FALSE
    if (scheme == "inwards") revbbox = TRUE
    # calculate area and genetic diversity in each bounding boxes
    outlist = lapply(bboxlist, mar:::.mutdiv.gridded, gm = gm, gmarea = gmarea, revbbox = revbbox)
    outdf = do.call(dplyr::bind_rows, lapply(outlist, as.data.frame))
    # return bounding boxes as well
    outdf$extent = paste0(unlist(lapply(bboxlist, paste0, collapse = ';')))
    if (revbbox) {outdf$extent = paste0('-', outdf$extent)} # mark reverse selections
    # set outdf as a marsamp class
    class(outdf) <- "marsamp" # marsampling output class
    attr(outdf, 'scheme') <- scheme
    return(outdf)
}

# sampling methods
# old method but fix the minimum bug.
.random_MARsampling <- function(latrange, lonrange, sidesize, quorum = FALSE) {
    bblist = lapply(sidesize, function(ss) {
        r1 <- sample(1:(latrange - ss + 1), size = 1); r2 <- r1 + ss - 1
        c1 <- sample(1:(lonrange - ss + 1), size = 1); c2 <- c1 + ss - 1
        return(c(r1, r2, c1, c2))
    })
    return(bblist)
}

# southnorth / northsouth sampling
.pole_MARsampling <- function(latrange, lonrange, sidesize, from = c('N', 'S')) {
    bblist = lapply(sidesize, function(ss) {
        # TODO: assumes row = 1, col = 1 of raster is northwest corner. Add row, go south. Add col, go east.
        # use a geometric distribution to scale probability
        myprob = dgeom(0:(latrange - ss), prob = 0.5)
        myprob = myprob/sum(myprob)
        if (from == 'S') {myprob = rev(myprob)}
        # same selection of r1 c1, except that r1 is sampled with a probability highest in the south/north
        r1 <- sample(1:(latrange - ss + 1), size = 1, prob = myprob); r2 <- r1 + ss - 1
        c1 <- sample(1:(lonrange - ss + 1), size = 1); c2 <- c1 + ss - 1
        return(c(r1, r2, c1, c2))
    })
    return(bblist)
}

# inwards / outwards sampling
.point_MARsampling <- function(latrange, lonrange, sidesize, r0c0) {
    bblist = lapply(sidesize, function(ss) {
        # TODO: assumes row = 1, col = 1 of raster is northwest corner. Add row, go south. Add col, go east.
        # use a geometric distribution to scale probability
        # maximum probability at r0+1/2 and c0+1/2.
        rvars = 1:(latrange - ss + 1)
        cvars = 1:(lonrange - ss + 1)
        rprob = dgeom(abs(rvars - (r0c0[1,1] + 0.5 - 0.5*ss))*2, prob = 0.5)
        cprob = dgeom(abs(cvars - (r0c0[1,2] + 0.5 - 0.5*ss))*2, prob = 0.5)
        rprob = rprob/sum(rprob)
        cprob = cprob/sum(cprob)
        # same selection of r1 c1, except that r1 is sampled with a probability highest in the south/north
        r1 <- sample(rvars, size = 1, prob = rprob); r2 <- r1 + ss - 1
        c1 <- sample(cvars, size = 1, prob = cprob); c2 <- c1 + ss - 1
        return(c(r1, r2, c1, c2))
    })
    return(bblist)
}

# animate the sampling results
.animate_MARsampling <- function(gm, bblist, pause = 0.2) {
    grDevices::dev.flush()
    plot(gm$samplemap)
    for (ii in seq_along(bblist)) {
        plot(rowcol_extent(gm, bblist[[ii]]), add = T, col = 'black')
        Sys.sleep(pause)
    }
}

