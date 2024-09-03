allowed_schemes = c("random", "inwards", "outwards", "southnorth", "northsouth")

# match with old MARsampling
# TODO: Add not gridded sampling
MARsampling <- function(gm, scheme = allowed_schemes, nrep = 10, gridded = TRUE, plot = FALSE, myseed = NULL, debug = FALSE) {
    # set seed if specified
    if (!is.null(myseed)) {
        set.seed(myseed)
    }
    # calculate and store raster area in the given gm$samplemap
    gmarea = areaofraster(gm$samplemap)
    # the x and y number of cells in gm$samplemap
    lonrange <- dim(gm$samplemap)[1]
    latrange <- dim(gm$samplemap)[2]
    minrange <- min(lonrange, latrange) # the maximum size box can become
    # match schemes (default to random)
    scheme = match.arg(scheme)
    # find right stepsize
    mystep = ifelse(minrange > 100, ceiling(minrange/100), 1)
    sidesize = seq(1, minrange, by = mystep)
    # differences btw different schemes are in the bounding boxes selected for diversity calculations
    bboxlist = vector('list', length = nrep)
    for (ii in 1:nrep) {
        bboxlist[[ii]] = switch(scheme,
                                random = .random_MARsampling(lonrange, latrange, sidesize),
                                inwards = NULL,
        )
    }
    bboxlist = unlist(bboxlist, recursive = FALSE)
    # calculate area and genetic diversity in each bounding boxes
    outlist = lapply(bboxlist, mar:::.mutdiv.gridded, gm = gm, gmarea = gmarea)
    outdf = do.call(dplyr::bind_rows, lapply(outlist, as.data.frame))
    # return bounding boxes if debug
    if(debug) {
        outdf$extent = unlist(lapply(bboxlist, paste0, collapse = ';'))
    }
    return(outdf)
}

# sampling methods
# old method but fix the minimum bug.
.random_MARsampling <- function(lonrange, latrange, sidesize) {
    bblist = lapply(sidesize, function(ss) {
        # xmin <- sample(1:(lonrange - ss + 1), 1); xmax <- xmin + ss - 1
        # ymin <- sample(1:(latrange - ss + 1), 1); ymax <- ymin + ss - 1
        xmin <- sample(1:(lonrange - ss), 1); xmax <- xmin + ss
        ymin <- sample(1:(latrange - ss), 1); ymax <- ymin + ss
        return(c(xmin, xmax, ymin, ymax))
    })
    return(bblist)
}

