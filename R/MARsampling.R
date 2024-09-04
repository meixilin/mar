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
    latrange <- dim(gm$samplemap)[1] # y is Row is Lat. Selected by r1, r2.
    lonrange <- dim(gm$samplemap)[2] # x is Col is Lon. Selected by c1, c2.
    minrange <- min(latrange, lonrange) # the maximum size box can become
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
                                southnorth = .pole_MARsampling(latrange, lonrange, sidesize, from = 'S'),
                                northsouth = .pole_MARsampling(latrange, lonrange, sidesize, from = 'N'),
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
.random_MARsampling <- function(latrange, lonrange, sidesize) {
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

