
## SIM EXTINCT ##############################################################################

MARextinction_radial <- function(genemaps, xfrac = 0.01, centerfun = median, debug = FALSE) {
    require(raster)
    require(dplyr)
    raster_samples <- genemaps[[1]]
    raster_mutmaps <- genemaps[[2]]
    rest_mutmaps <- raster_mutmaps
    # get most abundant location, or median
    tmp <- xyFromCell(genemaps[[1]], which(values(genemaps[[1]]) > 0))
    startcoordextinction <- fn(apply(tmp, 2, centerfun)) %>%
        t() %>%
        data.frame() %>%
        rename(x = X1, y = X2)
    # extinction cell from which we will expand extinction
    id.cell <- raster::extract(genemaps[[1]], SpatialPoints(startcoordextinction), cellnumbers = TRUE)[1]
    startrowcol <- rowColFromCell(genemaps[[1]], id.cell)
    # get the  present locations
    gridpresent <- which(apply(values(raster_mutmaps), 1, function(x) any(!is.na(x))) == TRUE)
    A <- length(gridpresent)
    Astart <- A
    xstep <- ceiling(xfrac * A)
    # get the latlong of every cell in the grid
    locs <- raster::as.data.frame(raster_samples, xy = TRUE)
    # get distance to the extinction point
    alldist <- as.matrix(dist(rbind(startcoordextinction, locs[, 1:2]), method = "euclidean"))[1, -1] %>% fn()
    # iterate
    listres <- list()
    # calculate original diversity
    listres <- c(
        listres,
        list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps))
    )
    if (debug) print(listres)
    while (A > 1) { # change 0 for Astop if wanted to stop earlier
        if (debug) message("A ", A)
        # extinct some grids. get the top that are closest in distance
        # Modification: Change to <= instead of < so when xstep == 1, still doable
        # Date: Mon Mar 13 00:35:57 2023
        toextinct <- gridpresent[which(alldist[gridpresent] <= sort(alldist[gridpresent])[xstep])]
        # extinct those values
        values(raster_samples)[toextinct] <- NA
        values(raster_mutmaps)[toextinct, ] <- NA
        values(rest_mutmaps)[toextinct, ] <- NA
        # calculate diversity
        tmpdiv <- mutdiv(raster_samples, raster_mutmaps, rest_mutmaps)
        listres <- c(listres, list(tmpdiv))
        # recalculate area remaining
        gridpresent <- which(apply(
            values(raster_mutmaps), 1,
            function(x) {
                any(!is.na(x))
            }
        ))
        A <- A - xstep
        if (debug) message("asub ", tmpdiv$asub)
        # if(debug) plot(raster_samples>100) # ** for debug# for debug
    }
    ## end function
    res <- listres %>%
        do.call(rbind, .) %>%
        data.frame() %>%
        mutate(
            ax = 1 - (asub / max(asub, na.rm = T)),
            mx = 1 - (M / max(M, na.rm = T))
        )
    return(res)
}

MARextinction_sn <- function(genemaps, xfrac = 0.01, centerfun = median, debug = FALSE) {
    require(raster)
    require(dplyr)
    raster_samples <- genemaps[[1]]
    raster_mutmaps <- genemaps[[2]]
    rest_mutmaps <- raster_mutmaps
    # get the  present locations
    gridpresent <- which(apply(values(raster_mutmaps), 1, function(x) any(!is.na(x))) == TRUE)
    A <- length(gridpresent)
    Astart <- A
    xstep <- ceiling(xfrac * A)
    # get the latlong of every cell in the grid
    locs <- raster::as.data.frame(raster_samples, xy = TRUE)
    # distance to top latitude would lead to South to North prob extinction
    # fixed for situations where pole is south. get the location
    # no mater if negative or positive, that mathes the largest number
    northdist <- 90 - locs$y
    if (mean(northdist, na.rm = T) < 0) northdist <- northdist * (-1)
    # iterate object
    listres <- list()
    # calculate original diversity
    listres <- c(
        listres,
        list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps))
    )
    if (debug) {
        listext <- list()
    }
    while (A > 1) {
        # extinct some grids
        toextinct <- sample(gridpresent, xstep,
            replace = TRUE,
            # prob =normalize(northdist[gridpresent])^10)  # likely created a bug with normalize
            prob = rank(-northdist[gridpresent])^10
        ) # high prob with exponential decay
        if (debug) {
            listext <- c(listext, list(toextinct))
            print(toextinct)
        }
        values(raster_samples)[toextinct] <- NA
        values(raster_mutmaps)[toextinct, ] <- NA
        values(rest_mutmaps)[toextinct, ] <- NA
        # calculate diversity
        listres <- c(
            listres,
            list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps))
        )
        # plot(raster_samples>0, col='black') # ** for debug# for debug
        # recalculate area remaining
        gridpresent <- which(apply(values(raster_mutmaps), 1, function(x) any(!is.na(x))) == TRUE)
        A <- A - xstep
    }
    res <- listres %>%
        do.call(rbind, .) %>%
        data.frame() %>%
        mutate(
            ax = 1 - (asub / max(asub, na.rm = T)),
            mx = 1 - (M / max(M, na.rm = T))
        )
    if (debug) {
        return(list(res, listext))
    } else {
        return(res)
    }
}

MARextinction_in <- function(genemaps, xfrac = 0.01, centerfun = median, debug = FALSE) {
    require(raster)
    require(dplyr)
    raster_samples <- genemaps[[1]]
    raster_mutmaps <- genemaps[[2]]
    rest_mutmaps <- raster_mutmaps
    # get the  present locations
    gridpresent <- which(apply(values(raster_mutmaps), 1, function(x) any(!is.na(x))) == TRUE)
    A <- length(gridpresent)
    Astart <- A
    xstep <- ceiling(xfrac * A)
    # get the median coordinate of each Arabidopsis sample
    tmp <- raster::xyFromCell(raster_samples, which(values(raster_samples) > 0))
    midcoord <- fn(apply(tmp, 2, centerfun)) %>%
        t() %>%
        data.frame() %>%
        rename(x = X1, y = X2)
    # get the latlong of every cell in the grid
    locs <- raster::as.data.frame(raster_samples, xy = TRUE)
    # get distance center to each cell
    alldist <- as.matrix(dist(rbind(midcoord, locs[, 1:2]), method = "euclidean"))[1, -1]
    # iterate object
    listres <- list()
    # calculate original diversity
    listres <- c(
        listres,
        list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps))
    )
    while (A > 1) { # change 0 for Astop if want to stop earlier
        # extinct some grids
        toextinct <- sample(
            x = gridpresent, size = xstep, replace = TRUE,
            prob = normalize(alldist[gridpresent])^10
        ) # high prob with exponential decay
        values(raster_samples)[toextinct] <- NA
        values(raster_mutmaps)[toextinct, ] <- NA
        values(rest_mutmaps)[toextinct, ] <- NA
        # calculate diversity
        listres <- c(
            listres,
            list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps))
        )
        # for debug plot(raster_samples>0)
        # recalculate area remaining
        gridpresent <- which(apply(values(raster_mutmaps), 1, function(x) any(!is.na(x))) == TRUE)
        A <- A - xstep
    }
    res <- listres %>%
        do.call(rbind, .) %>%
        data.frame() %>%
        mutate(
            ax = 1 - (asub / max(asub, na.rm = T)),
            mx = 1 - (M / max(M, na.rm = T))
        )
    return(res)
}

MARextinction_random <- function(genemaps, xfrac = 0.01, centerfun = median, debug = FALSE) {
    require(raster)
    require(dplyr)
    raster_samples <- genemaps[[1]]
    raster_mutmaps <- genemaps[[2]]
    rest_mutmaps <- raster_mutmaps
    # get the  present locations
    gridpresent <- which(apply(values(raster_mutmaps), 1, function(x) any(!is.na(x))) == TRUE)
    A <- length(gridpresent)
    Astart <- A
    xstep <- ceiling(xfrac * A)
    # iterate
    listres <- list()
    # calculate original diversity
    listres <- c(
        listres,
        list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps))
    )
    while (A > 1) {
        # extinct some grids
        toextinct <- sample(gridpresent, xstep, replace = TRUE)
        values(raster_mutmaps)[toextinct, ] <- NA
        values(raster_samples)[toextinct] <- NA
        values(rest_mutmaps)[toextinct, ] <- NA
        # calculate diversity
        listres <- c(
            listres,
            list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps))
        )
        # recalculate area remaining
        gridpresent <- which(apply(values(raster_mutmaps), 1, function(x) any(!is.na(x))) == TRUE)
        A <- A - xstep
    }
    res <- listres %>%
        do.call(rbind, .) %>%
        data.frame() %>%
        mutate(
            ax = 1 - (asub / max(asub, na.rm = T)),
            mx = 1 - (M / max(M, na.rm = T))
        )
    return(res)
}


MARextinction_sim <- function(genemaps, scheme = "random",
                              samples = 10, xfrac = 0.01,
                              centerfun = median, debug = FALSE) {
    require(raster)
    require(dplyr)
    ## Check conditions
    stopifnot(scheme %in% c("random", "inwards", "outwards", "southnorth", "radial"))
    stopifnot(is.function(centerfun))
    stopifnot(length(genemaps) == 2) # expecting density of samples and mutations
    ## End conditions
    # # some general variables
    # lonrange <- dim(genemaps[[1]])[1]
    # latrange <- dim(genemaps[[1]])[2]
    ## random #############################################################################
    ## sampling boxes of sizes 1 ... min range with 1 stride but at random (too many combinations)
    if (scheme == "random") {
        res <- MARextinction_random(genemaps = genemaps, xfrac = xfrac, centerfun = centerfun)
    } # end condition
    ## outwards ############################################################################
    ## sampling from center to edges
    else if (scheme == "outwards") {
        stop("outwards extinction not implemented yet")
    } # end condition
    ## inwards ############################################################################
    ## sampling from center to edges
    else if (scheme == "inwards") {
        res <- MARextinction_in(genemaps = genemaps, xfrac = xfrac, centerfun = centerfun)
    } # end condition
    ## southnorth ############################################################################
    ## Sampling grids from south to north
    else if (scheme == "southnorth") {
        res <- MARextinction_sn(genemaps = genemaps, xfrac = xfrac, centerfun = centerfun, debug = debug)
    } # end condition
    ## radial ############################################################################
    else if (scheme == "radial") {
        res <- MARextinction_radial(genemaps = genemaps, xfrac = xfrac, centerfun = centerfun)
    } # end condition
    # end of funciton
    return(res)
}

