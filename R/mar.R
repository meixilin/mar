require(dplyr)
################################################################################
# classic SFS
predSFS <- function(M, n = 1000) {
    # abundance = c 1/q
    # log(snps) = c
    # q=M/max(M,na.rm=TRUE) * 0.5 # because they are oriented to minor
    # q= round(q *100)/ 100
    # countsq=table(q)
    # lm(log(countsq) ~ log(fn(names(countsq))))
    samp <- sample(
        x = seq(1, 49, 0.01) / 100, size = length(M),
        prob = 1 / (seq(1, 49, 0.01) / 100), replace = TRUE
    )
    # samp=samp * max(M,na.rm=TRUE) #* 2 # because I used 0.49 as maximum frequency
    samp <- ceiling(samp * n) # to make this a
    return(fn(samp))
}
# sum of squares
pseudoR2 <- function(y, x) {
    mse <- sqrt(mean((y - x)^2)) # / sum(diff(range(x)))
    1 - mse / var(y)
}
# likelihood classic SFS

likelihoodSFSfold <- function(M) {
    # fit models
    m.ln <- fitsad(M, "lnorm")
    m.ls <- fitsad(M, "ls")
    m.sfs <- predSFS(M)
    # Predict expectation
    Mexpect_sfs <- m.sfs %>% ceiling()
    Mexpect_ln <- sads::radpred(m.ln)$abund %>% ceiling()
    Mexpect_ls <- sads::radpred(m.ls)$abund %>% ceiling()

    # Create tables
    tomerge_obs <- data.frame(possiblevalues = names(table(M)), Ofreq = fn(table(M)))
    tomerge_sfs <- data.frame(possiblevalues = names(table(Mexpect_sfs)), SFSfreq = fn(table(Mexpect_sfs)))
    tomerge_ln <- data.frame(possiblevalues = names(table(Mexpect_ln)), LNfreq = fn(table(Mexpect_ln)))
    tomerge_ls <- data.frame(possiblevalues = names(table(Mexpect_ls)), LSfreq = fn(table(Mexpect_ls)))

    # Merge
    tosfs <- merge(tomerge_obs, tomerge_sfs, by = "possiblevalues")
    tosfs <- merge(tosfs, tomerge_ln, by = "possiblevalues")
    tosfs <- merge(tosfs, tomerge_ls, by = "possiblevalues")

    # Calculate likelihood
    # if obs_i is the number of sites that have i derived alleles (so it is an integer) and exp_i is the proportion of sites that have i derived alleles under the model,
    # then the log likelihood is: # constant + \sum_i obs_i * log(exp_i)
    liksfs <- sum(tosfs$Ofreq * log(tosfs$SFSfreq / sum(tosfs$SFSfreq)))
    likln <- sum(tosfs$Ofreq * log(tosfs$LNfreq / sum(tosfs$LNfreq)))
    likls <- sum(tosfs$Ofreq * log(tosfs$LSfreq / sum(tosfs$LSfreq)))

    # AIC 2k - 2logL, but sum because the above is negative log likelihood, right?
    AICsfs <- 2 * 0 + 2 * sum(tosfs$Ofreq * log(tosfs$SFSfreq)) # 0 degrees of freedom
    AICln <- 2 * 2 + 2 * sum(tosfs$Ofreq * log(tosfs$LNfreq)) # 2 degrees of freedom
    AICls <- 2 * 1 + 2 * sum(tosfs$Ofreq * log(tosfs$LSfreq)) # 1 degrees of freedom

    return(list("sfs" = liksfs, "ln" = likln, "ls" = likls))
}
################################################################################
# Fast computation of SFS
# sourceCpp("R/genemaps.cpp")
# Create function to automatically compute SFS from subset of accessions
sfs <- function(bed, accessions = 1:nrow(bed), snps = 1:ncol(bed), folded = TRUE, nonzero = TRUE) {
    # tmp<-sfsC(bed,accessions,snps)
    # tmp<-apply(bed[accessions,snps],2,sum)
    tmp <- sfsmatC(bed[accessions, snps])
    # remove zeroes
    if (nonzero) tmp <- tmp[tmp != 0]
    # fold
    if (folded) tmp <- ifelse(tmp < length(accessions), tmp, length(accessions) * 2 - tmp)
    return(tmp)
}
sample_and_sfs <- function() {}

## OBJECT AND MAR BUILD ##############################################################################


# Spatial sampling
create_gene_maps <- function(coord, geno, bystep = 1, names = NULL) {
    require(raster)
    # stop condition
    stopifnot(dim(coord)[1] == dim(geno)[1])
    # create edges of the raster
    xmn <- range(coord[, 1])[1] - diff(range(coord[, 1]))*0.01
    xmx <- range(coord[, 1])[2] + diff(range(coord[, 1]))*0.01
    ymn <- range(coord[, 2])[1] - diff(range(coord[, 2]))*0.01
    ymx <- range(coord[, 2])[2] + diff(range(coord[, 2]))*0.01
    # create raster
    r <- raster(resolution = bystep, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)
    ## project coordinates in raster to count samples per grid point
    r_samples <- rasterize(coord, r, fun = "count")
    ## raster per mutation
    r_mutmaps <- lapply(
        # iterate through SNPs
        1:ncol(geno), function(j) {
            # get which SNPs have a hit in a genome and extract coords
            presentmut <- which(geno[, j] > 0)
            if (length(presentmut) > 0) {
                tmpcoords <- matrix(coord[presentmut, ], ncol = 2)
                tmpras <- rasterize(tmpcoords, r, fun = "count")
            } else {
                tmpras <- r
            }
            # rename if names provided
            if (!is.null(names)) names(tmpras) <- names[j]
            return(tmpras)
        }
    )
    r_mutmaps <- stack(r_mutmaps)
    # finish
    return(
        list(r_samples, r_mutmaps)
    )
}
get_total_area <- function(coord) {
    # create edges of the raster
    xmn <- range(coord[, 1])[1]
    xmx <- range(coord[, 1])[2]
    ymn <- range(coord[, 2])[1]
    ymx <- range(coord[, 2])[2]
    # create raster of 1 cell using the extent
    r <- raster(xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, nrows = 1, ncols = 1)
    # get the raster area
    values(raster::area(r))
}

# areapolygonofraster<-function(myraster){
#   tmp<-xyFromCell(myraster, which(values(myraster) >0))
#   tmp.mcp <- adehabitatHR::mcp(SpatialPoints(tmp,proj4string = myraster@crs),percent = 100,unout="km2")
#   tmp.rast<-hr.rast(tmp.mcp,myraster)
#   plot()
#   tmp.mcp$area
#   mcp.are
# }

mutVa <- function(raster_samples, raster_mutmaps, rest_mutmaps, betas) {
    stopifnot(nlayers(raster_mutmaps) == length(betas))
    require(raster)
    # Get the number of samples
    N <- sum(fn(values(raster_samples)), na.rm = T)
    # freqs
    P <- fn(apply(values(raster_mutmaps), 2, function(cells) {
        sum(cells > 0, na.rm = T)
    }))
    P <- P / N
    # compute Va
    Va <- sum(betas^2 * P * (1 - P))
    Suma2 <- sum(betas[P > 0]^2)
    # at least one cell has to hopave a presence of a mutation, and sum over
    M_ <- fn(apply(values(raster_mutmaps), 2, function(cells) any(cells > 0)))
    M_[is.na(M_)] <- 0
    M <- sum(M_)
    # find endemisms
    E_ <- apply(values(rest_mutmaps), 2, function(cells) any(cells > 0))
    E_[is.na(E_)] <- 0
    table(M_, E_)
    E <- sum(M_ & !E_)
    # Sum samples across cells
    N <- sum(values(raster_samples), na.rm = TRUE)
    # Get the number of SNPs for the sample
    L <- dim(raster_mutmaps)[3]
    # compute diversity, Theta Waterson
    if (N > 0 & M > 0) {
        theta <- M / (Hn(N) * L)
        thetapi <- sum(2 * P * (1 - P), na.rm = T) / L
    } else {
        theta <- 0
        thetapi <- 0
    }
    # area taking into account only grid cells with data
    # asub= sum(raster_samples[] > 0, na.rm = T) * (res(raster_samples)[1]*res(raster_samples)[2])
    asub <- areaofraster(raster_samples)

    # area based on simple square
    a <- dim(raster_samples)[1] * res(raster_samples)[1] * dim(raster_samples)[2] * res(raster_samples)[2]
    # return
    return(data.frame(
        thetaw = theta, pi = thetapi, M = M, E = E,
        N = N, a = a, asub = asub, Va = Va, Suma2 = Suma2
    ))
}


MARsampling_old <- function(genemaps, scheme = "random", samples = 10, centerfun = median) {
    require(raster)
    require(dplyr)
    ## Check conditions
    stopifnot(scheme %in% c("random", "inwards", "outwards", "southnorth"))
    stopifnot(is.function(centerfun))
    stopifnot(length(genemaps) == 2) # expecting density of samples and mutations
    ## End conditions

    # some general variables
    lonrange <- dim(genemaps[[1]])[1]
    latrange <- dim(genemaps[[1]])[2]
    ## random ############################################################################
    ## sampling boxes of sizes 1 ... min range with 1 stride but at random (too many combinations)
    if (scheme == "random") {
        minrange <- min(lonrange, latrange)
        # get minimum of lat or lon range to make square boxes for simplicity
        # iterate
        res <- lapply(1:minrange, function(sidesize) {
            lapply(1:samples, function(dum) {
                # create box of study
                # ML: small bug here when minrange reaches max. it sample(c(0,1)). Besides, the minimum sample size in this is 4 cells.
                xstart <- sample(1:(lonrange - sidesize), 1)
                ystart <- sample(1:(latrange - sidesize), 1)
                tmpextent <- extent(genemaps[[1]], xstart, xstart + sidesize, ystart, ystart + sidesize)
                # subset rasters of samples and mutation maps
                raster_samples <- crop(stack(genemaps[[1]]), tmpextent)
                raster_mutmaps <- crop(genemaps[[2]], tmpextent)
                # create a raster of reminder areas (for endemism study)
                rest_mutmaps <- genemaps[[2]]
                rest_mutmaps <- mask(rest_mutmaps, as(tmpextent, "SpatialPolygons"), inverse = TRUE)
                # compute diversity metrics and areas
                res = mutdiv(raster_samples, raster_mutmaps, rest_mutmaps)
                # DEBUG: add the extent stats
                res$extent = paste0(c(xstart, xstart + sidesize, ystart, ystart + sidesize), collapse = ';')
                return(res)
            }) %>% do.call(rbind, .)
        }) %>% do.call(rbind, .)
    }
    ## outwards ############################################################################
    ## sampling boxes of sizes 1 ... min range starting in the lower-left (low lat long) and moving up
    else if (scheme == "outwards") {
        # get most abundant location, or median
        tmp <- xyFromCell(genemaps[[1]], which(values(genemaps[[1]]) > 0))
        midcoord <- fn(apply(tmp, 2, centerfun)) %>% t()
        # id the central cell from which we will expand
        id.cell <- raster::extract(genemaps[[1]], SpatialPoints(midcoord), cellnumbers = TRUE)[1]
        startrowcol <- rowColFromCell(genemaps[[1]], id.cell)
        # iterate
        listres <- list()
        xstart <- startrowcol[1]
        xend <- startrowcol[1] + 1
        ystart <- startrowcol[2]
        yend <- startrowcol[2] + 1
        while (xstart > 0 | ystart > 0 | xend < lonrange | yend < yend) {
            # create box of study
            tmpextent <- extent(genemaps[[1]], xstart, xend, ystart, yend)
            # subset rasters of samples and mutation maps
            raster_samples <- crop(stack(genemaps[[1]]), tmpextent)
            raster_mutmaps <- crop(genemaps[[2]], tmpextent)
            # create a raster of reminder areas (for endemism study)
            rest_mutmaps <- mask(genemaps[[2]], as(tmpextent, "SpatialPolygons"), inverse = TRUE)
            # compute diversity metrics and areas
            listres <- c(listres, list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps)))
            # expand window
            xstart <- xstart - 1
            if (xstart < 0) xstart <- 0
            xend <- xend + 1
            if (xend > lonrange) xend <- lonrange
            ystart <- ystart - 1
            if (ystart < 0) ystart <- 0
            yend <- yend + 1
            if (yend > latrange) yend <- latrange
        }
        res <- listres %>% do.call(rbind, .)
    }
    ## inwards ############################################################################
    ## sampling boxes of sizes 1 ... min range starting in center of the distribution
    else if (scheme == "inwards") {
        # inwards is an identical process as outwards, but we flip the in and out objects
        # get most abundant location, or median
        tmp <- xyFromCell(genemaps[[1]], which(values(genemaps[[1]]) > 0))
        midcoord <- fn(apply(tmp, 2, centerfun)) %>% t()
        # id the central cell from which we will expand
        id.cell <- raster::extract(genemaps[[1]], SpatialPoints(midcoord), cellnumbers = TRUE)[1]
        startrowcol <- rowColFromCell(genemaps[[1]], id.cell)
        # iterate
        listres <- list()
        xstart <- startrowcol[1]
        xend <- startrowcol[1] + 1
        ystart <- startrowcol[2]
        yend <- startrowcol[2] + 1
        while (xstart > 0 | ystart > 0 | xend < lonrange | yend < yend) {
            # create box of study
            tmpextent <- extent(genemaps[[1]], xstart, xend, ystart, yend)
            # mask rasters of samples and mutation maps in the interior, has not been sampled yet as this is inwards
            raster_samples <- mask(genemaps[[1]], as(tmpextent, "SpatialPolygons"), inverse = TRUE)
            raster_mutmaps <- mask(genemaps[[2]], as(tmpextent, "SpatialPolygons"), inverse = TRUE)
            # create a raster of reminder areas (for endemism study). in this case inside distribution
            rest_mutmaps <- genemaps[[2]]
            rest_mutmaps <- mask(rest_mutmaps, as(tmpextent, "SpatialPolygons"), inverse = FALSE)
            # compute diversity metrics and areas
            listres <- c(listres, list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps)))
            # expand window
            xstart <- xstart - 1
            if (xstart < 0) xstart <- 0
            xend <- xend + 1
            if (xend > lonrange) xend <- lonrange
            ystart <- ystart - 1
            if (ystart < 0) ystart <- 0
            yend <- yend + 1
            if (yend > latrange) yend <- latrange
        }
        res <- listres %>% do.call(rbind, .)
    }
    ## southnorth ############################################################################
    ## sampling latitudinally
    else if (scheme == "southnorth") {
        # not implemented
        # # get most abundant location, or median
        # tmp<-xyFromCell(genemaps[[1]], which(values(genemaps[[1]]) >0))
        # midcoord<-fn(apply(tmp,2, centerfun)) %>% t()
        # # id the central cell from which we will expand
        # id.cell <- extract(genemaps[[1]],SpatialPoints(midcoord), cellnumbers=TRUE)[1]
        # startrowcol<-rowColFromCell(genemaps[[1]], id.cell)
        # # iterate
        # listres<-list()
        # xstart<-startrowcol[1]
        # xend<-startrowcol[1]+1
        # ystart<-startrowcol[2]
        # yend<-startrowcol[2]+1
        # while(xstart>0 | ystart>0 | xend<lonrange | yend<yend){
        #   # create box of study
        #   tmpextent<-extent(genemaps[[1]], xstart, xend, ystart, yend)
        #   # mask rasters of samples and mutation maps in the interior, has not been sampled yet as this is inwards
        #   raster_samples=mask(genemaps[[1]],as(tmpextent, 'SpatialPolygons'), inverse=TRUE)
        #   raster_mutmaps=mask(genemaps[[2]],as(tmpextent, 'SpatialPolygons'), inverse=TRUE)
        #   # create a raster of reminder areas (for endemism study). in this case inside distribution
        #   # rest_mutmaps=crop(genemaps[[2]],tmpextent)
        #   # rest_mutmaps = genemaps[[2]]
        #   rest_mutmaps= mask(rest_mutmaps,as(tmpextent, 'SpatialPolygons'), inverse=FALSE)
        #   # compute diversity metrics and areas
        #   listres<-c(listres,list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps)))
        #   # expand window
        #   xstart=xstart-1 ; if(xstart<0) xstart=0
        #   xend=xend+1     ; if(xend>lonrange) xend=lonrange
        #   ystart=ystart-1 ; if(ystart<0) ystart=0
        #   yend=yend+1     ; if(yend>latrange) yend=latrange
        # }
        # res<- listres %>% do.call(rbind,.)
    }
    # return any result
    return(res)
}


## THEORETICAL ##############################################################################
# Functions
MARextinction_ <- function(A, ax, marmod) {
    # A = original area. dymensions uninmporant unless rounding errors happen
    # ax = fraction of area extinct
    # marmod = a sar model
    Mnow <- sar_pred(marmod, A)
    Mfut <- sar_pred(marmod, round(A * (1 - ax)))
    (Mnow$Prediction - Mfut$Prediction) / Mnow$Prediction
}
MARextinction <- function(A, ax, marmod) {
    if (is.vector(ax)) {
        res <-
            sapply(ax, function(x) {
                MARextinction_(A, x, marmod)
            })
        data.frame(ax = ax, mx = res)
    } else {
        res <- MARextinction_(A, ax, marmod)
    }
    data.frame(ax = ax, mx = res)
}

## including Va ##############################################################################
MARextinction_Va <- function(bedfile, mapfile, famfile, latlonfile, effectfile, xfrac, xmod = "sn") {
    stopifnot(xmod %in% c("sn", "rand"))
    betas <- data.table::fread(effectfile)$beta
    fam <- data.table::fread(famfile)
    map <- data.table::fread(mapfile)
    bed <- readbed(bedfile,
        N = nrow(fam), p = nrow(map),
        myrows = 1:nrow(fam), mycols = 1:nrow(map)
    )
    ## Latitude and Longitude information
    latlon <- utils::read.csv(latlonfile)
    stopifnot(all(c("id", "latitude", "longitude") %in% colnames(latlon)))
    coords <- cbind(latlon$longitude, latlon$latitude)
    coords <- coords[!is.na(latlon$latitude) & !is.na(latlon$longitude), ]
    genomes <- bed[
        # remove those that do not have long/lat, previously parsed to match fam
        !is.na(latlon$latitude) & !is.na(latlon$longitude),
    ]
    # create gene maps
    genemaps <- create_gene_maps(coords, genomes, by = 1, names = map$V2)
    raster_samples <- genemaps[[1]]
    raster_mutmaps <- genemaps[[2]]
    rest_mutmaps <- genemaps[[2]]
    if (dim(raster_mutmaps)[[3]] != length(betas)) {
        names(betas) <- paste0("X", map$V2)
        betas <- fn(betas[names(raster_mutmaps)])
    }

    # extinction
    # MARextinction_sn<-function(genemaps,xfrac=0.01,centerfun=median, debug=FALSE){
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
        # list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps)))
        list(mutVa(raster_samples, raster_mutmaps, rest_mutmaps, betas))
    )
    while (A > 1) {
        # extinct some grids
        if (xmod == "sn") {
            toextinct <- sample(gridpresent, xstep,
                replace = TRUE,
                prob = rank(-northdist[gridpresent])^10
            ) # high prob with exponential decay
        } else if (xmod == "rand") {
            toextinct <- sample(gridpresent, xstep, replace = TRUE) # random probability
        }
        values(raster_samples)[toextinct] <- NA
        values(raster_mutmaps)[toextinct, ] <- NA
        values(rest_mutmaps)[toextinct, ] <- NA
        # calculate diversity
        listres <- c(
            listres,
            # list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps)))
            list(mutVa(raster_samples, raster_mutmaps, rest_mutmaps, betas))
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
            mx = 1 - (M / max(M, na.rm = T)),
            Vax = 1 - (Va / max(Va, na.rm = T))
        )
    return(res)
}

getextinction <- function(ax = 0.38, z = 0.3) {
    1 - (1 - ax)^z
}
