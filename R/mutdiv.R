# combine area and genetic diversity calculation in a gridded bounding box
# bbox should be c(xmin, xmax, ymin, ymax).
.mutdiv.gridded <- function(gm, gmarea, bbox) {
    stopifnot(length(bbox) == 4)
    xmin = bbox[1]; xmax = bbox[2]; ymin = bbox[3]; ymax = bbox[4]
    nrow = xmax - xmin + 1; ncol = ymax - ymin + 1
    resrow = res(gm$samplemap)[1]; rescol = res(gm$samplemap)[2]
    # calculate area by size of bounding box
    Asq <- areaofsquare(nrow, ncol, resrow, rescol)
    # locate cellids from bbox
    cellids <- rowcol_cellid(gm, bbox)
    # if no cells
    if (length(cellids) == 0) {
        out = list(Asq = Asq)
    } else {
        # calculate area by filled raster
        A <- sum(gmarea[cellids])
        # get samples
        sampleids <- cellid_sample(gm, cellids)
        # calculate genetic diversity
        out <- append(.calc_theta(gm, sampleids),
                      list(A = A,
                           Asq = Asq)
        )
    }
    return(out)
}

# genetic diversity estimator (use `gm$genotype` matrix)
# TODO: maybe needs to consider diploids?
.calc_theta <- function(gm, sampleid = NULL) {
    # get length of genome
    L <- attr(gm, "genolen")
    # subset ids
    if (is.null(sampleid)) {
        sampleid = 1:(dim(gm$genotype)[2])
    }
    # ingeno <- gm$genotype[sampleid, , drop = FALSE]
    # outgeno <- gm$genotype[-sampleid, , drop = FALSE]

    # number of samples
    N <- length(sampleid) # dim(ingeno)[1]

    AC <- matrixStats::rowSums2(gm$genotype, cols = sampleid)
    oAC <- matrixStats::rowSums2(gm$genotype, cols = -sampleid)

    # segregating sites
    M <- sum(AC > 0)
    # TODO: ploid. allele frequency
    P <- AC/N
    # compute diversity, Theta Waterson and Theta Pi (pairwise)
    if (N > 1 & M > 0) {
        thetaw <- M / (Hn(N-1) * L)
        thetapi <- (N/(N-1)) * sum(2 * P * (1 - P), na.rm = T) / L
    } else {
        thetaw <- 0
        thetapi <- 0
    }
    # endemic segregating sites
    E <- sum(AC > 0 & oAC == 0)

    # return a list
    out <- list(N = N,
                M = M,
                E = E,
                thetaw = thetaw,
                thetapi = thetapi)
    return(out)
}

# genetic diversity estimator
mutdiv_old <- function(raster_samples, raster_mutmaps, rest_mutmaps) {
    require(raster)
    # Get the number of samples
    N <- sum(fn(values(raster_samples)), na.rm = T)
    # freqs
    P <- fn(apply(values(raster_mutmaps), 2, function(cells) sum(cells, na.rm = T))) / N
    # at least one cell has to have a presence of a mutation, and sum over
    M_ <- fn(apply(values(raster_mutmaps), 2, function(cells) any(cells > 0)))
    M_[is.na(M_)] <- 0
    M <- sum(M_)
    # find endemisms
    E_ <- apply(values(rest_mutmaps), 2, function(cells) any(cells > 0))
    E_[is.na(E_)] <- 0
    table(M_, E_)
    E <- sum(M_ & !E_)
    # Get the number of SNPs for the sample
    L <- dim(raster_mutmaps)[3]
    # compute diversity, Theta Waterson and Theta Pi (pairwise)
    if (N > 1 & M > 0) {
        theta <- M / (Hn(N-1) * L)
        thetapi <- (N/(N-1)) * sum(2 * P * (1 - P), na.rm = T) / L
    } else {
        theta <- 0
        thetapi <- 0
    }
    # area taking into account only grid cells with data
    # asub= sum(raster_samples[] > 0, na.rm = T) * (res(raster_samples)[1]*res(raster_samples)[2])
    asub <- areaofraster(raster_samples, cached = FALSE, na.rm = TRUE)

    # area based on simple square
    a <- dim(raster_samples)[1] * res(raster_samples)[1] * dim(raster_samples)[2] * res(raster_samples)[2]
    # return
    return(data.frame(thetaw = theta, pi = thetapi, M = M, E = E, N = N, a = a, asub = asub))
}
