# combine area and genetic diversity calculation in a gridded bounding box
# bbox should be c(r1, r2, c1, c2).
mutdiv.gridded <- function(gm, gmarea, bbox, revbbox = FALSE) {
    stopifnot(length(bbox) == 4)
    r1 = bbox[1]; r2 = bbox[2]; c1 = bbox[3]; c2 = bbox[4]
    nrow = r2 - r1 + 1; ncol = c2 - c1 + 1
    resrow = res(gm$samplemap)[1]; rescol = res(gm$samplemap)[2]
    # calculate area by size of bounding box
    Asq <- areaofsquare(nrow, ncol, resrow, rescol)
    # if reverse bounding box
    if(revbbox) {
        Asq <- areaofsquare(dim(gm$samplemap)[1], dim(gm$samplemap)[2], resrow, rescol) - Asq
    }
    # locate cellids from bbox
    cellids <- rowcol_cellid(gm, bbox, revbbox = revbbox)
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
