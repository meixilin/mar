#' Calculate Genetic Diversity in a Gridded Bounding Box
#'
#' @param gm A genomaps object containing genetic data and geographic informations, created by [genomaps()] function.
#' @param gmarea Raster file that contains area size of each cell
#' @param bbox Numeric vector of length 4 specifying the bounding box coordinates c(r1, r2, c1, c2)
#' @param revbbox Logical, whether to reverse/invert the bounding box selection. Default FALSE
#'
#' @return A list containing:
#'   \item{N}{Number of samples}
#'   \item{M}{Number of segregating sites}
#'   \item{E}{Number of endemic segregating sites}
#'   \item{thetaw}{Watterson's theta estimate}
#'   \item{thetapi}{Pi (pairwise) diversity estimate}
#'   \item{A}{Total area with data}
#'   \item{Asq}{Area of the bounding box square}
#' @export
#'
#' @examples
#' # Calculate mutation diversity in a 10x10 grid region
#' gmarea <- mar:::.areaofraster(gm1001g$maps$samplemap)
#' div <- mutdiv.gridded(gm1001g, gmarea, bbox=c(1,10,2,8))
mutdiv.gridded <- function(gm, gmarea, bbox, revbbox = FALSE) {
    stopifnot(length(bbox) == 4)
    r1 = bbox[1]; r2 = bbox[2]; c1 = bbox[3]; c2 = bbox[4]
    nrow = r2 - r1 + 1; ncol = c2 - c1 + 1
    resrow = raster::res(gm$maps$samplemap)[1]; rescol = raster::res(gm$maps$samplemap)[2]
    # calculate area by size of bounding box
    Asq <- .areaofsquare(nrow, ncol, resrow, rescol)
    # if reverse bounding box
    if(revbbox) {
        Asq <- .areaofsquare(dim(gm$maps$samplemap)[1], dim(gm$maps$samplemap)[2], resrow, rescol) - Asq
    }
    # locate cellids from bbox
    cellids <- .rowcol_cellid(gm$maps, bbox, revbbox = revbbox)
    out <- .mutdiv.cellids(gm, gmarea, cellids, Asq)
    return(out)
}

#' Calculate Genetic Diversity for Specified Cells
#'
#' Calculates genetic diversity metrics for a specific set of cells. This function is particularly
#' useful for extinction simulations.
#'
#' @param gm A genomaps object containing genetic data and geographic informations, created by [genomaps()] function.
#' @param gmarea Raster file that contains area size of each cell
#' @param cellids Vector of cell IDs to analyze
#'
#' @return A list containing:
#'   \item{N}{Number of samples}
#'   \item{M}{Number of segregating sites}
#'   \item{E}{Number of endemic segregating sites}
#'   \item{thetaw}{Watterson's theta estimate}
#'   \item{thetapi}{Pi (pairwise) diversity estimate}
#'   \item{A}{Total area with data}
#'   \item{Asq}{Area of the cells being analyzed}
#' @export
#'
#' @examples
#' # Calculate genetic diversity for a specific set of cells
#' gmarea <- mar:::.areaofraster(gm1001g$maps$samplemap)
#' cell_ids <- c(613,726,727)
#' div <- mutdiv.cells(gm1001g, gmarea, cellids=cell_ids)
mutdiv.cells <- function(gm, gmarea, cellids) {
    resrow = raster::res(gm$maps$samplemap)[1]; rescol = raster::res(gm$maps$samplemap)[2]
    # calculate area by size of bounding box
    Asq <- .areaofsquare(length(cellids), ncol = 1, resrow, rescol)
    out <- .mutdiv.cellids(gm, gmarea, cellids, Asq)
    return(out)
}

.mutdiv.cellids <- function(gm, gmarea, cellids, Asq) {
    # if no cells
    if (length(cellids) == 0) {
        out <- list(N = NA, M = NA, E = NA, thetaw = NA, thetapi = NA, A = NA, Asq = Asq)
    } else {
        # calculate area by filled raster
        A <- sum(gmarea[cellids])
        # get samples
        sampleids <- .cellid_sample(gm$maps, cellids)
        # calculate genetic diversity
        out <- append(.calc_theta(gm, sampleids),
                      list(A = A,
                           Asq = Asq)
        )
    }
    return(out)
}

# genetic diversity estimator (use `gm$genotype` matrix)
# ploidy does not matter here. although > diploid is not well-defined.
# TODO: allow L calculations
.calc_theta <- function(gm, sampleid = NULL) {
    # get length of genome
    L <- attr(gm, "genolen")
    ploidy <- .get_genodata(gm$geno, "ploidy")
    # subset ids
    if (is.null(sampleid)) {
        sampleid = 1:(dim(gm$geno$genotype)[2])
    }
    # ingeno <- gm$genotype[sampleid, , drop = FALSE]
    # outgeno <- gm$genotype[-sampleid, , drop = FALSE]

    # number of samples (need to scale by ploidy)
    N <- length(sampleid) # dim(ingeno)[1]
    xN <- N * ploidy

    AC <- matrixStats::rowSums2(gm$geno$genotype, cols = sampleid)
    oAC <- matrixStats::rowSums2(gm$geno$genotype, cols = -sampleid)

    # segregating sites
    M <- sum(AC > 0)
    # allele frequency
    P <- AC/xN
    # compute diversity, Theta Waterson and Theta Pi (pairwise)
    if (xN > 1 & M > 0) {
        thetaw <- M / (.Hn(xN-1) * L)
        thetapi <- (xN/(xN-1)) * sum(2 * P * (1 - P), na.rm = T) / L
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
