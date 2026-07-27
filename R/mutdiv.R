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
#' div <- mutdiv.gridded(gm1001g, gmarea, bbox = c(1, 10, 2, 8))
mutdiv.gridded <- function(gm, gmarea, sm, bbox, revbbox = FALSE) {
    stopifnot('Bounding box should be length of 4' = length(bbox) == 4)
    # calculate area by bounding box using gmarea (so units are interpretable at km)
    Asq <- sum(gmarea[bbox[1]:bbox[2], bbox[3]:bbox[4]])
    # if reverse bounding box
    if (revbbox) {
        Asq <- terra::global(gmarea, fun = "sum", na.rm = TRUE)[1, 1] - Asq
    }
    # locate cellids from bbox
    cellids <- .rowcol_cellid(sm, bbox, revbbox = revbbox)
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
#' @export
#'
#' @examples
#' # Calculate genetic diversity for a specific set of cells
#' gmarea <- mar:::.areaofraster(gm1001g$maps$samplemap)
#' cell_ids <- c(613, 726, 727)
#' div <- mutdiv.cells(gm1001g, gmarea, cellids = cell_ids)
mutdiv.cells <- function(gm, gmarea, cellids) {
    # Not calculating Asq as it is very similar to A in MARextinction settings
    out <- .mutdiv.cellids(gm, gmarea, cellids, NULL)
    return(out)
}

.mutdiv.cellids <- function(gm, gmarea, cellids, Asq) {
    # if no cells
    if (length(cellids) == 0) {
        out <- list(N = NA, M = NA, E = NA, thetaw = NA, thetapi = NA, A = NA)
    } else {
        # calculate area by filled raster
        A <- sum(gmarea[cellids])
        # get samples
        sampleids <- .cellid_sample(gm$maps, cellids)
        # calculate genetic diversity
        out <- append(
            .calc_theta(gm, sampleids),
            list(A = A))
    }
    # Asq will be appended when called from mutdiv.gridded
    if (!is.null(Asq)) {
        out <- append(out, list(Asq = Asq))
    }
    return(out)
}

# genetic diversity estimator (use `gm$genotype` matrix)
# ploidy does not matter here. although > diploid is not well-defined.
# TODO: allow L calculations
.calc_theta <- function(gm, sampleid = NULL) {
    ploidy <- gm$geno$ploidy
    # subset ids
    if (is.null(sampleid)) {
        sampleid <- 1:(dim(gm$geno$genotype)[2])
    }

    # number of samples (need to scale by ploidy)
    N <- length(sampleid) # dim(ingeno)[1]

    AC <- matrixStats::rowSums2(gm$geno$genotype, cols = sampleid, na.rm = TRUE)
    oAC <- gm$geno$allele_count - AC

    # number of called alleles (allows missing data now)
    xN <- (N - matrixStats::rowCounts(gm$geno$genotype, cols = sampleid, value = NA_integer_)) * ploidy

    # number of mutations (can be all alternative mutations)
    M <- sum(AC > 0)
    # segregating sites (polymorphic in the sample)
    M_seg <- sum(AC > 0 & AC < xN)
    # compute diversity, Theta Waterson and Theta Pi (pairwise)
    if (any(xN > 1) & M > 0) {
        # total pairwise difference / total pairwise comparison
        thetapi <- sum(2 * AC * (xN - AC)) / sum(xN * (xN - 1))
        # Segregating sites / sum of all possible harmonic numbers of xN (subset to never evaluate .Hn(-1))
        thetaw <- M_seg / sum(.Hn(xN[xN > 0] - 1))
    } else {
        thetaw <- 0
        thetapi <- 0
    }
    # endemic segregating sites
    E <- sum(AC > 0 & oAC == 0)

    # return a list
    out <- list(
        N = N,
        M = M,
        E = E,
        thetaw = thetaw,
        thetapi = thetapi
    )
    return(out)
}
