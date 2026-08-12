# SFS operations
.foldsfs <- function(vect) {
    flen <- floor(length(vect) / 2)
    fvect <- (vect + rev(vect))[1:flen]
    if (length(vect) %% 2 == 1) {
        fvect <- c(fvect, vect[flen + 1])
    }
    return(fvect)
}

# create a class called sfs, and handle folding and zeros
.new_sfs <- function(vect, folded, nozero) {
    stopifnot(class(vect) %in% c("numeric", "integer"))
    # not allowing NAs or less than zero values
    stopifnot(all(vect >= 0, !is.na(vect)))
    if (folded) {
        vect <- .foldsfs(vect)
    }
    if (nozero) {
        # remove the 0 and N*ploidy (all homref or all homalt sites)
        vect <- vect[-1] # always fold first
        names(vect) <- 1:length(vect)
    } else {
        names(vect) <- 0:(length(vect) - 1)
    }
    class(vect) <- c(class(vect), "sfs")
    attr(vect, "folded") <- folded
    attr(vect, "nozero") <- nozero
    return(vect)
}

#' Calculate the observed site frequency spectrum
#'
#' Tabulates the allele counts stored in a [genomaps()] object into a site
#' frequency spectrum (SFS): the number of sites carrying `k` alternative
#' alleles, for `k` from 0 to `N * ploidy`.
#'
#' @param gm A [genomaps()] object created by [genomaps()].
#' @param folded Logical, whether to fold the spectrum, collapsing the `k` and
#'   `N * ploidy - k` bins for sites whose ancestral allele is unknown.
#'   Default FALSE.
#' @param nozero Logical, whether to drop the `k = 0` bin (sites with no
#'   alternative allele in the sample). Default TRUE.
#'
#' @return An object of class `sfs`: a numeric vector of site counts, named by
#'   the allele count `k`, carrying the `folded` and `nozero` attributes.
#' @seealso [expsfs()] for the neutral expectation and [plot.sfs()] for plotting.
#' @export
#'
#' @examples
#' # unfolded spectrum of the 1001 genomes example
#' obssfs <- sfs(gm1001g)
#' head(obssfs)
#' plot(obssfs)

#' # folded spectrum, for data without a known ancestral allele
#' head(sfs(gm1001g, folded = TRUE))
# set folded = FALSE so the MAR sampling theory could work
sfs <- function(gm, folded = FALSE, nozero = TRUE) {
    AC <- gm$geno$allele_count
    N <- length(gm$maps$sample.id)
    ploidy <- gm$geno$ploidy

    xN <- N * ploidy
    if (any(AC > xN)) {
        warning(paste0(sum(AC > xN), " SNPs had allele counts exceeding N*ploidy"))
    }
    vect <- sapply(0:xN, function(x) sum(AC == x))
    # add class and handle folding
    vect <- .new_sfs(vect, folded, nozero)
    return(vect)
}

#' Generate the expected site frequency spectrum under neutrality
#'
#' Builds the standard neutral expectation `theta / k` for a [genomaps()]
#' object, scaling `theta` by Watterson's estimator so that the expected and
#' observed spectra contain the same number of segregating sites and are
#' directly comparable.
#'
#' @inheritParams sfs
#'
#' @return An object of class `sfs`, in the same format as the output of [sfs()],
#'   holding the expected rather than the observed number of sites per bin.
#' @seealso [sfs()] for the observed spectrum and [plot.sfs()] for plotting.
#' @export
#'
#' @examples
#' expsfs(gm1001g)[1:5]
#'
#' # compare the observed spectrum against the neutral expectation
#' plot(sfs(gm1001g), expected = expsfs(gm1001g), log = "x")
expsfs <- function(gm, folded = FALSE, nozero = TRUE) {
    N <- length(gm$maps$sample.id)
    ploidy <- gm$geno$ploidy
    AC <- gm$geno$allele_count
    xN <- N * ploidy
    # get segregating sites
    M_seg <- sum(AC > 0 & AC < xN)
    theta <- M_seg / .Hn(xN) # scale theta
    expsfs <- c(0, theta / (1:xN)) # need to add 0 as xN is the same as zero when folded
    expsfs <- .new_sfs(expsfs, folded, nozero)
    return(expsfs)
}
