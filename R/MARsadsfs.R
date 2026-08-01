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

#' Calculate Site Frequency Spectrum
#'
#' @param AC Vector of allele counts
#' @param N Number of samples
#' @param ploidy Ploidy level of the organism
#' @param folded Logical, whether to fold the spectrum. Default TRUE
#' @param nozero Logical, whether to remove zero counts. Default TRUE
#'
#' @return An object of class "sfs" containing the site frequency spectrum
#' @export
#'
#' @examples
#' # Calculate SFS from allele counts
#' allele_counts <- c(1, 1, 0, 2, 0, 1, 1, 0, 0, 2, 30)
#' sfs_result <- sfs(allele_counts, N = 50, ploidy = 2)
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

#' Generate Expected Site Frequency Spectrum
#'
#' @param lenAC Length of allele count vector. In other words, total number of SNPs surveyed.
#' @param N Number of samples
#' @param ploidy Ploidy level of the organism
#' @param folded Logical, whether to fold the spectrum. Default TRUE
#' @param nozero Logical, whether to remove zero counts. Default TRUE
#'
#' @return An object of class "sfs" containing the expected site frequency spectrum
#' @export
#'
#' @examples
#' # Generate expected SFS
#' exp_sfs <- expsfs(lenAC = 1000, N = 100, ploidy = 2)
expsfs <- function(gm, folded = FALSE, nozero = TRUE) {
    N <- length(gm$maps$sample.id)
    ploidy <- gm$geno$ploidy
    xN <- N * ploidy
    lenAC <- nrow(gm$geno$genotype)
    theta <- lenAC / .Hn(xN) # scale theta
    expsfs <- c(0, theta / (1:xN)) # need to add 0 as xN is the same as zero when folded
    expsfs <- .new_sfs(expsfs, folded, nozero)
    return(expsfs)
}
