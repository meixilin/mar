# Predict the mutations-area relationship (MAR) from a site frequency spectrum (SFS)
#
# The theory: under the
# assumptions that individuals are randomly distributed in space and that
# mutations are unlinked, the expected number of segregating sites captured in a
# subsample of `n` chromosomes is a deterministic function of the full SFS.

# Variable lists:
# xN = total sequenced chromosomes (N * ploidy)
# L = total length of sequence (the number of trials)
# k = number of alleles in xN samples per site (SFS x-axis)
# p_K = The probability of sampling at least one mutation in $i$-th position in the $n$ samples given that this location has $k$ alleles in the $N$ total chromosomes

# ---- math helpers (ported verbatim from the SFS->MAR theory) --------------

# ratio of the choose functions C(xN-k,n)/C(xN,n) (hypergeometric "miss" prob).
# Done on the log scale: lchoose is vectorized over k and avoids the overflow
# and error accumulation of the running product. n = 0 gives 1, and k > xN - n
# gives 0 (lchoose returns -Inf), both of which are the intended limits.
.choose_Nkn <- function(xN, k, n) {
    result <- exp(lchoose(xN - k, n) - lchoose(xN, n))
    return(result)
}

# ratio of the choose functions C(k,n)/C(xN,n): the probability that every
# sampled chromosome carries the derived allele (the site looks fixed in the
# subsample). Zero when k < n, again via lchoose(k, n) = -Inf.
.choose_kn <- function(xN, k, n) {
    result <- exp(lchoose(k, n) - lchoose(xN, n))
    return(result)
}

# Sampling is modelled over chromosomes. This assumes a randomly breeding
# population: under Hardy-Weinberg an individual is ploidy independent gene
# copies, so drawing n = ploidy * j chromosomes matches drawing j individuals,
# and the scheme carries over to any ploidy. It does NOT hold when the copies
# within an individual are correlated (inbred lines, selfers): there the
# chromosome hypergeometric overstates the chance of catching a mutation and
# M / thetaw / thetapi are all overpredicted at small sample sizes.

.p_K <- function(xN, k, n, pk) {
    miss <- .choose_Nkn(xN, k, n)
    pk * (1 - miss)
}

# same, but the site must be SEGREGATING in the subsample (1 <= K <= n-1), so
# both fixed outcomes are excluded: no derived copy sampled, or only derived
# copies sampled.
.p_seg <- function(xN, k, n, pk) {
    miss <- .choose_Nkn(xN, k, n)
    allderiv <- .choose_kn(xN, k, n)
    pk * (1 - miss - allderiv)
}

# expected number of mutations E[S] captured at subsample size n: a site counts
# as soon as one derived copy lands in the subsample.
# kpdf is named 0:xN, so position i holds the allele count k = i - 1.
# k = 0 contributes nothing (miss = 1), it is kept in the sum for clarity.
.E_M <- function(xN, kpdf, n, L) {
    kk <- 0:xN
    sumpk <- sum(.p_K(xN = xN, k = kk, n = n, pk = kpdf[kk + 1]))
    out <- L * sumpk
    return(out)
}

# expected Watterson's theta E[thetaw] at subsample size n. Watterson's
# estimator is defined on segregating sites, so .p_seg is used rather than
# .p_K: counting sites that are fixed-derived in the subsample would leave
# thetaw non-flat even for a neutral (1/k) SFS, where it must equal theta.
.E_thetaw <- function(xN, kpdf, n) {
    kk <- 0:xN
    sumseg <- sum(.p_seg(xN = xN, k = kk, n = n, pk = kpdf[kk + 1]))
    out <- sumseg / .Hn(n - 1)
    return(out)
}

# expected nucleotide diversity E[thetapi], which does not depend on n:
#
#   E[pi_n] = sum_k kpdf_k * 2k(xN-k) / (xN(xN-1))   for every n >= 2
#
# pi is an average over PAIRS, and every pair drawn in the subsample is a
# uniformly random pair from the panel, so E[#discordant] = choose(n,2) *
# (panel discordant fraction) and the choose(n,2) cancels when dividing. The
# expected pairwise difference is therefore left untouched by subsampling, for
# any SFS -- no model assumption, just exchangeability. That makes the
# hypergeometric re-projection of the SFS unnecessary: this is O(xN) instead of
# O(n * xN) per n, and agrees with the projection to floating point (< 1e-16).
.E_thetapi <- function(xN, kpdf) {
    kk <- 0:xN
    out <- sum(kpdf[kk + 1] * 2 * kk * (xN - kk) / (xN * (xN - 1)))
    return(out)
}

# expected diversity at a single subsample size n. n is first input here for lapply
.mutdiv.theory <- function(n, xN, kpdf, L, ploidy) {
    M <- .E_M(xN = xN, kpdf = kpdf, n = n, L = L)
    # thetaw / thetapi are undefined at n = 1 (H_0 = 0 / division by n - 1)
    if (n < 2) {
        thetaw <- NA_real_
        thetapi <- NA_real_
    } else {
        thetaw <- .E_thetaw(xN = xN, kpdf = kpdf, n = n)
        thetapi <- .E_thetapi(xN = xN, kpdf = kpdf)
    }
    # N is the number of individuals, not the chromosomes anymore
    # (comparable to MARsampling's N column)
    out <- list(N = n / ploidy, M = M, thetaw = thetaw, thetapi = thetapi)
    return(out)
}

# ---- user-facing function -------------------------------------------------

#' Predict the mutations-area relationship (MAR) from a site frequency spectrum
#'
#' `MARtheory()` is the theoretical counterpart to [MARsampling()]. Instead of
#' repeatedly sub-sampling the map, it computes the *expected* accumulation of
#' mutations (and Watterson's theta / nucleotide diversity) with sample size
#' directly from the site frequency spectrum, under the assumptions that
#' individuals are randomly distributed in space and that mutations are unlinked.
#'
#' The returned data frame shares the `N`, `M`, `thetaw` and `thetapi` columns
#' of the [MARsampling()] / [MARextinction()] output, so it can be overlaid on
#' the empirical sampling points, or passed to [MARcalc()] with `Atype = "N"` to
#' fit the power-law exponent against sample size.
#'
#' @param gm A [genomaps()] object created by [genomaps()].
#' @param maxind Maximum number of individuals to predict up to. Default to all
#'   individuals available in genomaps()
#'
#' @return A data frame with class `martheory`, with a column `N` (number of
#'   individuals, from 1 to `maxind`) and columns `M`, `thetaw`, and `thetapi`.
#' @seealso [sfs()] for the spectrum the prediction is built from,
#'   [MARsampling()] for the empirical counterpart, and [plot.martheory()].
#' @export
#'
#' @examples
#' # Predicted MAR for the 1001 genomes example
#' expmardf <- MARtheory(gm1001g)
#' head(expmardf)
#'
#' # overlay the prediction on the empirical sampling points
#' mardf <- MARsampling(gm1001g, nrep = 2, xfrac = 0.1, myseed = 42)
#' plot(mardf, Atype = "N", theory = expmardf)
#'
#' # fit the predicted power-law exponent against sample size
#' MARcalc(expmardf, Mtype = "M", Atype = "N")
MARtheory <- function(gm, maxind = NULL) {
    stopifnot("gm must be a `genomaps` object" = inherits(gm, "genomaps"))

    # use all the information, not folded and allow monomorphic sites
    sfsvec <- sfs(gm, folded = FALSE, nozero = FALSE)
    ploidy <- gm$geno$ploidy
    N <- length(gm$maps$sample.id)
    xN <- as.integer(names(sfsvec)[length(sfsvec)]) # in number of chromosomes
    stopifnot("SFS entries not matching sample size" = (xN == N * ploidy))
    L <- sum(sfsvec)
    pdfvec <- sfsvec / L

    # chromosome range to predict over
    maxind <- if (is.null(maxind)) N else min(maxind, N)
    nchr <- (1:maxind) * ploidy

    # expected diversity at each subsample size (one list per row of output)
    outlist <- lapply(
        nchr,
        .mutdiv.theory,
        xN = xN,
        kpdf = pdfvec,
        L = L,
        ploidy = ploidy
    )
    # use rbind to avoid importing dplyr (same as MARsampling / MARextinction)
    outdf <- do.call(rbind, lapply(outlist, as.data.frame, stringsAsFactors = FALSE))

    # set outdf as a martheory class
    class(outdf) <- c("martheory", class(outdf)) # MARtheory output class
    return(outdf)
}
