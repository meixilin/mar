# fit sad models (sorted alphabetically)
.sad_models <- c("bs", "geom", "lnorm", "ls", "mzsm", "weibull")

#' Fit Species Abundance Distribution (SAD) Models to Genotype Data
#'
#' @param gm A genotype matrix object containing genetic data
#' @param sad_models Vector of SAD model names to fit. Default uses internal mar:::.sad_models
#' @param predict Logical, whether to predict SFS from fitted SAD models. Default TRUE
#' @param folded Logical, whether to fold the SFS. Default TRUE
#'
#' @return A list with S3 class "marsad" containing:
#'   \item{sadms}{List of fitted SAD models}
#'   \item{AICtabs}{AIC table comparing model fits}
#'   \item{sadsfss}{Predicted site frequency spectra if predict=TRUE}
#' @export
#'
#' @examples
#' # Fit SAD models to genotype data in gm1001g
#' sad_fit <- MARsad(gm1001g)
MARsad <- function(gm, sad_models = .sad_models, predict = TRUE, folded = TRUE) {
    AC <- .get_AC(gm$geno)
    N <- length(gm$maps$sample.id)
    ploidy <- gm$geno$ploidy
    sadms <- lapply(sad_models, function(x) sads::fitsad(AC, x))
    names(sadms) <- sad_models
    AICtabs <- bbmle::AICtab(sadms, base = TRUE, logLik = TRUE, mnames = sad_models)
    if (predict) {
        sadsfss <- lapply(sadms, function(x) .sadpred(x, N, ploidy, folded))
    } else {
        sadsfss <- NULL
    }
    output <- list(sadms = sadms, AICtabs = AICtabs, sadsfss = sadsfss)
    class(output) <- c(class(output), "marsad")
    return(output)
}

# TODO: I used pfunc() instead of qfunc() as described in sads package as the interest is in SFS.
# Not sure if it makes sense
.sadpred <- function(sadm, N, ploidy, folded) {
    sad <- sadm@sad
    xN <- N*ploidy
    S <- length(sadm@data$x) # length of AC
    J <- sum(sadm@data$x) # total number of individuals
    mycoef <- as.list(bbmle::coef(sadm))
    psad <- switch(sad,
        bs = sads::pbs,
        geom = stats::pgeom,
        lnorm = stats::plnorm,
        ls = sads::pls,
        mzsm = sads::pmzsm,
        weibull = stats::pweibull
    )
    plist <- do.call(psad, c(list(q = 1:(xN-1)), mycoef)) # q = 1:(xN-1)
    pbins <- c(plist,1) - c(0,plist)
    stopifnot(sum(pbins) == 1) # sanity check
    raw_sfs <- c(0, pbins * S) # need to add zero as xN is the same as zero when folded
    sadsfs <- .new_sfs(raw_sfs, folded, nozero = TRUE)
    return(sadsfs)
}

.pipe_sadsfs <- function(gm, marsad, genosfs, folded) {
    AC <- .get_AC(gm$geno)
    N <- length(gm$maps$sample.id)
    ploidy <- gm$geno$ploidy
    neutralsfs <- expsfs(lenAC = length(AC), N = N, ploidy = ploidy, folded = folded)
    allsfs <- list(genosfs, neutralsfs)
    names(allsfs) <- c("data", "neutral")
    # if SAD predicted
    if (!is.null(marsad$sadsfss)) {
        allsfs <- c(allsfs, marsad$sadsfss)
    }
    # compare by logLik
    ll_list <- sapply(allsfs, function(model) ll_sfs(model = model, data = genosfs))
    statdf <- data.frame(model = names(allsfs),
                         logLik = unname(ll_list),
                         stringsAsFactors = FALSE)
    # return list of statdf and allsfs
    output <- list(sfs = allsfs, # list of sfs class objects
                   statdf = statdf)
    return(output)
}

# generate per-cell genotype table
.genotype_bycell <- function(gm) {
    cellids = gm$maps$cellid
    ploidy = gm$geno$ploidy
    geno = sapply(unique(cellids), function(g) {
        matrixStats::rowMaxs(gm$geno$genotype, cols = which(cellids == g))
    })
    return(geno)
}

# allele counts (revert to earlier)
.get_AC <- function(gg) {
    AC <- matrixStats::rowSums2(gg$genotype)
    # stop if there are any NAs or warn if fully zero ACs (not a SNP in this dataset)
    stopifnot(all(!is.na(AC)))
    if (any(AC == 0)) {
        warning(paste0("There are ", sum(AC == 0)," invariant sites in the genotype matrix"))
        AC <- AC[AC != 0]
    }
    return(AC)
}

# SFS operations
.foldsfs <- function(vect) {
    flen = floor(length(vect)/2)
    fvect = (vect + rev(vect))[1:flen]
    if(length(vect) %% 2 == 1) {
        fvect = c(fvect, vect[flen+1])
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
        names(vect) <- 0:(length(vect)-1)
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
#' allele_counts <- c(1,1,0,2,0,1,1,0,0,2,30)
#' sfs_result <- sfs(allele_counts, N=50, ploidy=2)

sfs <- function(AC, N, ploidy, folded = TRUE, nozero = TRUE) {
    xN = N*ploidy
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
#' exp_sfs <- expsfs(lenAC=1000, N=100, ploidy=2)

expsfs <- function(lenAC, N, ploidy, folded = TRUE, nozero = TRUE) {
    xN = N*ploidy
    theta = lenAC / .Hn(xN) # scale theta
    expsfs = c(0, theta/(1:xN)) # need to add 0 as xN is the same as zero when folded
    expsfs <- .new_sfs(expsfs, folded, nozero)
    return(expsfs)
}

# sfs list to dataframe
.sfsl2df <- function(sfsl) {
    outdf = data.frame(matrix(nrow = length(sfsl[[1]]), ncol = length(sfsl) + 1), stringsAsFactors = FALSE)
    colnames(outdf) = c("AC", names(sfsl))
    outdf[,1] = as.integer(names(sfsl[[1]]))
    for (ii in seq_along(sfsl)) {
        outdf[, ii+1] = as.vector(sfsl[[ii]])
    }
    return(outdf)
}

# adapt it from dadi.Inference.ll (Original function available at Gutenkust et al. 2009)
..ll_per_bin <- function(model, data) {
    if (data == 0 | model == 0) {
        out = 0
    } else {
        out = - model + log(model) * data - lgamma(data + 1)
    }
    return(out)
}

.ll_per_bin <- Vectorize(..ll_per_bin)

#' Calculate Log-Likelihood Between Model and Data SFS
#'
#' @param model An object of class "sfs" containing model predictions
#' @param data An object of class "sfs" containing observed data
#' @param missing_model_cutoff Threshold for warning about SFS entries that are zeros in the model but non-zero in the data. Default 1e-6.
#'
#' @return Log-likelihood value comparing model to data
#' @export
#'
#' @examples
#' # Calculate log-likelihood between model and data SFS
#' model_sfs <- expsfs(lenAC=1000, N=50, ploidy=2)
#' data_sfs <- sfs(AC=c(1,1,0,2,0,1,1,0,0,2,30), N=50, ploidy=2)
#' ll <- ll_sfs(model_sfs, data_sfs)
ll_sfs <- function(model, data, missing_model_cutoff = 1e-6) {
    stopifnot("sfs" %in% c(class(model), class(data)))
    stopifnot(length(model) == length(data)) # same length
    stopifnot(all(names(model) == names(data))) # same entries
    # check for zeros in models
    d0 = data[which(model == 0)]
    if (sum(d0)/sum(data) > missing_model_cutoff) {
        warning(paste0("In ", sprintf("%.2f%%", 100 * sum(d0) / sum(data)), " of data. Model is 0 where data is neither masked nor 0."))
    }
    ll <- sum(.ll_per_bin(model, data))
    return(ll)
}
