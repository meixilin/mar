# fit sad models (sorted alphabetically)
.sad_models <- c("bs", "geom", "lnorm", "ls", "mzsm", "weibull")

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

# allele counts
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

# not allowing filter at this stage. All filters should be done at step "gm"
#' Title
#'
#' @param AC
#' @param N
#' @param ploidy
#' @param folded
#' @param nozero
#'
#' @return
#' @export
#'
#' @examples
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

# generate expected SFS for given SNPs surveyed (lenAC) and number of samples (N)
#' Title
#'
#' @param lenAC
#' @param N
#' @param ploidy
#' @param folded
#' @param nozero
#'
#' @return
#' @export
#'
#' @examples
expsfs <- function(lenAC, N, ploidy, folded = TRUE, nozero = TRUE) {
    xN = N*ploidy
    theta = lenAC / Hn(xN) # scale theta
    expsfs = c(0, theta/(1:xN)) # need to add 0 as xN is the same as zero when folded
    expsfs <- .new_sfs(expsfs, folded, nozero)
    return(expsfs)
}

# convert a SFS back to AC vector (rank abundance data to some extent)
# NOT USED AND BUGGY
.sfs2AC <- function(vect) {
    M = round(sum(vect))
    oAC <- vector()
    for(ii in seq_along(vect)) {
        binii = as.integer(names(vect)[ii])
        nii = round(vect[ii])
        if (nii > 0) {
            oAC <- c(oAC, rep(binii, times = nii))
        }
    }

    # fill in AC if needed
    if (length(oAC) < M) {
        iis = 1
        warning(any(round(vect) == 0)) # should have at least some entried rounded out
        while(length(oAC) < M) {
            svect <- sort(vect[round(vect) == 0], decreasing = TRUE)
            binii_s <- as.integer(names(svect)[iis])
            oAC <- c(oAC, binii_s)
            iis = iis + 1
        }
    }
    return(oAC)
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

# < 0 and NA test have been tested in the sfs data construction phase
#' Title
#'
#' @param model
#' @param data
#' @param missing_model_cutoff
#'
#' @return
#' @export
#'
#' @examples
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
