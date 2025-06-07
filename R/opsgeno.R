# matrix validators. Also validates that there are no invariant sites
.valid_genotype <- function(genotype, ploidy) {
    # if (ploidy != 2) {
    #     warning(paste0("ploidy = ", ploidy, " != 2, be careful with the interpretation of the results"))
    # }
    # 0 = HomRef, 1 = Het, 2 = HomAlt
    valid_vars <- seq(0, ploidy)
    all_vars <- unique(as.vector(genotype))
    AC <- matrixStats::rowSums2(genotype)
    # number of sampled lineages
    xN = ncol(genotype)*ploidy
    if (!is.matrix(genotype)) {
        stop("genotype must be a matrix")
    }
    if (!all(all_vars %in% valid_vars)) {
        badvars = all_vars[!all_vars %in% valid_vars]
        stop(paste0("genotype contains ", length(badvars), " unique invalid values: ", toString(head(badvars)), "..."))
    }
    if (any(AC == 0)) {
        badac = which(AC == 0)
        stop(paste0("genotype is homozygous reference at ", length(badac), " locations: ", toString(head(badac)), "..."))
    }
    if (any(AC == xN)) {
        badac = which(AC == xN)
        stop(paste0("genotype is homozygous alternative at ", length(badac), " locations: ", toString(head(badac)), "..."))
    }
    return(invisible())
}

# genotype operations that works on both margeno and SeqArray objects
.get_genodata <- function(x, what = c("sample.id", "variant.id", "position", "chromosome","genotype","ploidy", "num.variant")) {
    what <- match.arg(what)
    if ("margeno" %in% class(x)) {
        return(switch(what,
            "num.variant" = length(x$variant.id),
            x[[what]]
        ))
    } else if (class(x) == "SeqVarGDSClass") {
        # TODO: nsnps debug
        return(switch(what,
            "num.variant" = SeqArray::seqSummary(x, verbose = FALSE)$num.variant,
            "ploidy" = SeqArray::seqSummary(x, verbose = FALSE)$ploidy,
            "genotype" = .seqgeno2mat(SeqArray::seqGetData(x, what), SeqArray::seqSummary(x, verbose = FALSE)$ploidy),
            SeqArray::seqGetData(x, what)
        ))
    } else {
        stop("x must be either a SeqArray or a margeno object")
    }
}

# convert seqGetData output to genotype matrix in margeno
# in SeqArray, dim1 = ploidy, dim2 = number of sample, dim3 = number of variant
# in margeno$genotype, dim1 = number of variant (row), dim2 = number of sample (col), only diploid or pseudo-haploid allowed
.seqgeno2mat <- function(gg, ploidy) {
    stopifnot(class(gg) == "array")
    stopifnot(length(dim(gg)) == 3 & dim(gg)[1] == ploidy)
    # sum up the allele counts
    df <- apply(gg, c(3,2), sum)
    attr(df, 'dimnames') <- NULL
    .valid_genotype(df, ploidy)
    return(df)
}

# convert seqArray object to margeno object
.seqarray2margeno <- function(x) {
    stopifnot(class(x) == "SeqVarGDSClass")

    # create margeno object
    out <- .new_margeno(
        sample.id = .get_genodata(x, "sample.id"),
        variant.id = .get_genodata(x, "variant.id"),
        position = .get_genodata(x, "position"),
        chromosome = .get_genodata(x, "chromosome"),
        genotype = .get_genodata(x, "genotype"),
        ploidy = .get_genodata(x, "ploidy")
    )
    return(out)
}

# # convert genotype matrix to SeqArray format (3D array)
# .mat2seqgeno <- function(gg) {
#     stopifnot(class(gg) == "matrix")
#     stopifnot(length(dim(gg)) == 2)
#     df <- array(dim=c(2, rev(dim(gg))))
#     df[1,,] <- apply(gg, 1, function(xx) as.integer(xx > 1))
#     df[2,,] <- apply(gg, 1, function(xx) as.integer(xx > 0))
#     return(df)
# }

# # convert margeno object to SeqArray object (TODO: debug)
# # inspired by: https://github.com/zhengxwen/SeqArray/issues/62
# .margeno2seqarray <- function(x, gds.fn = NULL, opengds = FALSE) {
#     stopifnot(class(x) == "margeno")

#     # Create a new GDS file
#     if(is.null(gds.fn)) {
#         gds.fn <- tempfile(fileext = ".gds")
#     }

#     gds <- gdsfmt::createfn.gds(gds.fn)

#     # Add the data
#     gdsfmt::add.gdsn(gds, "sample.id", x$sample.id)
#     gdsfmt::add.gdsn(gds, "variant.id", x$variant.id)
#     gdsfmt::add.gdsn(gds, "position", x$position)
#     gdsfmt::add.gdsn(gds, "chromosome", x$chromosome)
#     gdsfmt::add.gdsn(gds, "allele", rep("N,N", length(x$variant.id)))
#     geno_array <- .mat2seqgeno(x$genotype)
#     gdsfmt::add.gdsn(gds, "genotype", geno_array)

#     gdsfmt::put.attr.gdsn(gds$root, "FileFormat", "SEQ_ARRAY")
#     gdsfmt::put.attr.gdsn(gds$root, "FileVersion", "v1.0")
#     gdsfmt::addfolder.gdsn(gds, "description")

#     # Close the GDS file
#     gdsfmt::closefn.gds(gds)

#     # Open as SeqVarGDSClass
#     if (opengds) {
#         seqfile <- SeqArray::seqOpen(gds.fn)
#         return(seqfile)
#     } else {
#         return(gds.fn)
#     }
# }

.filter_genosample <- function(x, sample.id) {
    stopifnot(all(sample.id %in% .get_genodata(x, "sample.id")))
    if ("margeno" %in% class(x)) {
        idx = match(sample.id, x$sample.id)
        x$sample.id = x$sample.id[idx]
        x$genotype = as.matrix(x$genotype[,idx])
        return(x)
    } else if (class(x) == "SeqVarGDSClass") {
        SeqArray::seqSetFilter(x, sample.id = sample.id)
        return(x)
    } else {
        stop("x must be either a SeqArray or a margeno object")
    }
}

.filter_genovariant <- function(x, variant.id) {
    stopifnot(all(variant.id %in% .get_genodata(x, "variant.id")))
    if ("margeno" %in% class(x)) {
        idx = match(variant.id, x$variant.id)
        x$variant.id = x$variant.id[idx]
        x$position = x$position[idx]
        x$chromosome = x$chromosome[idx]
        x$genotype = as.matrix(x$genotype[idx,])
        return(x)
    } else if (class(x) == "SeqVarGDSClass") {
        SeqArray::seqSetFilter(x, variant.id = variant.id)
        return(x)
    } else {
        stop("x must be either a SeqArray or a margeno object")
    }
}

#' Harmonic number
#'
#'
#' @title Harmonic number
#' @param n "non-negative numeric" Non integer is also allowed.
#' @return "numeric"
#' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
#' @note Created: Wed May 30 10:51:32 CEST 2012;
#'       Current: Wed May 30 10:51:37 CEST 2012.
#'       Source: https://github.com/feng-li/flutils/blob/master/R/math__hamonic.R
.Hn <- function(n)
{
    out <- digamma(n+1) - digamma(1)
    return(out)
}

