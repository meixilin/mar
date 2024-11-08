# matrix validators
.valid_genotype <- function(genotype, ploidy) {
    # if (ploidy != 2) {
    #     warning(paste0("ploidy = ", ploidy, " != 2, be careful with the interpretation of the results"))
    # }
    # 0 = HomRef, 1 = Het, 2 = HomAlt
    valid_vars <- seq(0, ploidy)
    stopifnot(is.matrix(genotype))
    stopifnot(all(unique(as.vector(genotype)) %in% valid_vars))
    return(invisible())
}

# genotype operations that works on both margeno and SeqArray objects
.get_genodata <- function(x, what = c("sample.id", "variant.id", "position", "chromosome","genotype","ploidy", "num.variant")) {
    what <- match.arg(what)
    if (class(x) == "margeno") {
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

# convert genotype matrix to SeqArray format (3D array)
.mat2seqgeno <- function(gg) {
    stopifnot(class(gg) == "matrix")
    stopifnot(length(dim(gg)) == 2)
    df <- array(dim=c(2, rev(dim(gg))))
    df[1,,] <- apply(gg, 1, function(xx) as.integer(xx > 1))
    df[2,,] <- apply(gg, 1, function(xx) as.integer(xx > 0))
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

# convert margeno object to SeqArray object (TODO: debug)
# inspired by: https://github.com/zhengxwen/SeqArray/issues/62
.margeno2seqarray <- function(x, gds.fn = NULL, opengds = FALSE) {
    stopifnot(class(x) == "margeno")

    # Create a new GDS file
    if(is.null(gds.fn)) {
        gds.fn <- tempfile(fileext = ".gds")
    }

    gds <- gdsfmt::createfn.gds(gds.fn)

    # Add the data
    gdsfmt::add.gdsn(gds, "sample.id", x$sample.id)
    gdsfmt::add.gdsn(gds, "variant.id", x$variant.id)
    gdsfmt::add.gdsn(gds, "position", x$position)
    gdsfmt::add.gdsn(gds, "chromosome", x$chromosome)
    gdsfmt::add.gdsn(gds, "allele", rep("N,N", length(x$variant.id)))
    geno_array <- .mat2seqgeno(x$genotype)
    gdsfmt::add.gdsn(gds, "genotype", geno_array)

    gdsfmt::put.attr.gdsn(gds$root, "FileFormat", "SEQ_ARRAY")
    gdsfmt::put.attr.gdsn(gds$root, "FileVersion", "v1.0")
    gdsfmt::addfolder.gdsn(gds, "description")

    # Close the GDS file
    gdsfmt::closefn.gds(gds)

    # Open as SeqVarGDSClass
    if (opengds) {
        seqfile <- SeqArray::seqOpen(gds.fn)
        return(seqfile)
    } else {
        return(gds.fn)
    }
}

.filter_genosample <- function(x, sample.id) {
    stopifnot(all(sample.id %in% .get_genodata(x, "sample.id")))
    if (class(x) == "margeno") {
        idx = match(sample.id, x$sample.id)
        x$sample.id = x$sample.id[idx]
        x$genotype = as.matrix(x$genotype[,idx])
        return(x)
    } else if (class(x) == "SeqVarGDSClass") {
        return(SeqArray::seqSetFilter(x, sample.id = sample.id))
    } else {
        stop("x must be either a SeqArray or a margeno object")
    }
}

.filter_genovariant <- function(x, variant.id) {
    stopifnot(all(variant.id %in% .get_genodata(x, "variant.id")))
    if (class(x) == "margeno") {
        idx = match(variant.id, x$variant.id)
        x$variant.id = x$variant.id[idx]
        x$position = x$position[idx]
        x$chromosome = x$chromosome[idx]
        x$genotype = as.matrix(x$genotype[idx,])
        return(x)
    } else if (class(x) == "SeqVarGDSClass") {
        return(SeqArray::seqSetFilter(x, variant.id = variant.id))
    } else {
        stop("x must be either a SeqArray or a margeno object")
    }
}


# # set sample and variant filters
# seqSetFilter(f, sample.id=samp.id[c(2,4,6,8)])
# set.seed(100)
# seqSetFilter(f, variant.id=sample(variant.id, 5))

# # get genotypic data
# dim(seqGetData(f, "genotype"))
# seqGetData(f, "genotype")[1:2,1:5, 1:10]

# ## OR
# # set sample and variant filters
# seqSetFilter(f, sample.sel=c(2,4,6,8))
# set.seed(100)
# seqSetFilter(f, variant.sel=sample.int(length(variant.id), 5))

# # get genotypic data
# seqGetData(f, "genotype")

# ploidy <- seqSummary(f)$ploidy

# ## set the intersection

# seqResetFilter(f)
# seqSetFilterChrom(f, 10L)
# seqSummary(f, "genotype", check="none")

# AF <- seqAlleleCount(f)
# AN <- seqNumAllele(f)
# plot(AF)

# table(AF <= 0.9)

# seqSetFilter(f, variant.sel=(AF<=0.9), action="intersect")
# seqSummary(f, "genotype", check="none")
