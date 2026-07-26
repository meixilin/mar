# matrix validators. allows NA and invariant sites.
.valid_genotype <- function(genotype, ploidy) {
    if (!is.matrix(genotype)) {
        stop("genotype must be a matrix")
    }

    valid_vars <- seq(0, ploidy)

    all_vars <- unique(as.vector(genotype))

    # allow NA, but reject invalid allele states
    badvars <- all_vars[!is.na(all_vars) & !(all_vars %in% valid_vars)]

    if (length(badvars) > 0) {
        stop(paste0(
            "genotype contains ",
            length(badvars),
            " unique invalid values: ",
            toString(head(badvars)), "..."
        ))
    }

    return(invisible())
}

.filter_genosample <- function(x, sample.id) {
    stopifnot(all(sample.id %in% .get_genodata(x, "sample.id")))
    if ("margeno" %in% class(x)) {
        idx <- match(sample.id, x$sample.id)
        x$sample.id <- x$sample.id[idx]
        x$genotype <- as.matrix(x$genotype[, idx])
        return(x)
    } else if (is.data.frame(x)) {
        sample_cols <- names(x)[-(1:9)]
        idx <- match(sample.id, sample_cols)
        return(x[, c(1:9, 9 + idx)])
    } else {
        stop("x must be either a data frame (VCF) or a margeno object")
    }
}

.filter_genovariant <- function(x, variant.id) {
    stopifnot(all(variant.id %in% .get_genodata(x, "variant.id")))
    if ("margeno" %in% class(x)) {
        idx <- match(variant.id, x$variant.id)
        x$variant.id <- x$variant.id[idx]
        x$position <- x$position[idx]
        x$chromosome <- x$chromosome[idx]
        x$genotype <- as.matrix(x$genotype[idx, ])
        return(x)
    } else if (is.data.frame(x)) {
        idx <- match(variant.id, x$ID)
        return(x[idx, ])
    } else {
        stop("x must be either a data frame (VCF) or a margeno object")
    }
}

# author Feng Li, Department of Statistics, Stockholm University, Sweden.
# https://github.com/feng-li/flutils/blob/master/R/math__hamonic.R
.Hn <- function(n) {
    out <- digamma(n + 1) - digamma(1)
    return(out)
}
