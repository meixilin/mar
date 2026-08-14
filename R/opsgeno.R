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
            toString(utils::head(badvars)), "..."
        ))
    }

    return(invisible())
}

# author Feng Li, Department of Statistics, Stockholm University, Sweden.
# https://github.com/feng-li/flutils/blob/master/R/math__hamonic.R
.Hn <- function(n) {
    out <- digamma(n + 1) - digamma(1)
    return(out)
}
