.new_genemaps <- function(lonlat,
                          sample,
                          genotype,
                          genoinfo) {
    obj <-
        list(
            lonlat = lonlat,
            sample = sample,
            genotype = genotype,
            genoinfo = genoinfo
        )
    class(obj) <- "genemaps"

    # add attribute for the starting extent
    lonlatr <- c(
        range(lonlat[,1]),
        range(lonlat[,2])
    )
    attr(obj, "extent") <- lonlatr
    return(obj)
}

# validate and helper functions for genemaps class
genemaps <- function(lonlat,
                     genotype,
                     sample = NULL,
                     genoinfo = NULL) {
    # convert lonlat and genotype if needed
    # TODO: big genotype data manipulations
    if (class(lonlat) != "matrix")
        lonlat = as.matrix(lonlat)
    if (class(genotype) != "matrix")
        genotype = as.matrix(genotype)

    # check dimensions
    stopifnot(dim(lonlat)[2] == 2 & dim(lonlat) > 0)
    # not allowing NA values in lonlat
    stopifnot(!any(is.na(lonlat)))
    # set variables
    nsamples = dim(lonlat)[1]
    nsnps = dim(genotype)[2]
    logger::log_info("number of samples: ", nsamples)
    logger::log_info("number of genomic sites: ", nsnps)

    # fill in sample and genoinfo
    if (is.null(sample)) {
        sample = data.frame(id = 1:nsamples)
    }
    if (is.null(genoinfo)) {
        genoinfo = data.frame(chr = NA, pos = 1:nsnps)
    }

    # check for number of sample dimensions
    stopifnot(nsamples == dim(genotype)[1] &
                  nsamples == dim(sample)[1])
    stopifnot(nsnps == dim(genotype)[2] & nsnps == dim(genoinfo)[1])

    # assemble genemaps
    output <- .new_genemaps(
        lonlat = lonlat,
        sample = sample,
        genotype = genotype,
        genoinfo = genoinfo
    )
    return(output)
}
