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
    lonlatr <- apply(lonlat, 2, range)
    attr(obj, "extent") <- lonlatr
    # add attribute to recommend operation resolution
    attr(obj, "res") <- .lonlat_res(lonlat, lonlatr)
    print(attributes(obj))
    return(obj)
}

# resolution determination inspired from: https://www.sciencedirect.com/science/article/pii/S0098300405002657 (eq. 12)
.lonlat_res <- function(lonlat, lonlatr) {
    aa = apply(lonlatr, 2, diff)
    area = aa[1]*aa[2]
    res = 0.5*sqrt(area/nrow(lonlat))
    # round to first non-zero digits
    outres = as.numeric(sprintf('%.e', res))
    return(outres)
}

# these two functions can be used to create `raster_samples` equivalent
# TODO: currently decided to not use this framework (constrains flexibility and speed)
.raster_lonlat <- function(lonlat, lonlatr, res = NULL, padding = 0.01) {
    # if resolution not provided, determine resolution
    if(is.null(res)) {
        res = .lonlat_res(lonlat, lonlatr)
    }
    lonlatrp <- t(apply(lonlatr, 2, function(xx){xx + c(-1,1)*diff(xx)*padding}))
    baser <- raster(resolution = res, extent(lonlatrp))
    rr <- rasterize(lonlat, baser, fun = "count")
    print(rr)
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

