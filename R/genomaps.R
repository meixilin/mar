# margeno class
.new_margeno <- function(sample.id, variant.id, position, chromosome, genotype, ploidy) {
    # create object (order determined by SeqArray output)
    obj <- list(
        sample.id = sample.id,
        variant.id = variant.id,
        position = position,
        chromosome = chromosome,
        genotype = genotype,
        ploidy = ploidy
    )
    class(obj) <- "margeno"
    return(obj)
}

margeno <- function(sample.id, variant.id, position, chromosome, genotype, ploidy) {
    # validate data class
    stopifnot(class(sample.id) %in% c("character", "integer"))
    stopifnot(class(variant.id) == "integer")
    stopifnot(class(position) == "integer" | is.null(position))
    stopifnot(class(chromosome) %in% c("character", "integer") | is.null(position))
    stopifnot(class(genotype) == "matrix")
    stopifnot(class(ploidy) == "numeric")

    # validate data dimensions
    stopifnot(anyDuplicated(sample.id) == 0)
    stopifnot(anyDuplicated(variant.id) == 0)
    stopifnot(length(variant.id) == dim(genotype)[1])
    stopifnot(length(sample.id) == dim(genotype)[2])

    # validate genotype
    .valid_genotype(genotype, ploidy)

    # create object
    output <- .new_margeno(
        sample.id = sample.id,
        variant.id = variant.id,
        position = position,
        chromosome = chromosome,
        genotype = genotype,
        ploidy = ploidy
    )

    message("number of samples: ", length(sample.id))
    message("number of genomic sites: ", length(variant.id))

    return(output)
}


# marmaps class
# create marmaps object
# resolution determination inspired from: https://www.sciencedirect.com/science/article/pii/S0098300405002657 (eq. 12)
.lonlat_res <- function(lonlat, lonlatr) {
    aa = apply(lonlatr, 2, diff)
    area = aa[1]*aa[2]
    mapres = 0.5*sqrt(area/nrow(lonlat))
    # round to first non-zero digits
    out = as.numeric(sprintf('%.e', mapres))
    return(out)
}

# these two functions can be used to create `raster_samples` equivalent
.lonlat_raster <- function(lonlat, lonlatr, mapres, mapcrs) {
    # raster ext respects only xmin and ymax when resolution is specified
    # calculate the extent of the raster
    xmin = min(lonlatr[,1])-0.5*mapres
    ymax = max(lonlatr[,2])+0.5*mapres
    xmax = xmin + (ceiling(diff(lonlatr[,1])/mapres)+1)*mapres
    ymin = ymax - (ceiling(diff(lonlatr[,2])/mapres)+1)*mapres
    baser <- raster::raster(xmn = xmin, xmx = xmax,
                            ymn = ymin, ymx = ymax,
                            res = mapres,
                            crs = mapcrs)
    rr <- raster::rasterize(lonlat, baser, fun = "count")
    return(rr)
}

# define the plotting method for marmaps
plot.marmaps <- function(marmaps) {
    plot(raster::extent(marmaps$samplemap), xlab = 'lon', ylab = 'lat')
    plot(marmaps$samplemap, add = T)
    points(marmaps$lonlat)
    return(invisible())
}

# Constructor for marmaps class
.new_marmaps <- function(sample.id, lonlat, samplemap, cellid) {
    obj <- list(
        sample.id = sample.id,
        lonlat = lonlat,
        samplemap = samplemap,
        cellid = cellid
    )
    class(obj) <- "marmaps"
    return(obj)
}

# Main marmaps function
marmaps <- function(lonlatdf, mapres, mapcrs) {
    # unpack lonlatdf
    stopifnot(class(lonlatdf) == "data.frame" & ncol(lonlatdf) == 3)
    sample.id <- lonlatdf[[1]]
    lonlat <- as.matrix(lonlatdf[,2:3])

    # Validate inputs
    stopifnot(class(sample.id) %in% c("character", "integer"))
    stopifnot(is.matrix(lonlat))
    .valid_lonlat(lonlat)
    stopifnot(length(sample.id) == nrow(lonlat))
    stopifnot(is.null(mapres) | is.numeric(mapres))
    stopifnot(is.character(mapcrs))

    # Calculate map resolution if not provided
    lonlatr <- apply(lonlat, 2, range)
    if (is.null(mapres)) {
        mapres <- .lonlat_res(lonlat, lonlatr)
    }

    # Create sample map
    samplemap <- .lonlat_raster(lonlat, lonlatr, mapres, mapcrs)

    # Get cell IDs
    cellid <- raster::cellFromXY(samplemap, lonlat)
    stopifnot(!any(is.na(cellid))) # stop if lonlat outside of raster

    # Create object using constructor
    output <- .new_marmaps(
        sample.id = sample.id,
        lonlat = lonlat,
        samplemap = samplemap,
        cellid = cellid
    )

    message("number of samples: ", length(sample.id))
    message("longitude range: [", lonlatr[1,1], ", ", lonlatr[2,1], "]")
    message("latitude range: [", lonlatr[1,2], ", ", lonlatr[2,2], "]")
    message("map resolution: ", mapres)

    return(output)
}

# genomaps class

# combine the margeno and marmaps objects
.new_genomaps <- function(geno, maps) {
    obj <- list(geno = geno,
                maps = maps)
    class(obj) <- "genomaps"
    attr(obj, "genolen") <- .get_genodata(geno, "num.variant")
    return(obj)
}

genomaps <- function(geno, maps) {
    # validate inputs
    stopifnot(class(geno) == "margeno")
    stopifnot(class(maps) == "marmaps")
    stopifnot(all(geno$sample.id == maps$sample.id))

    # create object
    obj <- .new_genomaps(geno, maps)
    return(obj)
}

