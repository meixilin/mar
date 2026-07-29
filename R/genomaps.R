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
    class(obj) <- c(class(obj), "margeno")
    return(obj)
}

# marmaps class
# create marmaps object
# resolution determination inspired from: https://www.sciencedirect.com/science/article/pii/S0098300405002657 (eq. 12)
.lonlat_res <- function(lonlat, lonlatr) {
    aa <- apply(lonlatr, 2, diff)
    area <- aa[1] * aa[2]
    mapres <- 0.5 * sqrt(area / nrow(lonlat))
    # round to first non-zero digits
    out <- as.numeric(sprintf("%.e", mapres))
    return(out)
}

# these two functions can be used to create `raster_samples` equivalent
.lonlat_raster <- function(lonlat, lonlatr, mapres, mapcrs) {
    # raster ext respects only xmin and ymax when resolution is specified
    # calculate the extent of the raster
    xmin <- min(lonlatr[, 1]) - 0.5 * mapres
    ymax <- max(lonlatr[, 2]) + 0.5 * mapres
    xmax <- xmin + (ceiling(diff(lonlatr[, 1]) / mapres) + 1) * mapres
    ymin <- ymax - (ceiling(diff(lonlatr[, 2]) / mapres) + 1) * mapres
    baser <- terra::rast(
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        resolution = mapres,
        crs = mapcrs
    )
    pts <- terra::vect(lonlat, type = "points", crs = mapcrs)
    rr <- terra::rasterize(pts, baser, fun = length)
    return(rr)
}

# Constructor for marmaps class
.new_marmaps <- function(sample.id, lonlat, samplemap, cellid) {
    obj <- list(
        sample.id = sample.id,
        lonlat = lonlat,
        samplemap = samplemap,
        cellid = cellid
    )
    class(obj) <- c(class(obj), "marmaps")
    return(obj)
}

# genomaps class
# combine the margeno and marmaps objects
.new_genomaps <- function(geno, maps) {
    obj <- list(
        geno = geno,
        maps = maps
    )
    class(obj) <- c(class(obj), "genomaps")
    attr(obj, "genolen") <- .get_genodata(geno, "num.variant")
    return(obj)
}

#' Create a margeno object
#'
#' @param sample.id Character or integer or numeric vector of unique sample IDs
#' @param variant.id Integer vector of unique variant IDs
#' @param position Integer vector of variant positions. Or NULL.
#' @param chromosome Character or integer vector of chromosome IDs. Or NULL.
#' @param genotype Matrix of genotypes, where rows represent samples and columns represent variants. Each value represents the number of alternative alleles.
#' @param ploidy Numeric value of ploidy. *Warning* ploidy other than 2 can be used, but result interpretation is not guaranteed.
#'
#' @return A margeno object
#' @export
#'
#' @examples
#' sample_id <- c("sample1", "sample2", "sample3")
#' variant_id <- 1:4
#' position <- as.integer(c(100, 200, 300, 400)) # position has to be integer
#' chromosome <- c("chr1", "chr1", "chr1", "chr2")
#' position <- NULL
#' chromosome <- NULL # if no position or chromosome
#' genotype <- matrix(c(1, 2, 0, 0, 1, 1, 2, 0, 1, 2, 2, 1), nrow = 4, ncol = 3)
#' ploidy <- 2
#' margeno(sample_id, variant_id, position, chromosome, genotype, ploidy)
margeno <- function(sample.id, variant.id, position, chromosome, genotype, ploidy) {
    # validate data class
    stopifnot("sample id must be character, integer, or numeric" = class(sample.id) %in% c("character", "integer", "numeric"))
    stopifnot("variant id must be an integer" = class(variant.id) == "integer")
    stopifnot("postition must be an integer or null" = class(position) == "integer" | is.null(position))
    stopifnot("chromosome must be chacter, integer, or null" = class(chromosome) %in% c("character", "integer") | is.null(chromosome))
    stopifnot("genotype must be a matrix" = "matrix" %in% class(genotype))
    stopifnot("ploidy must be a numeric" = class(ploidy) == "numeric")

    # validate data dimensions
    stopifnot("sample ids must be unique" = anyDuplicated(sample.id) == 0)
    stopifnot("variant ids must be unique" = anyDuplicated(variant.id) == 0)
    stopifnot("variant id dimensions do not match genotype matrix" = length(variant.id) == dim(genotype)[1])
    stopifnot("sample id dimensions do not match genotype matrix" = length(sample.id) == dim(genotype)[2])

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

#' Create a marmaps object
#'
#' @param lonlatdf A data frame with three columns: sample ID, longitude, and latitude.
#' @param mapres Optional numeric value for the map resolution. If not provided (mapres = NULL), it will be automatically calculated.
#' @param mapcrs A character string specifying the coordinate reference system (CRS) for the map.
#'
#' @return A marmaps object.
#' @export
#'
#' @examples
#' lonlatdf <- data.frame(
#'     id = c("sample1", "sample2", "sample3"),
#'     longitude = c(-122.1, -121.9, -122.0),
#'     latitude = c(37.8, 37.7, 37.9),
#'     stringsAsFactors = FALSE
#' )
#' marmaps(lonlatdf, mapres = NULL, mapcrs = "EPSG:4326")
marmaps <- function(lonlatdf, mapres, mapcrs) {
    # unpack lonlatdf
    stopifnot(class(lonlatdf) == "data.frame" & ncol(lonlatdf) == 3)
    sample.id <- lonlatdf[[1]]
    lonlat <- as.matrix(lonlatdf[, 2:3])

    # Validate inputs
    stopifnot("sample ids must be a character, integer, or numeric" = class(sample.id) %in% c("character", "integer", "numeric"))
    stopifnot("lonlat must be a matrix" = is.matrix(lonlat))
    .valid_lonlat(lonlat)
    stopifnot("sample id dimensions must match lonlat dimensions" = length(sample.id) == nrow(lonlat))
    stopifnot("mapres must be a numeric or null" = is.null(mapres) | is.numeric(mapres))
    stopifnot("mapcrs must be a character" = is.character(mapcrs))

    # Calculate map resolution if not provided
    lonlatr <- apply(lonlat, 2, range)
    if (is.null(mapres)) {
        mapres <- .lonlat_res(lonlat, lonlatr)
    }

    # Create sample map
    samplemap <- .lonlat_raster(lonlat, lonlatr, mapres, mapcrs)

    # Get cell IDs
    cellid <- terra::cellFromXY(samplemap, lonlat)
    stopifnot("lonlat outside range" = !any(is.na(cellid))) # stop if lonlat outside of raster

    # Create object using constructor
    output <- .new_marmaps(
        sample.id = sample.id,
        lonlat = lonlat,
        samplemap = samplemap,
        cellid = cellid
    )

    message("number of samples: ", length(sample.id))
    message("longitude range: [", lonlatr[1, 1], ", ", lonlatr[2, 1], "]")
    message("latitude range: [", lonlatr[1, 2], ", ", lonlatr[2, 2], "]")
    message("map resolution: ", mapres)

    return(output)
}

#' Create a genomaps object
#'
#' The genomaps object is the central data structure for the mar package. It combines the genetic data
#' from the margeno object and the spatial information from the marmaps object, providing a
#' comprehensive representation of the genetic and geographic data, connected by the sample IDs.
#'
#' @param geno A margeno object.
#' @param maps A marmaps object.
#'
#' @return A genomaps object.
#' @export
#'
#' @examples
#' sample_id <- c("sample1", "sample2", "sample3")
#' variant_id <- 1:4
#' position <- as.integer(c(100, 200, 300, 400)) # position has to be integer
#' chromosome <- c("chr1", "chr1", "chr1", "chr2")
#' position <- NULL
#' chromosome <- NULL # if no position or chromosome
#' genotype <- matrix(c(1, 2, 0, 0, 1, 1, 2, 0, 1, 2, 2, 1), nrow = 4, ncol = 3)
#' ploidy <- 2
#' geno <- margeno(sample_id, variant_id, position, chromosome, genotype, ploidy)
#'
#' lonlatdf <- data.frame(
#'     id = sample_id,
#'     longitude = c(-122.1, -121.9, -122.0),
#'     latitude = c(37.8, 37.7, 37.9),
#'     stringsAsFactors = FALSE
#' )
#' maps <- marmaps(lonlatdf, mapres = NULL, mapcrs = "EPSG:4326")
#'
#' genomaps(geno, maps)
genomaps <- function(geno, maps) {
    # validate inputs
    stopifnot("geno must be a margeno object" = "margeno" %in% class(geno))
    stopifnot("maps must be a marmaps object" = "marmaps" %in% class(maps))
    stopifnot("geno and maps ids must match" = all(geno$sample.id == maps$sample.id))

    # create object
    obj <- .new_genomaps(geno, maps)
    return(obj)
}
