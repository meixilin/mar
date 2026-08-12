# margeno class
.new_margeno <- function(sample.id, variant.id, position, chromosome, genotype, ploidy) {
    # create object (order determined by SeqArray output)
    # 2026 update: add an allele count entry so sfs, .calc_theta can all use (if all NA, AC is 0)
    AC = matrixStats::rowSums2(genotype, na.rm = TRUE)
    obj <- list(
        sample.id = sample.id,
        variant.id = variant.id,
        position = position,
        chromosome = chromosome,
        genotype = genotype,
        allele_count = AC,
        ploidy = ploidy
    )
    class(obj) <- c(class(obj), "margeno")
    return(obj)
}

# marmaps class
# create marmaps object

# reproject coordinates from the input CRS (incrs) into the map CRS (mapcrs)
# so that resolution, extent, and cell assignment are all computed in map units
.reproject_lonlat <- function(lonlat, incrs, mapcrs) {
    # skip when the coordinates are already in the target CRS
    if (identical(incrs, mapcrs)) {
        return(lonlat)
    }
    pts <- terra::vect(lonlat, type = "points", crs = incrs)
    pts <- terra::project(pts, mapcrs)
    out <- terra::crds(pts)
    colnames(out) <- colnames(lonlat)
    # guard against coordinates that fall outside the target projection
    stopifnot(all(is.finite(out)))
    return(out)
}

# resolution determination inspired from: https://www.sciencedirect.com/science/article/pii/S0098300405002657 (eq. 12)
# NOTE: `lonlat`/`lonlatr` must already be expressed in the map CRS units (see
# .reproject_lonlat), otherwise the returned resolution will be in the wrong
# units once the map projection differs from the input coordinate CRS.
.lonlat_res <- function(lonlat, lonlatr) {
    aa <- apply(lonlatr, 2, diff)
    area <- aa[1] * aa[2]
    mapres <- 0.5 * sqrt(area / nrow(lonlat))
    # round to first non-zero digits
    out <- as.numeric(sprintf("%.e", mapres))
    return(out)
}

# create raster based on lonlat object
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
    wrr <- terra::wrap(rr, proxy = FALSE)
    return(wrr)
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

# accessor for the samplemap raster stored in a marmaps object.
# samplemap is stored wrapped (PackedSpatRaster) so the object can be
# serialized (save/saveRDS); unwrap it back to a live SpatRaster on read.
# tolerant of already-unwrapped rasters for backward compatibility.
.get_samplemap <- function(mm) {
    sm <- mm$samplemap
    if (inherits(sm, "PackedSpatRaster")) terra::unwrap(sm) else sm
}

# genomaps class
# combine the margeno and marmaps objects
.new_genomaps <- function(geno, maps) {
    obj <- list(
        geno = geno,
        maps = maps
    )
    class(obj) <- c(class(obj), "genomaps")
    return(obj)
}

#' Create a margeno object
#'
#' @param genotype Matrix of genotypes, where rows represent samples and columns represent variants. Each value represents the number of alternative alleles.
#' @param ploidy Integer. *Warning* ploidy other than 2 can be used, but result interpretation is not guaranteed.
#' @param sample.id Character or integer or numeric vector of unique sample IDs.
#'   Or NULL (the default), in which case samples are numbered from the genotype columns.
#' @param variant.id Integer vector of unique variant IDs.
#'   Or NULL (the default), in which case variants are numbered from the genotype rows.
#' @param position Integer vector of variant positions. Or NULL (the default).
#' @param chromosome Character or integer vector of chromosome IDs. Or NULL (the default).
#'
#' @return A margeno object. In addition to the inputs, it stores `allele_count`,
#'   the number of alternative alleles observed at each site across all samples.
#' @seealso [text_parser()] and [vcf_parser()], which build a margeno object
#'   directly from a file, and [genomaps()] to pair it with a [marmaps()] object.
#' @export
#'
#' @examples
#' # the genotypes behind the 1001 genomes example object
#' geno <- gm1001g$geno
#'
#' # only genotype and ploidy are required
#' margeno(geno$genotype, geno$ploidy)
#'
#' # sample.id, variant.id, position and chromosome are all optional
#' margeno(geno$genotype, geno$ploidy,
#'     sample.id = geno$sample.id,
#'     variant.id = geno$variant.id,
#'     position = geno$position,
#'     chromosome = geno$chromosome
#' )
margeno <- function(genotype, ploidy, sample.id = NULL, variant.id = NULL, position = NULL, chromosome = NULL) {
    # validate data class
    stopifnot(class(sample.id) %in% c("character", "integer", "numeric") | is.null(sample.id))
    stopifnot(class(variant.id) == "integer" | is.null(variant.id))
    stopifnot(class(position) == "integer" | is.null(position))
    stopifnot(class(chromosome) %in% c("character", "integer") | is.null(chromosome))
    stopifnot("matrix" %in% class(genotype))
    stopifnot(class(ploidy) == "numeric")

    # if not provided, fill in sample.id and variant.id
    if (is.null(sample.id)) sample.id <- seq_len(ncol(genotype))
    if (is.null(variant.id)) variant.id <- seq_len(nrow(genotype))
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
#' Rasterizes the sampling coordinates of a species onto a regular grid. The grid
#' is the unit of area for every downstream MAR calculation, so the map's coordinate
#' reference system should be an equal-area projection whenever the samples span a wide latitudinal
#' range.
#'
#' @param lonlatdf A data frame with three columns: sample ID, longitude, and latitude. Column names are not checked.
#' @param mapcrs A character string specifying the coordinate reference system (CRS) for the output map.
#' @param mapres Optional numeric value for the output map resolution, in the units of `mapcrs`.
#'   If not provided (mapres = NULL), it will be automatically calculated.
#' @param incrs A character string specifying the CRS of the input coordinates in `lonlatdf`.
#'   Defaults to "EPSG:4326" (WGS84 longitude/latitude). When `incrs` differs from `mapcrs`,
#'   the coordinates are reprojected into `mapcrs` before the resolution, extent, and cell
#'   assignment are computed, so those stay consistent with the map units.
#'   Set both `incrs` and `mapcrs` to `""` for coordinates on an arbitrary plane,
#'   such as those of a spatial simulation.
#'
#' @return A marmaps object, holding the sample IDs, the coordinates in `mapcrs`,
#'   the sample map raster (wrapped with `terra::wrap()` so that it can be saved,
#'   unwrap it with `terra::unwrap()` before use), and the map cell each sample
#'   falls in.
#' @seealso [lonlat_parser()], which builds a marmaps object directly from a
#'   file, and [genomaps()] to pair it with a [margeno()] object.
#' @export
#'
#' @examples
#' # coordinates of the 1001 genomes accessions, in WGS84 lon/lat
#' lonlat.fn <- system.file("extdata", "1001g_lonlat.txt", package = "mar")
#' lonlatdf <- read.delim(lonlat.fn)
#' head(lonlatdf)
#'
#' # reprojected onto an equal-area map, resolution chosen automatically
#' marmaps(lonlatdf, mapcrs = "EPSG:8857", incrs = "EPSG:4326")
#'
#' # the same map at a coarser, fixed resolution of 200 km
#' marmaps(lonlatdf, mapcrs = "EPSG:8857", mapres = 2e5)
#'
#' # coordinates already in the target projection are not reprojected
#' marmaps(lonlatdf, mapcrs = "EPSG:4326", incrs = "EPSG:8857")
marmaps <- function(lonlatdf, mapcrs, mapres = NULL, incrs = "EPSG:4326") {
    # unpack lonlatdf
    stopifnot(class(lonlatdf) == "data.frame" & ncol(lonlatdf) == 3)
    sample.id <- lonlatdf[[1]]
    lonlat <- as.matrix(lonlatdf[, 2:3])

    # Validate inputs
    stopifnot("sample ids must be a character, integer, or numeric" = class(sample.id) %in% c("character", "integer", "numeric"))
    stopifnot("lonlat must be a matrix" = is.matrix(lonlat))
    .valid_lonlat(lonlat)
    stopifnot(length(sample.id) == nrow(lonlat))
    stopifnot(is.null(mapres) | is.numeric(mapres))
    stopifnot(is.character(mapcrs))
    stopifnot(is.character(incrs))

    # Reproject coordinates into the map CRS so that resolution, extent, and
    # cell assignment are all computed consistently in the map units
    lonlat <- .reproject_lonlat(lonlat, incrs, mapcrs)

    # Calculate map resolution if not provided
    lonlatr <- apply(lonlat, 2, range)
    if (is.null(mapres)) {
        mapres <- .lonlat_res(lonlat, lonlatr)
    }

    # Create sample map
    samplemap <- .lonlat_raster(lonlat, lonlatr, mapres, mapcrs)

    # Get cell IDs
    sm <- terra::unwrap(samplemap)
    cellid <- terra::cellFromXY(sm, lonlat)
    stopifnot(!any(is.na(cellid))) # stop if lonlat outside of raster
    # check that the cellid is the same as terra::cells
    stopifnot('SpatRaster cellid not matching' = all(sort(unique(cellid)) == terra::cells(sm)))

    # Create object using constructor
    output <- .new_marmaps(
        sample.id = sample.id,
        lonlat = lonlat,
        samplemap = samplemap,
        cellid = cellid
    )

    message("number of samples: ", length(sample.id))
    message("x range (", mapcrs, "): [", lonlatr[1, 1], ", ", lonlatr[2, 1], "]")
    message("y range (", mapcrs, "): [", lonlatr[1, 2], ", ", lonlatr[2, 2], "]")
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
#' @return A genomaps object, a list of the `geno` and `maps` inputs.
#' @seealso [margeno()] and [marmaps()] for the two components, and
#'   [MARsampling()], [MARextinction()] and [MARtheory()] for using genomaps.
#' @export
#'
#' @examples
#' # gm1001g is itself a genomaps object, rebuild it from its two components
#' geno <- gm1001g$geno
#' maps <- gm1001g$maps
#' gm <- genomaps(geno, maps)
#'
#' print(gm)
#' # genomaps() will fail if sample names are not matching
#' maps_bad <- gm1001g$maps
#' maps_bad$sample.id <- 1:1004
#' gm <- genomaps(geno, maps_bad)
genomaps <- function(geno, maps) {
    # validate inputs
    stopifnot("geno must be a margeno object" = "margeno" %in% class(geno))
    stopifnot("maps must be a marmaps object" = "marmaps" %in% class(maps))
    stopifnot("geno and maps ids must match" = all(geno$sample.id == maps$sample.id))

    # create object
    obj <- .new_genomaps(geno, maps)
    return(obj)
}
