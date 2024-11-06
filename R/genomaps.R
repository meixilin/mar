
# matrix validators
.valid_genotype <- function(genotype) {
    # 0 = HomRef, 1 = Het, 2 = HomAlt
    valid_vars <- c(0, 1, 2)
    stopifnot(is.matrix(genotype))
    stopifnot(all(unique(as.vector(genotype)) %in% valid_vars))
    return(invisible())
}

.valid_lonlat <- function(lonlat) {
    stopifnot(is.matrix(lonlat))
    stopifnot(ncol(lonlat) == 2 & nrow(lonlat) > 0)
    stopifnot(!any(is.na(lonlat)))
    return(invisible())
}

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
    stopifnot(class(position) == "integer")
    stopifnot(class(chromosome) %in% c("character", "integer"))
    stopifnot(class(genotype) == "matrix")
    stopifnot(class(ploidy) == "numeric")

    # validate data dimensions
    stopifnot(anyDuplicated(sample.id) == 0)
    stopifnot(anyDuplicated(variant.id) == 0)
    stopifnot(length(variant.id) == dim(genotype)[1])
    stopifnot(length(sample.id) == dim(genotype)[2])

    # validate genotype
    .valid_genotype(genotype)

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
marmaps <- function(sample.id, lonlat, mapres, mapcrs) {
    # Validate inputs
    stopifnot(class(sample.id) == "character")
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

# combine the margeno/SeqArray and marmaps objects
.new_genomaps <- function(geno, maps, subs, useSeqArray) {
    obj <- list(geno = geno,
                maps = maps,
                subs = subs)
    class(obj) <- "genomaps"
    attr(obj, "useSeqArray") <- useSeqArray
    return(obj)
}

genomaps <- function(geno, maps, subs, useSeqArray) {
    # validate inputs
    stopifnot(file.exists(geno))
    stopifnot(file.exists(maps))
    stopifnot(is.list(subs))
    stopifnot(length(subs) == 2)
    stopifnot(is.logical(useSeqArray))

    # create object
    obj <- .new_genomaps(geno, maps, subs, useSeqArray)

    # test open the object and check that the sample IDs are the same
    objo <- gmOpen(obj)
    geno_sampid <- .get_genodata(objo$geno, "sample.id")
    maps_sampid <- objo$maps$sample.id
    stopifnot(all(geno_sampid == maps_sampid))
    gmClose(objo)

    # return the unopened object
    return(obj)
}

# utility function to open a genomaps object
gmOpen <- function(genomaps) {
    stopifnot(class(genomaps) == "genomaps")
    if (attr(genomaps, "useSeqArray")) {
        geno = SeqArray::seqOpen(genomaps$geno)
    } else {
        geno = get(load(genomaps$geno))
    }
    maps = get(load(genomaps$maps))

    # if subsets are needed
    if (!all(sapply(genomaps$subs, is.null))) {
        if (!is.null(genomaps$subs[[1]])) {
            geno = .filter_genosample(geno, genomaps$subs[[1]])
        }
        if (!is.null(genomaps$subs[[2]])) {
            geno = .filter_genovariant(geno, genomaps$subs[[2]])
        }
    }
    obj <- list(geno = geno,
                maps = maps)
    class(obj) <- "genomapso" # opened genomaps
    return(obj)
}

gmClose <- function(genomapso) {
    if (class(genomapso$geno) == "SeqVarGDSClass") {
        SeqArray::seqClose(genomapso$geno)
    }
    return(invisible())
}

