.valid_lonlat <- function(lonlat) {
    stopifnot(is.matrix(lonlat))
    stopifnot(ncol(lonlat) == 2 & nrow(lonlat) > 0)
    stopifnot(!any(is.na(lonlat)))
    return(invisible())
}

# calculate area of a given raster
.areaofraster <- function(rr, na.rm = FALSE, tol = 1, cached = TRUE) {
    asub <- raster::area(rr, na.rm = na.rm)
    # if cells' coefficient of variation too large (different area per cell)
    cv_asub = raster::cv(raster::values(asub), na.rm = T)
    # if cv > tol/100, warn the variation
    if (cv_asub > tol) {
        warning(paste('Area of raster CV =', round(cv_asub, 1), '%'))
    }
    if (cached) {
        # return the area raster
        return(asub)
    } else {
        return(raster::cellStats(asub, 'sum'))
    }
    return(asub)
}

# TODO: only works in lonlat system
.areaofsquare <- function(nrow, ncol, resrow, rescol) {
    a <- nrow * ncol * resrow * rescol
    return(a)
}

# subset samples by cellids
.cellid_sample <- function(mm, cellid) {
    stopifnot(length(cellid) > 0)
    # use numeric indexing as genotype matrix has no dimnames
    sampleid <- which(mm$cellid %in% cellid)
    stopifnot(length(sampleid) > 0)
    return(sampleid)
}

# find cellids by row and column list
# bbox should be c(r1, r2, c1, c2).
# TODO: speed TBD with just using extent_sample function
.rowcol_cellid <- function(mm, bbox, revbbox = FALSE) {
    stopifnot(length(bbox) == 4)
    # get the cells
    # cellFromRowColCombine returns the cell numbers obtained by the combination of all row and
    # column numbers supplied as arguments
    cells <- raster::cellFromRowColCombine(mm$samplemap, bbox[1]:bbox[2], bbox[3]:bbox[4])
    # reverse the cells if revbbox
    if (revbbox) {
        cells <- setdiff(1:raster::ncell(mm$samplemap), cells)
    }
    cellsnotna <- intersect(mm$cellid, cells)
    return(cellsnotna)
}

.rowcol_extent <- function(mm, bbox) {
    stopifnot(length(bbox) == 4)
    # create an extent from mm$samplemap
    # When x is a Raster* object, you can pass four additional arguments to crop the
    # extent: r1, r2, c1, c2, representing the first and last row and column number
    out = raster::extent(mm$samplemap, bbox[1], bbox[2], bbox[3], bbox[4])
    return(out)
}
