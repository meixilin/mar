# calculate area of a given raster
areaofraster <- function(rr, na.rm = FALSE, tol = 1, cached = TRUE) {
    asub <- raster::area(rr, na.rm = na.rm)
    # if cells' coefficient of variation too large (different area per cell)
    cv_asub = raster::cv(values(asub), na.rm = T)
    # if cv > tol/100, warn the variation
    if (cv_asub > tol) {
        warning(paste('Area of raster CV =', round(cv_asub, 1), '%'))
    }
    if (cached) {
        # return the area raster
        return(asub)
    } else {
        return(cellStats(asub, 'sum'))
    }
    return(asub)
}

# TODO: only works in lonlat system
areaofsquare <- function(nrow, ncol, resrow, rescol) {
    a <- nrow * ncol * resrow * rescol
    return(a)
}

# subset samples by cellids
cellid_sample <- function(gm, cellid) {
    stopifnot(length(cellid) > 0)
    sampleid <- gm$sample[gm$sample$cellid %in% cellid, 'id']
    return(sampleid)
}

# find cellids by row and column list
# bbox should be c(r1, r2, c1, c2).
# TODO: speed TBD with just using extent_sample function
rowcol_cellid <- function(gm, bbox, revbbox = FALSE) {
    stopifnot(length(bbox) == 4)
    # get the cells
    # cellFromRowColCombine returns the cell numbers obtained by the combination of all row and
    # column numbers supplied as arguments
    cells <- raster::cellFromRowColCombine(gm$samplemap, bbox[1]:bbox[2], bbox[3]:bbox[4])
    # reverse the cells if revbbox
    if (revbbox) {
        cells <- setdiff(1:ncell(gm$samplemap), cells)
    }
    cellsnotna <- intersect(gm$sample$cellid, cells)
    return(cellsnotna)
}

rowcol_extent <- function(gm, bbox) {
    stopifnot(length(bbox) == 4)
    # create an extent from gm$samplemap
    # When x is a Raster* object, you can pass four additional arguments to crop the
    # extent: r1, r2, c1, c2, representing the first and last row and column number
    out = raster::extent(gm$samplemap, bbox[1], bbox[2], bbox[3], bbox[4])
    return(out)
}

# TODO
extent_sample <- function() {
    # as polygon
    # return sampleids
}
