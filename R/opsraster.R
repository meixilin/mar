# calculate area of a given raster
areaofraster <- function(rr, na.rm = FALSE, tol = 1, cached = TRUE) {
    asub <- raster::area(rr, na.rm = na.rm)
    # if cells' coefficient of variation too large (different area per cell)
    cv_asub = raster::cv(values(asub), na.rm = T)
    # if cv > tol/100, stop the program
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
# TODO: speed TBD with just using extent_sample function
rowcol_cellid <- function(gm, xmin, xmax, ymin, ymax) {
    # get the cells
    cells <- raster::cellFromRowColCombine(gm$samplemap, xmin:xmax, ymin:ymax)
    cellsnotna <- intersect(gm$sample$cellid, cells)
    return(cellsnotna)
}

# TODO: orphaned function, not used anywhere yet.
rowcol_extent <- function(gm, xmin, xmax, ymin, ymax) {
    # create an extent from gm$samplemap
    out = raster::extent(gm$samplemap, xmin, xmax, ymin, ymax)
    return(out)
}

# TODO
extent_sample <- function() {
    # as polygon
    # return sampleids
}
