.valid_lonlat <- function(lonlat) {
    stopifnot(is.matrix(lonlat))
    stopifnot(ncol(lonlat) == 2 & nrow(lonlat) > 0)
    stopifnot(!any(is.na(lonlat)))
    return(invisible())
}

# calculate area of a given raster
.areaofraster <- function(rr, na.rm = FALSE, tol = 1, cached = TRUE) {
    if (inherits(rr, "Raster")) {
        rr <- terra::rast(rr)
    }

    asub <- terra::cellSize(rr, unit = "km")
    # if cells' coefficient of variation too large (different area per cell)
    vals <- terra::values(asub, na.rm = TRUE)
    cv_asub <- sd(vals) / mean(vals) * 100
    # if cv > tol/100, warn the variation
    if (cv_asub > tol) {
        warning(paste("Area of raster CV =", round(cv_asub, 1), "%"))
    }
    if (cached) {
        # return the area raster
        return(asub)
    } else {
        return(terra::global(asub, fun = "sum", na.rm = TRUE)[1, 1])
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
    # create all combinations of rows and columns
    row_seq <- bbox[1]:bbox[2]
    col_seq <- bbox[3]:bbox[4]
    rc_grid <- expand.grid(row = row_seq, col = col_seq)

    # get cell indices
    cells <- terra::cellFromRowCol(mm$samplemap, rc_grid$row, rc_grid$col)
    # reverse the cells if revbbox
    if (revbbox) {
        cells <- setdiff(1:terra::ncell(mm$samplemap), cells)
    }
    cellsnotna <- intersect(mm$cellid, cells)
    return(cellsnotna)
}

.rowcol_extent <- function(mm, bbox) {
    stopifnot(length(bbox) == 4)
    # create an extent from mm$samplemap
    # When x is a Raster* object, you can pass four additional arguments to crop the
    # extent: r1, r2, c1, c2, representing the first and last row and column number
    out <- terra::ext(mm$samplemap, cells = terra::cellFromRowCol(mm$samplemap, bbox[1:2], bbox[3:4]))
    return(out)
}
