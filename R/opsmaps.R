.valid_lonlat <- function(lonlat) {
    stopifnot("lonlat must be a matrix" = is.matrix(lonlat))
    stopifnot("lonlat must be exactly 2 columns" = ncol(lonlat) == 2)
    stopifnot("lonlat must not be empty" = nrow(lonlat) > 0)
    stopifnot("lonlat cannot have missing values" = !any(is.na(lonlat)))
    return(invisible())
}

# calculate area of a given raster
.areaofraster <- function(rr, tol = 5) {
    aa <- terra::cellSize(rr, unit = "km")
    # check if cells' coefficient of variation too large (different area per cell)
    # not using any na.rm = TRUE here because aa should not have any NA in any cases
    cv <- terra::global(aa, "sd")[1,1] / terra::global(aa, "mean")[1,1] * 100
    # if cv > tol/100, warn the variation
    if (cv > tol) {
        warning(paste("Given projection distorts area. Area coefficient of variation =", round(cv, 1), "%"))
    }
    return(aa)
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
.rowcol_cellid <- function(sm, bbox, revbbox = FALSE) {
    stopifnot(length(bbox) == 4)
    # get the cells
    cells <- terra::cellFromRowColCombine(sm, bbox[1]:bbox[2], bbox[3]:bbox[4])
    # reverse the cells if revbbox
    if (revbbox) {
        cells <- setdiff(1:terra::ncell(sm), cells)
    }
    cellsnotna <- intersect(terra::cells(sm), cells)
    return(cellsnotna)
}

