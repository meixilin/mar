# old areaofraster
library(raster)
areaofraster_old <- function(myraster) {
    asub <- raster::area(myraster) # better with area function
    values(asub)[is.na(values(myraster))] <- NA
    asub <- sum(values(asub), na.rm = T)
    return(asub)
}

test_that("areaofraster works", {
    load('testdata/gm-joshua.rda')
    # load('testdata/gm-arabidopsis.rda')
    oldaa = areaofraster_old(gm$samplemap)
    newaa = areaofraster(gm$samplemap, cached = FALSE, na.rm = TRUE)
    expect_equal(oldaa, newaa)
    # test that it works when not removing NAs
    aa = areaofraster(gm$samplemap, cached = FALSE)
})

# visual inspections that points fall in rasters
test_that("cellid_sample works", {
    load('testdata/gm-joshua.rda')
    # select a list of cell ids that are not empty
    # which(!is.na(values(gm$samplemap)))
    cellids = c(12, 197, 820)
    sampleids = cellid_sample(gm, cellids)
    ## plot the selected sampleids
    testmap = gm$samplemap
    testmap[setdiff(1:ncell(testmap), cellids)] <- NA
    # plot(gm$samplemap)
    # plot(testmap)
    # points(x = gm$lonlat[sampleids,1], y = gm$lonlat[sampleids,2])
})

# revbbox
test_that("revbbox") {
    load('testdata/gm-joshua.rda')
    bbox = c(1,1,1,1)
    bbox = c(13,13,1,1)
    bbox = c(1,10,1,20)
    rowcol_cellid(gm, bbox)
    ids = rowcol_cellid(gm, bbox, revbbox = T)
    testmap = gm$samplemap
    values(testmap) = NA
    testmap[ids] <- 1
    plot(gm$samplemap)
    plot(rowcol_extent(gm, bbox), add = TRUE)
    plot(testmap, add = TRUE)
}

