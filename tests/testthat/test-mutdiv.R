# test that mutdiv calculation is correct
# load('testdata/genemaps-joshua.rda')

library(raster)
# old mutdiv functions
test_that("old mutdiv works", {
    load('testdata/genemaps_new-joshua.rda')
    genemaps = newmaps
    raster_samples <- genemaps[[1]]
    raster_mutmaps <- genemaps[[2]]
    rest_mutmaps <- raster_mutmaps
    outmut = mutdiv_old(raster_samples, raster_mutmaps, rest_mutmaps)

    # Alternative calculations
    require(raster)
    N = cellStats(raster_samples, "sum")
    P = cellStats(raster_mutmaps, "sum")/N
    M = sum(cellStats(raster_mutmaps, "sum") > 0)
    L <- dim(raster_mutmaps)[3]
    theta <- M / (Hn(N-1) * L)
    thetapi <- (N/(N-1)) * sum(2 * P * (1 - P), na.rm = T) / L
    expect_equal(outmut$thetaw, theta)
    expect_equal(outmut$pi, thetapi)
    expect_equal(outmut$M, M)
    expect_equal(outmut$N, N)
})

# new mutdiv
# 1 warning expected due to lack of CRS in the raster
test_that("new mutdiv works", {
    load('testdata/genemaps_new-joshua.rda')
    raster_samples <- newmaps[[1]]
    raster_mutmaps <- newmaps[[2]]
    rest_mutmaps <- raster_mutmaps
    outmut = mutdiv_old(raster_samples, raster_mutmaps, rest_mutmaps)
    # load new data
    load('testdata/gm-joshua.rda')
    gmarea = areaofraster(gm$samplemap)
    newmut = .mutdiv.gridded(gm, gmarea, bbox = c(1,dim(gm$samplemap)[1],1,dim(gm$samplemap)[2]))
    expect_equal(newmut$N, outmut$N)
    expect_equal(newmut$M, outmut$M)
    expect_equal(newmut$E, 100) # outmut$E used to be 0
    expect_equal(newmut$thetaw, outmut$thetaw)
    expect_equal(newmut$thetapi, outmut$pi)
    expect_equal(newmut$A, outmut$asub)
    expect_equal(newmut$Asq, outmut$a)
})
