# test that mutdiv calculation is correct
# load('testdata/genemaps-joshua.rda')

# create new genemaps
load('testdata/genemaps_new-joshua.rda')
genemaps = newmaps

test_that("mutdiv works", {
    raster_samples <- genemaps[[1]]
    raster_mutmaps <- genemaps[[2]]
    rest_mutmaps <- raster_mutmaps
    outmut = mutdiv(raster_samples, raster_mutmaps, rest_mutmaps)

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
