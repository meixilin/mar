load('testdata/mares_new-joshua.rda')

library(raster)

test_that("set seeds works", {
    load('testdata/genemaps_new-joshua.rda')
    set.seed(7)
    lonrange <- dim(newmaps[[1]])[1]
    latrange <- dim(newmaps[[1]])[2]
    minrange <- min(lonrange, latrange)
    data <- c()
    for (sidesize in 1:minrange) {
        for (ii in 1:10) {
            xstart <- sample(1:(lonrange - sidesize), 1)
            ystart <- sample(1:(latrange - sidesize), 1)
            data <- c(data, (paste0(c(xstart, xstart + sidesize, ystart, ystart + sidesize), collapse = ';')))
        }
    }
    expect_equal(mares$extent, data)
})

test_that("MARsampling random works", {
    load('testdata/gm-joshua.rda')
    outres <- mar::MARsampling(gm = gm, myseed = 7)
    # fill in zero for NA
    outres[is.na(outres)] = 0
    # check that it completely reproduces the results from mares
    expect_equal(outres$N, mares$N)
    expect_equal(outres$M, mares$M)
    expect_equal(outres$E, mares$E)
    expect_equal(outres$thetaw, mares$thetaw)
    expect_equal(outres$thetapi, mares$pi)
    expect_equal(outres$A, mares$asub)
    # There is a small bug when sidesize = 28 (max of the lon/latrange)
    expect_equal(outres$Asq, mares$a, tolerance = 0.001)
})

test_that("MARsampling random is fast", {
    load('testdata/gm-arabidopsis.rda')
    runtime = system.time(mar::MARsampling(gm = gm, myseed = 7))
    # runtime should be less than 10 seconds
    expect_lt(runtime['elapsed'], 10)
})

