# for manual tests
# load('tests/testthat/testdata/gm-joshua.rda')
# scheme = 'random'
# nrep = 10
# gridded = TRUE
# myseed = 7
# load('tests/testthat/testdata/mardf-joshua.rda')
# load('tests/testthat/testdata/mares_new-joshua.rda')

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
    mardf <- mar::MARsampling(gm = gm, myseed = 7)
    # # output mardf if needed
    # save(mardf, file = 'testdata/mardf-joshua.rda')
    # fill in zero for NA
    mardf[is.na(mardf)] = 0
    # check that it completely reproduces the results from mares
    expect_equal(mardf$N, mares$N)
    expect_equal(mardf$M, mares$M)
    expect_equal(mardf$E, mares$E)
    expect_equal(mardf$thetaw, mares$thetaw)
    expect_equal(mardf$thetapi, mares$pi)
    expect_equal(mardf$A, mares$asub)
    # There is a small bug when sidesize = 28 (max of the lon/latrange)
    expect_equal(mardf$Asq, mares$a, tolerance = 0.001)
})

test_that("MARsampling random is fast", {
    load('testdata/gm-arabidopsis.rda')
    mardf = mar::MARsampling(gm = gm, myseed = 7)
    # output mardf if needed
    save(mardf, file = 'testdata/mardf-arabidopsis.rda')
    runtime = system.time(mar::MARsampling(gm = gm, myseed = 7))
    # runtime should be less than 10 seconds
    expect_lt(runtime['elapsed'], 10)
})

test_that("MARsampling set seed is consistent") {
    load('testdata/gm-joshua.rda')
    mardf0 <- mar::MARsampling(gm = gm, myseed = 7, debug = TRUE)
    mardf <- mar::MARsampling(gm = gm, myseed = 7, debug = TRUE)
    expect_equal(mardf0, mardf)
}

test_that("MARsampling random update") {
    load('testdata/gm-joshua.rda')
    mardf <- mar::MARsampling(gm = gm, myseed = 7, debug = TRUE)
    # new mardf will not have the same order compared with mares
    # new mardf are grouped by replicates to allow for inward / outward sampling
    # confirm the areas are the same
    adf <- data.frame(newa = sort(mardf$Asq), olda = mares$a)
    expect_equal(adf[1:270, 1], adf[1:270, 2])
}

test_that("MARsampling random bugfix") {
    load('testdata/gm-joshua.rda')
    mardf <- mar::MARsampling(gm = gm, myseed = 7, debug = TRUE)
    # new mardf will not have the same order compared with mares
    # new mardf are grouped by replicates to allow for inward / outward sampling
    # new mardf will have area with smaller sizes
    adf <- data.frame(newa = sort(mardf$Asq), olda = mares$a)
    # plot(adf$newa, adf$olda)
    expect_equal(adf[11:280, 1], adf[1:270, 2])
    # the largest boundary in the old version is invalid
    # plot(rowcol_extent(gm, c(1,28,9,36), add = T)
    # plot(rowcol_extent(gm, c(0,28,9,36), col = 'red', add = T)
}

test_that("Probability-based sampling") {
    latrange = 28
    # ss = 1
    ss = 5
    myprob = 1:(latrange - ss + 1)
    if (from == 'S') {myprob = rev(myprob)}
    xx = sapply(1:100000, function(i) sample(1:(latrange - ss + 1), size = 1, prob = myprob))
    hist(xx, breaks = 0:28)
}

test_that("pole_MARsampling") {
    load('testdata/gm-joshua.rda')
    latrange <- dim(gm$samplemap)[1]
    lonrange <- dim(gm$samplemap)[2]
    sidesize <- c(1,2,10)
    bblist <- .pole_MARsampling(latrange, lonrange, sidesize, from ='N')
    plot(gm$samplemap)
    sapply(bblist, function(bb) plot(rowcol_extent(gm, bb), add = T))
}

test_that('Raster xy directions') {
    plot(gm$samplemap)
    plot(rowcol_extent(gm, c(1,1,1,1)), add = T)
    plot(rowcol_extent(gm, c(2,2,1,1)), add = T, col = 'red') # latitude changes on r1,r2
    plot(rowcol_extent(gm, c(1,1,2,2)), add = T, col = 'blue') # longitude changes on c1,c2
}

