# test that mutdiv calculation is correct
# load('testdata/genemaps-joshua.rda')

library(raster)

# old genetic diversity estimator
mutdiv_old <- function(raster_samples, raster_mutmaps, rest_mutmaps) {
    require(raster)
    # Get the number of samples
    N <- sum(fn(values(raster_samples)), na.rm = T)
    # freqs
    P <- fn(apply(values(raster_mutmaps), 2, function(cells) sum(cells, na.rm = T))) / N
    # at least one cell has to have a presence of a mutation, and sum over
    M_ <- fn(apply(values(raster_mutmaps), 2, function(cells) any(cells > 0)))
    M_[is.na(M_)] <- 0
    M <- sum(M_)
    # find endemisms
    E_ <- apply(values(rest_mutmaps), 2, function(cells) any(cells > 0))
    E_[is.na(E_)] <- 0
    table(M_, E_)
    E <- sum(M_ & !E_)
    # Get the number of SNPs for the sample
    L <- dim(raster_mutmaps)[3]
    # compute diversity, Theta Waterson and Theta Pi (pairwise)
    if (N > 1 & M > 0) {
        theta <- M / (Hn(N-1) * L)
        thetapi <- (N/(N-1)) * sum(2 * P * (1 - P), na.rm = T) / L
    } else {
        theta <- 0
        thetapi <- 0
    }
    # area taking into account only grid cells with data
    # asub= sum(raster_samples[] > 0, na.rm = T) * (res(raster_samples)[1]*res(raster_samples)[2])
    asub <- areaofraster(raster_samples, cached = FALSE, na.rm = TRUE)

    # area based on simple square
    a <- dim(raster_samples)[1] * res(raster_samples)[1] * dim(raster_samples)[2] * res(raster_samples)[2]
    # return
    return(data.frame(thetaw = theta, pi = thetapi, M = M, E = E, N = N, a = a, asub = asub))
}

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
