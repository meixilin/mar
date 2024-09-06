# for manual tests
# load('tests/testthat/testdata/gm-joshua.rda')
# scheme = 'random'
# nrep = 10
# gridded = TRUE
# myseed = 7
# load('tests/testthat/testdata/mardf-joshua.rda')
# load('tests/testthat/testdata/mares_new-joshua.rda')

options(warn = -1)
library(raster)

# output mardf
test_that('MARsampling random works', {
    myspecies = c('joshua', 'arabidopsis')
    for (ii in myspecies) {
        # output mardf if needed
        if (!file.exists(paste0('testdata/mardf-', ii, '.rda'))) {
            load(paste0('testdata/gm-', ii, '.rda'))
            mardf = mar::MARsampling(gm = gm, myseed = 7)
            expect_s3_class(mardf, 'marsamp')
            save(mardf, file = paste0('testdata/mardf-', ii, '.rda'))
        }
    }
})

test_that("set seeds gives reproducible results", {
    load('testdata/gm-joshua.rda')
    for (sc in mar:::.MARsampling_schemes) {
        dflist <- lapply(1:2, function(ii) {mar::MARsampling(gm = gm, myseed = 7, scheme = sc)})
        expect_equal(dflist[[1]], dflist[[2]])
        dflist <- lapply(1:2, function(ii) {mar::MARsampling(gm = gm, myseed = 7, scheme = sc, quorum = TRUE)})
        expect_equal(dflist[[1]], dflist[[2]])
    }
})

test_that("MARsampling is fast", {
    load('testdata/gm-arabidopsis.rda')
    for (sc in mar:::.MARsampling_schemes) {
        runtime = system.time(mardf <- mar::MARsampling(gm = gm, scheme = sc, quorum = TRUE))
        # runtime should be less than 15 seconds for 10000 SNPs
        expect_lt(runtime['elapsed'], 15)
    }
})

test_that("MARsampling area bug fixed", {
    load('testdata/gm-joshua.rda')
    # confirm with the old results now obselete at v0.0.6
    load('testdata/mares_new-joshua.rda')
    mardf <- mar::MARsampling(gm = gm, myseed = 7)
    # confirm the areas are the same except that
    # new mardf will have area with smaller sizes
    # the largest boundary in the old version is invalid
    expect_equal(mardf$Asq[11:280], mares$a[1:270])
})

# # Not run but available to check
# test_that("Probability-based sampling") {
#     load('testdata/gm-joshua.rda')
#     ss = 10
#     myscheme = 'inwards'
#     mynrep = 1000
#     r0c0 = matrix(c(5,5), nrow = 1, ncol = 2)
#     bblist <- mar:::.bblist_sample(ss, gm, scheme = myscheme, nrep = mynrep, quorum = FALSE,
#                                    latrange = 28, lonrange = 45, r0c0 = r0c0, revbbox = FALSE)
#     baseraster <- gm$samplemap
#     values(baseraster) <- 0
#     plot(baseraster)
#     for(ii in seq_along(bblist)) {
#         bbextents <- rasterize(rowcol_extent(gm, bblist[[ii]]), baseraster, field = 1, background = 0)
#         baseraster <- baseraster + bbextents
#     }
#     # this plot shows the sampling frequency for all the bounding boxes created
#     # not perfectly even as the bounds are less likely to be sampled
#     plot(baseraster)
#     baseraster[5,5] <- -5
#     plot(baseraster)
# }


