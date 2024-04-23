# add tolerance = 1e-07

load('testdata/genemaps-joshua.rda')
plot_diff <- function(outdf, testdf) {
    pp = par(mfrow = c(2,3))
    plot(testdf$thetaw)
    points(outdf$thetaw, col = 'blue')
    plot(testdf$pi,ylim = c(0.1,0.4))
    points(outdf$pi, col = 'blue')
    plot(testdf$thetaw, outdf$thetaw)
    plot(testdf$pi, outdf$pi)
    plot(outdf$thetaw, outdf$pi)
    par(pp)
    return(invisible())
}

test_that("MARextinction joshua-random", {
    # run MAR extinction
    testdf = readRDS('testdata/xsim-joshua-random.rds')
    set.seed(7)
    outdf = MARextinction_sim(genemaps, scheme = "random")
    expect_equal(outdf[,3:9], testdf[,3:9], tolerance = 1e-07)
    plot_diff(outdf, testdf)
})

test_that("MARextinction joshua-inwards", {
    # run MAR extinction
    testdf = readRDS('testdata/xsim-joshua-inwards.rds')
    set.seed(7)
    outdf = MARextinction_sim(genemaps, scheme = "inwards")
    expect_equal(outdf[,3:9], testdf[,3:9], tolerance = 1e-07)
    plot_diff(outdf, testdf)
})

test_that("MARextinction joshua-southnorth", {
    # run MAR extinction
    testdf = readRDS('testdata/xsim-joshua-southnorth.rds')
    set.seed(7)
    outdf = MARextinction_sim(genemaps, scheme = "southnorth")
    expect_equal(outdf[,3:9], testdf[,3:9], tolerance = 1e-07)
    plot_diff(outdf, testdf)
})

# after commit 6289, it is not expected to PASS with old test data and new data generated
test_that("MARextinction joshua-radial", {
    # run MAR extinction
    testdf = readRDS('testdata/xsim-joshua-radial.rds')
    set.seed(7)
    outdf = MARextinction_sim(genemaps, scheme = "radial")
    expect_equal(outdf[,3:9], testdf[,3:9], tolerance = 1e-07)
    plot_diff(outdf, testdf)
})

# test the debug options
test_that("MARextinction joshua-southnorth debug", {
    # run MAR extinction
    testdf = readRDS('testdata/xsim-joshua-southnorth.rds')
    set.seed(7)
    outlist = MARextinction_sim(genemaps, scheme = "southnorth", debug = TRUE)
    outdf = outlist[[1]]
    listext = outlist[[2]]
    expect_equal(outdf[,3:9], testdf[,3:9], tolerance = 1e-07)
})
