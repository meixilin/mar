# test myprob
scheme = 'inwards'
scheme = 'outwards'

rcprob <- switch(
    scheme,
    random = list(NULL, NULL),
    # no prob
    inwards = lapply(.point_prob(rvars, cvars, r0c0, ss=1), function(x) 1-x),
    outwards = .point_prob(rvars, cvars, r0c0, ss=1),
    southnorth = .pole_prob(rvars, from = 'S'),
    northsouth = .pole_prob(rvars, from = 'N')
)
myprob <- .rcprob2myprob(rcprob)
names(myprob) <- gridpresent
ss = gm$maps$samplemap
plot(ss)
ss[as.integer(names(myprob))] = myprob
plot(ss)





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


