load('testdata/genemaps-joshua.rda')

test_that("MARextinction joshua-random", {
    # run MAR extinction
    testdf = readRDS('testdata/xsim-joshua-random.rds')
    set.seed(7)
    outdf = MARextinction_sim(genemaps, scheme = "random")
    expect_equal(outdf, testdf)
})

test_that("MARextinction joshua-inwards", {
    # run MAR extinction
    testdf = readRDS('testdata/xsim-joshua-inwards.rds')
    set.seed(7)
    outdf = MARextinction_sim(genemaps, scheme = "inwards")
    expect_equal(outdf, testdf)
})

test_that("MARextinction joshua-southnorth", {
    # run MAR extinction
    testdf = readRDS('testdata/xsim-joshua-southnorth.rds')
    set.seed(7)
    outdf = MARextinction_sim(genemaps, scheme = "southnorth")
    expect_equal(outdf, testdf)
})

# this is not expected to pass as the radial extinction had
test_that("MARextinction joshua-radial", {
    # run MAR extinction
    testdf = readRDS('testdata/xsim-joshua-radial.rds')
    set.seed(7)
    outdf = MARextinction_sim(genemaps, scheme = "radial")
    expect_equal(outdf, testdf)
})


