test_that("harmonic number works", {
    expect_equal(Hn(1), 1)
    expect_equal(Hn(2), 1+1/2)
    expect_equal(Hn(100), sum(1/(1:100)))
    expect_equal(Hn(1000), sum(1/(1:1000)))
})


