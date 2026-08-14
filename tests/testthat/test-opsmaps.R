test_that(".valid_lonlat enforces matrix structure", {
    good <- matrix(c(1, 2, 3, 4), ncol = 2)
    expect_invisible(mar:::.valid_lonlat(good))

    expect_error(mar:::.valid_lonlat(list(1, 2)))
    expect_error(mar:::.valid_lonlat(matrix(1:3, ncol = 3)))
    expect_error(mar:::.valid_lonlat(matrix(c(1, NA), ncol = 2)))
})


