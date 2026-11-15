test_that("MARtheory works with gm1001g genomaps data", {
    data("gm1001g")

    result <- MARtheory(gm = gm1001g)

    expect_s3_class(result, "martheory")
    expect_equal(nrow(result), 1004)
    expect_true(all(c("N", "M", "thetaw", "thetapi") %in% names(result)))
    expect_false(any(is.na(result)))
    expect_equal(
        result[14, ],
        data.frame(
            N = 14,
            M = 2898.0563591,
            thetaw = 0.07350813118,
            thetapi = 0.05805656199
        ),
        ignore_attr = TRUE,
        tolerance = 1e-6
    )
    # at N = J every variable site is captured
    expect_equal(
        result[1004, ],
        data.frame(
            N = 1004,
            M = 10000,
            thetaw = 0.12222158079,
            thetapi = 0.05805656199
        ),
        ignore_attr = TRUE,
        tolerance = 1e-6
    )
})

test_that("MARtheory works with gmexp genomaps data", {
    data("gmexp")

    result <- MARtheory(gm = gmexp)

    expect_s3_class(result, "martheory")
    expect_equal(nrow(result), 50)
    expect_true(all(c("N", "M", "thetaw", "thetapi") %in% names(result)))
    expect_false(any(is.na(result)))
    expect_equal(
        result[14, ],
        data.frame(
            N = 14,
            M = 399.67770724,
            thetaw = 0.11598410737,
            thetapi = 0.05737236775
        ),
        ignore_attr = TRUE,
        tolerance = 1e-6
    )
    # at N = J every variable site is captured
    expect_equal(
        result[50, ],
        data.frame(
            N = 50,
            M = 885,
            thetaw = 0.19314797822,
            thetapi = 0.05737236775
        ),
        ignore_attr = TRUE,
        tolerance = 1e-6
    )
})
