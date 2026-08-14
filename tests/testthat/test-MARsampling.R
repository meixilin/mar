test_that("MARsampling works with gm1001g genomaps data", {
    data("gm1001g", package = "mar")

    result <- mar::MARsampling(
        gm = gm1001g,
        scheme = "random",
        nrep = 2,
        quorum = FALSE,
        myseed = 123
    )

    expect_s3_class(result, "marsamp")
    expect_equal(nrow(result), 116)
    expect_identical(attr(result, "scheme"), "random")
    expect_true(all(c("N", "M", "E", "thetaw", "thetapi", "A", "Asq", "extent") %in% names(result)))
    expect_false(any(is.na(result$Asq)))
    expect_equal(
        result[14, ],
        data.frame(
            N = 3,
            M = 1016,
            E = 12,
            thetaw = 0.03525547,
            thetapi = 0.04293333,
            A = 19999.79,
            Asq = 489994.313,
            extent = "14;20;43;49"
        ),
        ignore_attr = TRUE,
        tolerance = 1e-6
    )
})

test_that("MARsampling works with animate", {
    data("gm1001g", package = "mar")

    result <- mar::MARsampling(
        gm = gm1001g,
        scheme = "random",
        nrep = 2,
        quorum = FALSE,
        animate = TRUE,
        myseed = 123
    )

    expect_s3_class(result, "marsamp")
    expect_equal(nrow(result), 116)
    expect_identical(attr(result, "scheme"), "random")
    expect_true(all(c("N", "M", "E", "thetaw", "thetapi", "A", "Asq", "extent") %in% names(result)))
    expect_false(any(is.na(result$Asq)))
    expect_equal(
        result[14, ],
        data.frame(
            N = 3,
            M = 1016,
            E = 12,
            thetaw = 0.03525547,
            thetapi = 0.04293333,
            A = 19999.79,
            Asq = 489994.313,
            extent = "14;20;43;49"
        ),
        ignore_attr = TRUE,
        tolerance = 1e-6
    )
})
