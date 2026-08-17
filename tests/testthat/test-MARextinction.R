test_that(".rcprob2myprob works correctly", {
    gridpresent <- 1:5

    # Test with both probabilities NULL
    rcprob <- list(NULL, NULL)
    result <- .rcprob2myprob(rcprob, gridpresent)
    expect_null(result)

    # Test with one probability NULL
    rcprob <- list(NULL, c(0.2, 0.2, 0.2, 0.2, 0.2))
    result <- .rcprob2myprob(rcprob, gridpresent)
    expect_equal(length(result), 5)
    expect_equal(sum(result), 1)
    expect_equal(names(result), as.character(gridpresent))

    # Test with both probabilities present
    rcprob <- list(c(0.2, 0.2, 0.2, 0.2, 0.2), c(0.2, 0.2, 0.2, 0.2, 0.2))
    result <- .rcprob2myprob(rcprob, gridpresent)
    expect_equal(length(result), 5)
    expect_equal(sum(result), 1)
})

test_that(".rescale_prob works correctly", {
    # Test with non-NULL input
    prob <- c(1, 2, 3, 4)
    result <- .rescale_prob(prob)
    expect_equal(sum(result), 1)
    expect_equal(length(result), 4)

    # Test with NULL input
    expect_null(.rescale_prob(NULL))
})

test_that(".extsample works correctly", {
    gridpresent <- 1:10
    myprob <- rep(0.1, 10)
    names(myprob) <- gridpresent

    # Test with default step size
    result <- .extsample(gridpresent, myprob, mystep = 2)
    expect_type(result, "list")
    expect_equal(length(result[[1]]), 10) # First step has all cells
    expect_true(length(result[[length(result)]]) <= 2) # Last step has <= mystep cells

    # Test with NULL probability (random sampling)
    result <- .extsample(gridpresent, NULL, mystep = 2)
    expect_type(result, "list")
    expect_equal(length(result[[1]]), 10)
})

test_that("MARextinction basic functionality works", {
    data("gm1001g", package = "mar")

    # Test with default parameters
    result <- MARextinction(gm1001g, nrep = 2)
    expect_s3_class(result, "marextinct")
    expect_equal(attr(result, "scheme"), "random")

    # Test with different schemes
    schemes <- c("random", "inwards", "outwards", "southnorth", "northsouth")
    for (scheme in schemes) {
        result <- MARextinction(gm1001g, scheme = scheme, nrep = 2)
        expect_s3_class(result, "marextinct")
        expect_equal(attr(result, "scheme"), scheme)
    }

    # Test with seed
    result1 <- MARextinction(gm1001g, myseed = 123, nrep = 2)
    result2 <- MARextinction(gm1001g, myseed = 123, nrep = 2)
    expect_equal(result1, result2)
    expect_equal(nrow(result1), 168)
    expect_true(all(c("N", "M", "E", "thetaw", "thetapi", "A", "extl", "repid") %in% names(result1)))
    expect_false(any(is.na(result1)))
    expect_equal(
        result1[83, ],
        data.frame(
            N = 12,
            M = 1634,
            E = 15,
            thetaw = 0.04102518,
            thetapi = 0.04937101,
            A = 29999.67,
            extent = "1415;1539;2019",
            repid = 1
        ),
        ignore_attr = TRUE,
        tolerance = 1e-6
    )
})

test_that("MARextinction handles edge cases correctly", {
    gm <- create_test_genomaps()

    # Test with single replicate
    result <- MARextinction(gm, nrep = 1)
    expect_s3_class(result, "marextinct")
    expect_true(all(result$repid == 1))

    # Test with very small xfrac
    result <- MARextinction(gm, xfrac = 0.001, nrep = 2)
    expect_s3_class(result, "marextinct")

    # Test with very large xfrac
    result <- MARextinction(gm, xfrac = 0.5, nrep = 2)
    expect_s3_class(result, "marextinct")
})

test_that("MARextinction maintains data integrity", {
    gm <- create_test_genomaps()
    result <- MARextinction(gm, nrep = 2)

    # Check data frame structure
    expect_true(all(c("A", "M", "E", "thetaw", "thetapi", "extl", "repid") %in% colnames(result)))

    # Check numeric columns are numeric
    numeric_cols <- c("A", "M", "E", "thetaw", "thetapi")
    for (col in numeric_cols) {
        expect_type(result[[col]], "double")
    }

    # Check extinction list column format
    expect_type(result$extl, "character")

    # Check replicate ID format
    expect_type(result$repid, "integer")
    expect_true(all(result$repid %in% 1:2))
})

test_that(".animate_MARextinction works", {
    skip_if_not(interactive(), "Animation tests only run in interactive mode")

    gm <- create_test_genomaps()
    extlist <- .extlist_sample(gm, xfrac = 0.1, scheme = "random", nrep = 1, r0c0 = c(1,1))

    # Test animation function doesn't error
    expect_invisible(.animate_MARextinction(gm, extlist))
})
