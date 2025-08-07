create_test_genomaps <- function() {
    # Create sample data
    sample.id <- c("s1", "s2", "s3", "s4")
    variant.id <- as.integer(1:3)
    position <- as.integer(c(100, 200, 300))
    chromosome <- c("1", "1", "2")
    genotype <- matrix(c(0,1,2,1,0,2,2,1,0,1,2,0), nrow=3, byrow=TRUE)
    mg <- margeno(sample.id, variant.id, position, chromosome, genotype, ploidy=2)

    # Create spatial data
    lonlatdf <- data.frame(
        id = sample.id,
        longitude = c(-73.935, -73.934, -73.933, -73.932),
        latitude = c(40.730, 40.731, 40.732, 40.733)
    )
    mm <- marmaps(lonlatdf, mapres=0.001, mapcrs="+proj=longlat +datum=WGS84")

    # Create genomaps object
    gm <- genomaps(mg, mm)
    return(gm)
}

test_that("mutdiv.gridded basic functionality works", {
    gm <- create_test_genomaps()
    attr(gm, "genolen") <- 1000  # Set genome length for diversity calculations
    gmarea <- .areaofraster(gm$maps$samplemap)

    # Test with a simple 2x2 bounding box
    result <- mutdiv.gridded(gm, gmarea, bbox=c(1,2,1,2))

    # Check structure
    expect_type(result, "list")
    expect_equal(names(result), c("N", "M", "E", "thetaw", "thetapi", "A", "Asq"))

    # Check types
    expect_true(is.numeric(result$N))
    expect_true(is.numeric(result$M))
    expect_true(is.numeric(result$E))
    expect_true(is.numeric(result$thetaw))
    expect_true(is.numeric(result$thetapi))
    expect_true(is.numeric(result$A))
    expect_true(is.numeric(result$Asq))

    # Test with revbbox
    result_rev <- mutdiv.gridded(gm, gmarea, bbox=c(1,2,1,2), revbbox=TRUE)
    expect_false(identical(result, result_rev))

    # Test invalid inputs
    expect_error(mutdiv.gridded(gm, gmarea, bbox=c(1,2,3)))
    expect_error(mutdiv.gridded(gm, gmarea, bbox=c("a","b","c","d")))
})

test_that("mutdiv.cells basic functionality works", {
    gm <- create_test_genomaps()
    attr(gm, "genolen") <- 1000
    gmarea <- .areaofraster(gm$maps$samplemap)

    # Get actual cellids from the test data
    cellids <- unique(gm$maps$cellid)[1:2]

    # Test with valid cellids
    result <- mutdiv.cells(gm, gmarea, cellids)

    # Check structure
    expect_type(result, "list")
    expect_equal(names(result), c("N", "M", "E", "thetaw", "thetapi", "A", "Asq"))

    # Check that results are reasonable
    expect_true(result$N > 0)
    expect_true(result$N <= length(gm$maps$sample.id))
    expect_true(result$M <= nrow(gm$geno$genotype))

    # Test with empty cellids
    result_empty <- mutdiv.cells(gm, gmarea, character(0))
    expect_true(all(is.na(result_empty[c("N", "M", "E", "thetaw", "thetapi", "A")])))
    expect_false(is.na(result_empty$Asq))
})

test_that(".calc_theta produces valid results", {
    gm <- create_test_genomaps()
    attr(gm, "genolen") <- 1000

    # Test with all samples
    result_all <- .calc_theta(gm)

    # Check structure
    expect_type(result_all, "list")
    expect_equal(names(result_all), c("N", "M", "E", "thetaw", "thetapi"))

    # Check values are within expected ranges
    expect_equal(result_all$N, length(gm$maps$sample.id))
    expect_true(result_all$M <= nrow(gm$geno$genotype))
    expect_true(result_all$E <= result_all$M)
    expect_true(result_all$thetaw >= 0)
    expect_true(result_all$thetapi >= 0)

    # Test with subset of samples
    result_subset <- .calc_theta(gm, sampleid=1:2)
    expect_equal(result_subset$N, 2)
    expect_true(result_subset$M <= result_all$M)
})

test_that(".mutdiv.cellids handles various inputs correctly", {
    gm <- create_test_genomaps()
    attr(gm, "genolen") <- 1000
    gmarea <- .areaofraster(gm$maps$samplemap)
    Asq <- 1.0

    # Get actual cellids
    cellids <- unique(gm$maps$cellid)

    # Test with all cellids
    result <- .mutdiv.cellids(gm, gmarea, cellids, Asq)
    expect_type(result, "list")
    expect_equal(names(result), c("N", "M", "E", "thetaw", "thetapi", "A", "Asq"))

    # Test with single cellid
    result_single <- .mutdiv.cellids(gm, gmarea, cellids[1], Asq)
    expect_true(result_single$N <= result$N)

    # Test with empty cellids
    result_empty <- .mutdiv.cellids(gm, gmarea, character(0), Asq)
    expect_true(all(is.na(result_empty[c("N", "M", "E", "thetaw", "thetapi", "A")])))
    expect_equal(result_empty$Asq, Asq)
})

test_that("diversity calculations are consistent", {
    gm <- create_test_genomaps()
    attr(gm, "genolen") <- 1000
    gmarea <- .areaofraster(gm$maps$samplemap)

    # Calculate diversity using different methods
    cellids <- unique(gm$maps$cellid)
    bbox <- c(1, 2, 1, 2)

    result_cells <- mutdiv.cells(gm, gmarea, cellids)
    result_gridded <- mutdiv.gridded(gm, gmarea, bbox)
    result_theta <- .calc_theta(gm)

    # Check basic relationships
    expect_true(result_cells$M <= nrow(gm$geno$genotype))
    expect_true(result_gridded$M <= nrow(gm$geno$genotype))
    expect_equal(result_theta$M, sum(matrixStats::rowSums2(gm$geno$genotype) > 0))

    # Check that diversity measures are non-negative
    expect_true(all(c(result_cells$thetaw, result_cells$thetapi) >= 0))
    expect_true(all(c(result_gridded$thetaw, result_gridded$thetapi) >= 0))
    expect_true(all(c(result_theta$thetaw, result_theta$thetapi) >= 0))
})
