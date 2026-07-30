# Create test genomaps data
create_test_genomaps <- function() {
    # Create sample data
    sample.id <- c("s1", "s2", "s3", "s4")
    variant.id <- as.integer(1:3)
    position <- as.integer(c(100, 200, 300))
    chromosome <- c("1", "1", "2")
    genotype <- matrix(c(0, 1, 2, 1, 0, 2, 2, 1, 0, 1, 2, 0), nrow = 3, byrow = TRUE)
    mg <- margeno(genotype, ploidy = 2, sample.id, variant.id, position, chromosome)

    # Create spatial data
    lonlatdf <- data.frame(
        id = sample.id,
        longitude = c(-73.935, -73.934, -73.933, -73.932),
        latitude = c(40.730, 40.731, 40.732, 40.733)
    )
    mm <- marmaps(lonlatdf, mapres = 0.001, mapcrs = "+proj=longlat +datum=WGS84")

    # Create genomaps object
    gm <- genomaps(mg, mm)
    return(gm)
}


test_that("sfs functions work correctly", {
    # Test .foldsfs
    unfold <- c(1, 2, 3, 4, 3, 2, 1)
    fold <- .foldsfs(unfold)
    expect_equal(fold, c(2, 4, 6, 4))

    # Test .new_sfs
    sfs_obj <- .new_sfs(c(0, 1, 2, 3, 2, 1), folded = TRUE, nozero = TRUE)
    expect_s3_class(sfs_obj, "sfs")
    expect_true(attr(sfs_obj, "folded"))
    expect_true(attr(sfs_obj, "nozero"))

    # Test sfs function using genomaps data
    gm <- create_test_genomaps()
    AC <- gm$geno$allel_count
    sfs_result <- sfs(gm, folded = TRUE)
    expect_s3_class(sfs_result, "sfs")
    expect_true(all(sfs_result >= 0))
})
