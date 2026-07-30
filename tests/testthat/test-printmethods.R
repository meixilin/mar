# Helper constructors for print method tests
create_print_margeno <- function() {
    sample.id <- c("sample1", "sample2", "sample3")
    variant.id <- as.integer(1:2)
    position <- as.integer(c(100, 200))
    chromosome <- c("1", "2")
    genotype <- matrix(c(0, 1, 2, 1, 0, 2), nrow = 2, byrow = TRUE)
    ploidy <- 2
    margeno(genotype, ploidy, sample.id, variant.id, position, chromosome)
}

create_print_marmaps <- function() {
    lonlatdf <- data.frame(
        id = c("sample1", "sample2", "sample3"),
        longitude = c(-73.935242, -118.243683, -122.419416),
        latitude = c(40.730610, 34.052235, 37.774929),
        stringsAsFactors = FALSE
    )
    marmaps(lonlatdf, mapres = NULL, mapcrs = "+proj=longlat +datum=WGS84")
}

test_that("print.margeno summarizes geno data", {
    mg <- create_print_margeno()

    expect_output(print(mg), "margeno object")
    expect_output(print(mg), "number of samples:  3")
    expect_output(print(mg), "number of genomic sites:  2")
    expect_output(print(mg), "ploidy:  2")
    expect_output(print(mg), "head of variantid, position, chromosome, genotype:")
    expect_output(print(mg), "sample1")
    expect_invisible(print(mg))
})

test_that("print.marmaps summarizes map data", {
    mm <- create_print_marmaps()

    expect_output(print(mm), "marmaps object")
    expect_output(print(mm), "number of samples:  3")
    expect_output(print(mm), "longitude range:")
    expect_output(print(mm), "latitude range:")
    expect_output(print(mm), "samplemap raster layer:")
    expect_output(print(mm), "SpatRaster")
    expect_output(print(mm), "head of sampleid and lonlat:")
    expect_output(print(mm), "sample1")
    expect_invisible(print(mm))
})
