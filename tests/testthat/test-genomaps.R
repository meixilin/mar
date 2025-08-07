# Mock data creation functions
create_mock_margeno <- function() {
    sample.id <- c("sample1", "sample2", "sample3")
    variant.id <- as.integer(1:2)
    position <- as.integer(c(100, 200))
    chromosome <- c("1", "2")
    genotype <- matrix(c(0,1,2,1,0,2), nrow=2, byrow=TRUE)
    ploidy <- 2
    return(margeno(sample.id, variant.id, position, chromosome, genotype, ploidy))
}

create_mock_marmaps <- function() {
    lonlatdf <- data.frame(
        id = c("sample1", "sample2", "sample3"),
        longitude = c(-73.935242, -118.243683, -122.419416),
        latitude = c(40.730610, 34.052235, 37.774929),
        stringsAsFactors = FALSE
    )
    return(marmaps(lonlatdf, mapres = NULL, mapcrs = "+proj=longlat +datum=WGS84"))
}

test_that("matrix validators work correctly", {
    # Test .valid_genotype
    valid_matrix <- matrix(c(0,1,2,1,0,2), nrow=2)
    expect_invisible(.valid_genotype(valid_matrix))
    invalid_matrix <- matrix(c(0,1,3,1,0,2), nrow=2)
    expect_error(.valid_genotype(invalid_matrix))

    # Test .valid_lonlat
    valid_lonlat <- matrix(c(-73.93, 40.73, -118.24, 34.05), nrow=2)
    expect_invisible(.valid_lonlat(valid_lonlat))
    invalid_lonlat <- matrix(c(-73.93, NA, -118.24, 34.05), nrow=2)
    expect_error(.valid_lonlat(invalid_lonlat))
})

test_that("margeno class works correctly", {
    mg <- create_mock_margeno()
    expect_s3_class(mg, "margeno")
    expect_equal(mg$ploidy, 2)
    expect_equal(dim(mg$genotype), c(2,3))

    # Test validation
    expect_error(margeno(
        sample.id = c('sample1', 'sample1', 'sample2'),  # Duplicate IDs
        variant.id = as.integer(1:2),
        position = as.integer(c(100, 200)),
        chromosome = c("1", "2"),
        genotype = matrix(c(0,1,2,1,0,2), nrow=2),
        ploidy = 2
    ))
})

test_that("marmaps class works correctly", {
    mm <- create_mock_marmaps()
    expect_s3_class(mm, "marmaps")
    expect_equal(length(mm$sample.id), 3)
    expect_s4_class(mm$samplemap, "RasterLayer")

    # Test invalid inputs
    expect_error(marmaps(
        data.frame(id=c(1,1), lon=c(1,2), lat=c(1,2)),  # Duplicate IDs
        mapres=NULL,
        mapcrs="+proj=longlat +datum=WGS84"
    ))

    # Test auto resolution calculation
    lonlatdf <- data.frame(
        id = c("s1", "s2"),
        longitude = c(-73.93, -118.24),
        latitude = c(40.73, 34.05)
    )
    mm2 <- marmaps(lonlatdf, mapres=NULL, mapcrs="+proj=longlat +datum=WGS84")
    expect_type(mm2$samplemap@data@values, "numeric")
})

test_that("genomaps class works correctly", {
    mg <- create_mock_margeno()
    mm <- create_mock_marmaps()

    # Test constructor
    gm <- genomaps(mg, mm)
    expect_s3_class(gm, "genomaps")
    expect_equal(gm$geno$sample.id, gm$maps$sample.id)

    # Test error on mismatched sample IDs
    mg_diff <- margeno(
        sample.id = c("sample1", "sample2", "sample4"),
        variant.id = as.integer(1:2),
        position = as.integer(c(100, 200)),
        chromosome = c("1", "2"),
        genotype = matrix(c(0,1,2,1,0,2), nrow=2, byrow=TRUE),
        ploidy = 2
    )
    expect_error(genomaps(mg_diff, mm))
})
