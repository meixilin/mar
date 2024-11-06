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
    sample.id <- c("sample1", "sample2", "sample3")
    lonlat <- matrix(c(-73.935242, 40.730610,  # New York
                      -118.243683, 34.052235,  # Los Angeles
                      -122.419416, 37.774929), # San Francisco
                     nrow=3, byrow=TRUE)
    colnames(lonlat) <- c("LONGITUDE", "LATITUDE")
    return(marmaps(sample.id, lonlat, mapres = NULL, mapcrs = "+proj=longlat +datum=WGS84"))
}

# Test cases
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
    # Test constructor
    mg <- create_mock_margeno()
    expect_s3_class(mg, "margeno")
    expect_equal(mg$ploidy, 2)
    expect_equal(dim(mg$genotype), c(2,3))

    # Test validation
    expect_error(margeno(
        sample.id = c('sample1', 'sample1', 'sample2'),
        variant.id = as.integer(1:2),
        position = as.integer(c(100, 200)),
        chromosome = c("1", "2"),
        genotype = matrix(c(0,1,2,1,0,2), nrow=2),
        ploidy = 2
    ))

    expect_error(margeno(
        sample.id = c('sample1', 'sample3', 'sample2'),
        variant.id = as.integer(1:2),
        position = as.integer(c(100, 200)),
        chromosome = c("1", "2"),
        genotype = matrix(c(0,1,2,1,0,NA), nrow=2),
        ploidy = 2
    ))
})

test_that("marmaps class works correctly", {
    # Test constructor
    mm <- create_mock_marmaps()
    expect_s3_class(mm, "marmaps")
    expect_equal(length(mm$sample.id), 3)
    expect_s4_class(mm$samplemap, "RasterLayer")

    # Test plotting
    expect_invisible(plot(mm))

    # Test resolution calculation
    lonlat <- matrix(c(-73.93, 40.73, -118.24, 34.05), nrow=2)
    lonlatr <- apply(lonlat, 2, range)
    res <- .lonlat_res(lonlat, lonlatr)
    expect_type(res, "double")
})

# Create temporary files for testing
temp_geno <- tempfile(fileext = ".rda")
temp_maps <- tempfile(fileext = ".rda")
temp_genomaps <- tempfile(fileext = ".rda")
temp_genomaps2 <- tempfile(fileext = ".rda")

test_that("genomaps class works correctly", {
    # Save mock objects
    mg <- create_mock_margeno()
    mm <- create_mock_marmaps()
    save(mg, file = temp_geno)
    save(mm, file = temp_maps)

    # Test constructor
    gm <- genomaps(temp_geno, temp_maps,
                   subs=list(NULL, NULL),
                   useSeqArray=FALSE)
    expect_s3_class(gm, "genomaps")
    expect_false(attr(gm, "useSeqArray"))
    save(gm, file = temp_genomaps)
})

test_that("genomaps class works correctly for SeqArray", {
    # load mock margeno
    mg <- get(load(temp_geno))
    # convert to SeqVarGDSClass
    gds <- .margeno2seqarray(mg)
    gm2 <- genomaps(gds, temp_maps,
                    subs=list(NULL, NULL),
                    useSeqArray=TRUE)
    expect_s3_class(gm2, "genomaps")
    expect_true(attr(gm2, "useSeqArray"))
    save(gm2, file = temp_genomaps2)

    # Clean up GDS file
    seqClose(gds)
    unlink(gds@filename)
})

test_that("genomaps operations work correctly", {
    gm <- get(load(temp_genomaps))
    gmo <- gmOpen(gm)
    expect_s3_class(gmo, "genomapso")
    expect_equal(names(gmo), c("geno", "maps"))
    expect_equal(class(gmo$geno), "margeno")
    expect_equal(class(gmo$maps), "marmaps")
    expect_invisible(gmClose(gmo))

    gm2 <- get(load(temp_genomaps2))
    gm2o <- gmOpen(gm2)
    expect_s3_class(gm2o, "genomapso")
    expect_equal(names(gm2o), c("geno", "maps"))
    expect_equal(class(gm2o$geno), "SeqVarGDSClass", check.attributes = FALSE)
    expect_equal(class(gm2o$maps), "marmaps")
    expect_invisible(gmClose(gm2o))
})

test_that("genomaps throws error when sample IDs don't match", {
    # Create margeno with different sample IDs
    mg_diff <- margeno(
        sample.id = c("sample1", "sample2", "sample4"),  # Different from marmaps
        variant.id = as.integer(1:2),
        position = as.integer(c(100, 200)),
        chromosome = c("1", "2"),
        genotype = matrix(c(0,1,2,1,0,2), nrow=2, byrow=TRUE),
        ploidy = 2
    )

    # Create marmaps with original sample IDs
    mm <- create_mock_marmaps()

    # Save to temporary files
    temp_geno_diff <- tempfile(fileext = ".rda")
    temp_maps_diff <- tempfile(fileext = ".rda")
    save(mg_diff, file = temp_geno_diff)
    save(mm, file = temp_maps_diff)

    # Test that genomaps throws error due to mismatched sample IDs
    expect_error(
        genomaps(temp_geno_diff, temp_maps_diff,
                subs=list(NULL, NULL),
                useSeqArray=FALSE),
        "all.*==.*"  # Error message contains the stopifnot comparison
    )

    # Clean up
    unlink(temp_geno_diff)
    unlink(temp_maps_diff)
})

test_that("genomaps filtering works correctly", {
    # Create mock data with more samples and variants
    sample.id <- c("sample1", "sample2", "sample3", "sample4", "sample5")
    variant.id <- as.integer(1:4)
    position <- as.integer(c(100, 200, 300, 400))
    chromosome <- c("1", "1", "2", "2")
    genotype <- matrix(c(
        0,1,2,1,0,
        1,0,2,0,1,
        2,1,0,1,2,
        0,2,1,0,1
    ), nrow=4, byrow=TRUE)
    ploidy <- 2
    mg <- margeno(sample.id, variant.id, position, chromosome, genotype, ploidy)

    # Create corresponding marmaps
    lonlat <- matrix(c(
        -73.935242, 40.730610,  # New York
        -118.243683, 34.052235, # Los Angeles
        -122.419416, 37.774929, # San Francisco
        -87.632401, 41.883228,  # Chicago
        -71.058880, 42.360082   # Boston
    ), nrow=5, byrow=TRUE)
    colnames(lonlat) <- c("LONGITUDE", "LATITUDE")
    mm <- marmaps(sample.id, lonlat, mapres = NULL, mapcrs = "+proj=longlat +datum=WGS84")

    # Save to temporary files
    temp_geno <- tempfile(fileext = ".rda")
    temp_maps <- tempfile(fileext = ".rda")
    save(mg, file = temp_geno)
    save(mm, file = temp_maps)

    # Test sample filtering
    sample_subset <- c("sample1", "sample3", "sample5")
    gm_filtered_samples <- genomaps(temp_geno, temp_maps,
                                  subs=list(sample_subset, NULL),
                                  useSeqArray=FALSE)
    gmo <- gmOpen(gm_filtered_samples)
    expect_equal(gmo$geno$sample.id, sample_subset)
    expect_equal(dim(gmo$geno$genotype)[2], length(sample_subset))
    gmClose(gmo)

    # Test variant filtering
    variant_subset <- c(1L, 4L)
    gm_filtered_variants <- genomaps(temp_geno, temp_maps,
                                   subs=list(NULL, variant_subset),
                                   useSeqArray=FALSE)
    gmo <- gmOpen(gm_filtered_variants)
    expect_equal(gmo$geno$variant.id, variant_subset)
    expect_equal(dim(gmo$geno$genotype)[1], length(variant_subset))
    gmClose(gmo)

    # Test both sample and variant filtering
    gm_filtered_both <- genomaps(temp_geno, temp_maps,
                                subs=list(sample_subset, variant_subset),
                                useSeqArray=FALSE)
    gmo <- gmOpen(gm_filtered_both)
    expect_equal(gmo$geno$sample.id, sample_subset)
    expect_equal(gmo$geno$variant.id, variant_subset)
    expect_equal(dim(gmo$geno$genotype), c(length(variant_subset), length(sample_subset)))
    gmClose(gmo)

    # Clean up
    unlink(c(temp_geno, temp_maps))
})

# Clean up
unlink(c(temp_geno, temp_maps, temp_genomaps, temp_genomaps2))
closeAllConnections()
