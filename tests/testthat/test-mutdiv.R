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

create_test_genomaps_na <- function() {
    # Create sample data
    sample.id <- c("s1", "s2", "s3", "s4")
    variant.id <- as.integer(1:3)
    position <- as.integer(c(100, 200, 300))
    chromosome <- c("1", "1", "2")
    genotype <- matrix(c(0, 1, 2, 1, 0, NA, 2, 1, 0, 1, NA, 0), nrow = 3, byrow = TRUE)
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

# pixy comparison fixtures (Korunes & Samuk 2021, Fig. 1).
# Four haploid sequences over 12 sites; "-" in the figure is missing data.
#   site:  1  2  3  4  5  6  7  8  9 10 11 12
#   h1:    -  T  T  A  -  G  G  G  C  T  A  -
#   h2:    -  G  T  A  -  G  C  G  C  -  A  -
#   h3:    -  -  T  A  -  G  G  G  C  -  A  -
#   h4:    -  G  T  -  -  G  -  G  C  T  A  -
# Sites 2 (ALT T) and 7 (ALT C) are the only polymorphic ones. Sites 1, 5 and 12
# are missing in every sequence, which is why the figure reports length 9 not 12.
# ploidy is 1 so xN is the number of called sequences and .calc_theta's
# xN * (xN - 1) denominator reproduces the figure's pairwise comparison counts.

.pixy_ref <- c("A", "G", "T", "A", "A", "G", "G", "G", "C", "T", "A", "A")
.pixy_alt <- c(".", "T", ".", ".", ".", ".", "C", ".", ".", ".", ".", ".")

# gt is a 12 x 4 character matrix of GT strings ("0", "1" or "." for missing)
.pixy_write_vcf <- function(gt) {
    header <- c(
        "##fileformat=VCFv4.2",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        paste(c(
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", "h1", "h2", "h3", "h4"
        ), collapse = "\t")
    )
    rows <- vapply(seq_len(nrow(gt)), function(i) {
        paste(c(
            "1", i, paste0("v", i), .pixy_ref[i], .pixy_alt[i], "100", "PASS",
            ".", "GT", gt[i, ]
        ), collapse = "\t")
    }, character(1))

    fn <- tempfile(fileext = ".vcf")
    writeLines(c(header, rows), fn)
    return(fn)
}

.pixy_genomaps <- function(gt) {
    vcf.fn <- .pixy_write_vcf(gt)
    on.exit(unlink(vcf.fn), add = TRUE)

    mg <- vcf_parser(vcf.fn, ploidy = 1)

    lonlatdf <- data.frame(
        id = c("h1", "h2", "h3", "h4"),
        longitude = c(-73.935, -73.934, -73.933, -73.932),
        latitude = c(40.730, 40.731, 40.732, 40.733)
    )
    mm <- marmaps(lonlatdf, mapres = 0.001, mapcrs = "+proj=longlat +datum=WGS84")

    return(genomaps(mg, mm))
}

# CASE 1 : missing data assumed invariant. Every "-" is called as
# the reference allele, so all 4 sequences are compared at all 12 sites.
create_test_pixy_genomaps_full <- function() {
    gt <- rbind(
        c("0", "0", "0", "0"), #  1  A A A A
        c("1", "0", "0", "0"), #  2  T G G G
        c("0", "0", "0", "0"), #  3  T T T T
        c("0", "0", "0", "0"), #  4  A A A A
        c("0", "0", "0", "0"), #  5  A A A A
        c("0", "0", "0", "0"), #  6  G G G G
        c("0", "1", "0", "0"), #  7  G C G G
        c("0", "0", "0", "0"), #  8  G G G G
        c("0", "0", "0", "0"), #  9  C C C C
        c("0", "0", "0", "0"), # 10  T T T T
        c("0", "0", "0", "0"), # 11  A A A A
        c("0", "0", "0", "0") # 12  A A A A
    )
    return(.pixy_genomaps(gt))
}

# CASE 2: missing data excluded from both numerator and
# denominator. Called sequences per site are 0 3 4 3 0 4 3 4 4 2 4 0, giving the
# figure's comparison counts 0 3 6 3 0 6 3 6 6 1 6 0
create_test_pixy_genomaps_na <- function() {
    gt <- rbind(
        c(".", ".", ".", "."), #  1  - - - -
        c("1", "0", ".", "0"), #  2  T G - G
        c("0", "0", "0", "0"), #  3  T T T T
        c("0", "0", "0", "."), #  4  A A A -
        c(".", ".", ".", "."), #  5  - - - -
        c("0", "0", "0", "0"), #  6  G G G G
        c("0", "1", "0", "."), #  7  G C G -
        c("0", "0", "0", "0"), #  8  G G G G
        c("0", "0", "0", "0"), #  9  C C C C
        c("0", ".", ".", "0"), # 10  T - - T
        c("0", "0", "0", "0"), # 11  A A A A
        c(".", ".", ".", ".") # 12  - - - -
    )
    return(.pixy_genomaps(gt))
}

test_that("mutdiv.gridded basic functionality works", {
    gm <- gm1001g
    gmarea <- .areaofraster(gm$maps$samplemap)

    # Test with a simple 50x50 bounding box
    result <- mutdiv.gridded(gm, gmarea, bbox = c(1, 50, 1, 50))

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
    result_rev <- mutdiv.gridded(gm, gmarea, bbox = c(1, 2, 1, 2), revbbox = TRUE)
    expect_false(identical(result, result_rev))

    # Test invalid inputs
    expect_error(mutdiv.gridded(gm, gmarea, bbox = c(1, 2, 3)))
    expect_error(mutdiv.gridded(gm, gmarea, bbox = c("a", "b", "c", "d")))
})

test_that("mutdiv.cells basic functionality works", {
    gm <- gm1001g
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
    result_subset <- .calc_theta(gm, sampleid = 1:2)
    expect_equal(result_subset$N, 2)
    expect_true(result_subset$M <= result_all$M)
})


test_that(".calc_theta correctly handles fully missing sites", {
    gm <- create_test_genomaps()
    attr(gm, "genolen") <- 1000

    gm_na <- gm
    geno <- gm_na$geno$genotype

    # make entire site missing (row 1)
    geno[1, ] <- NA
    gm_na$geno$genotype <- geno

    res <- .calc_theta(gm_na)

    # should not crash
    expect_true(is.numeric(res$thetapi))
    expect_true(is.numeric(res$thetaw))

    # segregating sites must be <= total sites - 1
    expect_true(res$M <= (nrow(geno) - 1))
})

test_that("test that .calc_theta output is correct", {
    gm <- create_test_pixy_genomaps_full()
    gm_na <- create_test_pixy_genomaps_na()

    res_full <- .calc_theta(gm)
    res_na <- .calc_theta(gm_na)
    print(res_full$M, res_full$N, res_full$E)

    expect_equal(res_full$thetapi, 0.055)
    expect_equal(res_na$thetapi, 0.099)
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

#
# test_that("diversity calculations are consistent", {
#     gm <- gm1001g
#     gmarea <- .areaofraster(gm$maps$samplemap)
#
#     # Calculate diversity using different methods
#     cellids <- unique(gm$maps$cellid)
#     bbox <- c(1, 2, 1, 2)
#
#     result_cells <- mutdiv.cells(gm, gmarea, cellids)
#     result_gridded <- mutdiv.gridded(gm, gmarea, bbox)
#     result_theta <- .calc_theta(gm)
#
#     # Check basic relationships
#     expect_true(result_cells$M <= nrow(gm$geno$genotype))
#     expect_true(result_gridded$M <= nrow(gm$geno$genotype))
#     expect_equal(result_theta$M, sum(matrixStats::rowSums2(gm$geno$genotype) > 0))
#
#     # Check that diversity measures are non-negative
#     expect_true(all(c(result_cells$thetaw, result_cells$thetapi) >= 0))
#     expect_true(all(c(result_gridded$thetaw, result_gridded$thetapi) >= 0))
#     expect_true(all(c(result_theta$thetaw, result_theta$thetapi) >= 0))
# })
