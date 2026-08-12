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
