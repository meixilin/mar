# Create test genomaps data
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

test_that("MARsad basic functionality works", {
    gm <- create_test_genomaps()

    # Test with default parameters
    result <- MARsad(gm)
    expect_s3_class(result, "marsad")
    expect_type(result$sadms, "list")
    expect_true(all(.sad_models %in% names(result$sadms)))

    # Test with predict=FALSE
    result_no_pred <- MARsad(gm, predict=FALSE)
    expect_null(result_no_pred$sadsfss)

    # Test with subset of models
    result_subset <- MARsad(gm, sad_models=c("bs", "lnorm"))
    expect_equal(length(result_subset$sadms), 2)
})

test_that(".sadpred works correctly", {
    gm <- create_test_genomaps()
    sad_fit <- MARsad(gm)

    # Test prediction for each model
    for(model in names(sad_fit$sadms)) {
        pred <- .sadpred(sad_fit$sadms[[model]],
                        N=length(gm$maps$sample.id),
                        ploidy=gm$geno$ploidy,
                        folded=TRUE)
        expect_s3_class(pred, "sfs")
        expect_true(all(pred >= 0))  # SFS should be non-negative
        expect_true(sum(pred) > 0)   # SFS should not be all zeros
    }
})

test_that("sfs functions work correctly", {
    # Test .foldsfs
    unfold <- c(1,2,3,4,3,2,1)
    fold <- .foldsfs(unfold)
    expect_equal(fold, c(2,4,6,4))

    # Test .new_sfs
    sfs_obj <- .new_sfs(c(0,1,2,3,2,1), folded=TRUE, nozero=TRUE)
    expect_s3_class(sfs_obj, "sfs")
    expect_true(attr(sfs_obj, "folded"))
    expect_true(attr(sfs_obj, "nozero"))

    # Test sfs function using genomaps data
    gm <- create_test_genomaps()
    AC <- .get_AC(gm$geno)
    sfs_result <- sfs(AC, N=length(gm$maps$sample.id), ploidy=gm$geno$ploidy, folded=TRUE)
    expect_s3_class(sfs_result, "sfs")
    expect_true(all(sfs_result >= 0))
})

test_that(".get_AC works correctly", {
    gm <- create_test_genomaps()

    # Test normal case
    AC <- .get_AC(gm$geno)
    expect_type(AC, "double")
    expect_true(all(AC >= 0))

    # Test with zero counts
    gm_zero <- create_test_genomaps()
    gm_zero$geno$genotype[,1] <- 0
    expect_warning(.get_AC(gm_zero$geno))
})

test_that("ll_sfs works correctly", {
    gm <- create_test_genomaps()
    AC <- .get_AC(gm$geno)
    N <- length(gm$maps$sample.id)
    ploidy <- gm$geno$ploidy

    # Create actual data SFS
    data_sfs <- sfs(AC, N, ploidy, folded=TRUE)

    # Create expected SFS
    model_sfs <- expsfs(lenAC=length(AC), N=N, ploidy=ploidy, folded=TRUE)

    # Test basic functionality
    ll <- ll_sfs(model_sfs, data_sfs)
    expect_type(ll, "double")

    # Test warning for zero model values
    model_zero <- model_sfs
    model_zero[1] <- 0
    expect_warning(ll_sfs(model_zero, data_sfs))

    # Test error for mismatched lengths
    model_short <- model_sfs[-1]
    expect_error(ll_sfs(model_short, data_sfs))
})

test_that(".pipe_sadsfs works correctly", {
    gm <- create_test_genomaps()
    marsad <- MARsad(gm)
    genosfs <- sfs(.get_AC(gm$geno),
                   N=length(gm$maps$sample.id),
                   ploidy=gm$geno$ploidy,
                   folded=TRUE)

    result <- .pipe_sadsfs(gm, marsad, genosfs, folded=TRUE)

    expect_type(result, "list")
    expect_equal(names(result), c("sfs", "statdf"))
    expect_true(all(c("data", "neutral") %in% names(result$sfs)))
    expect_equal(colnames(result$statdf), c("model", "logLik"))
})
