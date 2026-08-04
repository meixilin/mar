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
