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
    sfs_result <- sfs(gm, folded = TRUE)
    expect_s3_class(sfs_result, "sfs")
    expect_true(all(sfs_result >= 0))
})

test_that("expsfs returns the neutral theta / k spectrum", {
    # 4 diploids (xN = 8), allele counts 4, 5, 3: all 3 sites segregating
    gm <- create_test_genomaps()
    theta <- 840 / 761 # 3 / H_8, with H_8 = 761 / 280

    # Test default: unfolded, k = 0 dropped
    exp_sfs <- expsfs(gm)
    expect_s3_class(exp_sfs, "sfs")
    expect_false(attr(exp_sfs, "folded"))
    expect_true(attr(exp_sfs, "nozero"))
    expect_equal(names(exp_sfs), as.character(1:8))
    expect_equal(as.vector(exp_sfs), theta / (1:8))
    expect_equal(sum(exp_sfs), 3) # matches the number of segregating sites

    # Test keeping the k = 0 bin
    exp_sfs <- expsfs(gm, nozero = FALSE)
    expect_equal(names(exp_sfs), as.character(0:8))
    expect_equal(as.vector(exp_sfs), c(0, theta / (1:8)))

    # Test folded: k and 8 - k bins collapsed, then k = 0 dropped
    exp_sfs <- expsfs(gm, folded = TRUE)
    expect_true(attr(exp_sfs, "folded"))
    expect_equal(names(exp_sfs), as.character(1:4))
    expect_equal(as.vector(exp_sfs), c(960, 560, 448, 210) / 761)
})
