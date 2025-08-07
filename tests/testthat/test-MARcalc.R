# Create test data
create_test_data <- function() {
    data.frame(
        A = 1:10,
        Asq = (1:10)^2,
        M = c(50, 75, 100, 125, 150, 175, 200, 225, 250, 275),
        E = c(0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.94, 0.95),
        thetaw = c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55),
        thetapi = c(12, 18, 24, 30, 36, 42, 48, 54, 60, 66)
    )
}

test_that("MARcalc basic functionality works", {
    mardf <- create_test_data()

    # Test with default parameters
    result <- MARcalc(mardf)
    expect_s3_class(result, "sars")
    expect_equal(result$model$name, "Power")

    # Test with different Mtype and Atype
    result <- MARcalc(mardf, Mtype = "E", Atype = "Asq")
    expect_s3_class(result, "sars")

    # Test handling of NA values
    mardf_na <- mardf
    mardf_na$M[1] <- NA
    result_na <- MARcalc(mardf_na)
    expect_s3_class(result_na, "sars")

    # Test handling of zero values
    mardf_zero <- mardf
    mardf_zero$M[1] <- 0
    result_zero <- MARcalc(mardf_zero)
    expect_s3_class(result_zero, "sars")
})

test_that(".marsummary works correctly", {
    mardf <- create_test_data()
    mar <- MARcalc(mardf)

    # Test successful summary
    result <- .marsummary(mar)
    expect_type(result, "list")
    expect_equal(names(result),
                c("model", "c", "z", "c_p", "z_p", "R2_adj"))
    expect_equal(result$model, "Power")
    expect_true(is.numeric(result$c))
    expect_true(is.numeric(result$z))

    # Test handling of NULL input
    null_result <- .marsummary(NULL)
    expect_true(all(is.na(unlist(null_result))))
})

test_that(".pipe_MARcalc works correctly", {
    mardf <- create_test_data()
    attr(mardf, "scheme") <- "test_scheme"

    # Test basic functionality
    result <- .pipe_MARcalc(mardf, verbose = FALSE)
    expect_s3_class(result, "marcalc")
    expect_equal(nrow(result), length(.Mtype) * length(.Atype))
    expect_equal(ncol(result), 8) # M, A, model, c, z, c_p, z_p, R2_adj

    # Check column names and types
    expect_equal(colnames(result),
                c("M", "A", "model", "c", "z", "c_p", "z_p", "R2_adj"))
    expect_true(all(result$M %in% .Mtype))
    expect_true(all(result$A %in% .Atype))
})

test_that("MARcalc handles invalid inputs correctly", {
    # Test invalid Mtype
    expect_error(MARcalc(create_test_data(), Mtype = "invalid"))

    # Test invalid Atype
    expect_error(MARcalc(create_test_data(), Atype = "invalid"))

    # Test empty data frame
    expect_warning(MARcalc(data.frame(A=numeric(), M=numeric())))

    # Test data frame with all NA values
    df_na <- data.frame(A=c(NA,NA), M=c(NA,NA))
    expect_warning(MARcalc(df_na))

    # Test data frame with all zero values
    df_zero <- data.frame(A=1:2, M=c(0,0))
    expect_warning(MARcalc(df_zero))
})

test_that("MARcalc produces expected power-law relationships", {
    mardf <- create_test_data()

    # Test with M (richness)
    mar_m <- MARcalc(mardf, Mtype = "M", Atype = "A")
    sum_m <- summary(mar_m)
    expect_true(sum_m$Parameters[2,"Estimate"] > 0) # z should be positive
    expect_true(sum_m$R2a > 0.9) # Should have good fit

    # Test with E (evenness)
    mar_e <- MARcalc(mardf, Mtype = "E", Atype = "A")
    sum_e <- summary(mar_e)
    expect_true(sum_e$Parameters[2,"Estimate"] > 0)
    expect_true(sum_e$R2a > 0.9)

    # Test with different area metrics
    mar_asq <- MARcalc(mardf, Mtype = "M", Atype = "Asq")
    sum_asq <- summary(mar_asq)
    expect_true(sum_asq$Parameters[2,"Estimate"] > 0)
    expect_true(sum_asq$R2a > 0.9)
})
