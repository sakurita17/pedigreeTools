
# Mendelian sampling variance (MS)

# Case 1: Without metafounders
test_that("case_1_msv", {

    ped <- pedigree(
        sire  = c(NA, NA, 1, 1, 4, 5),
        dam   = c(NA, NA, 2, NA, 3, 2),
        label = 1:6
    )

    DEst <- getD(ped)
    DExp <- c(1.00000, 1.00000, 0.50000, 0.75000, 0.50000, 0.46875)
    names(DExp) <- ped@label

    expect_equal(DEst, DExp, tolerance = .Machine$double.eps)
})


# Case 2: Single metafounder
test_that("case_2_msv", {

    ped_single <- pedigree(
        sire  = c(0, 1, 1, 2, 2, 4, 4),
        dam   = c(0, 1, 1, 3, 3, 5, 3),
        label = 1:7
    )

    gamma_single <- matrix(0.02)

    DEst <- getD(ped_single, gamma = gamma_single, vector = FALSE)
    DExp <- Matrix::Diagonal(
        n = length(ped_single@label),
        x = c(0.0, 0.99, 0.99, 0.495, 0.495, 0.495, 0.495))
    DExp[1, 1] <- gamma_single[1, 1]
    dimnames(DExp) <- list(ped_single@label, ped_single@label)

    expect_equal(DEst, DExp, tolerance = 1e-8)

})

# Case 3: Two metafounders
test_that("case_3_mvs", {

    ped_multi <- pedigree(
        sire  = c(0, 0, 1, 2, 3, 3, 5, 5),
        dam   = c(0, 0, 1, 2, 4, 4, 6, 4),
        label = 1:8
    )

    gamma_multi <- matrix(
        c(0.01, 0.05,
          0.05, 0.02),
        nrow = 2,
        byrow = TRUE
    )

    DEst <- getD(ped_multi, gamma = gamma_multi, vector = FALSE)
    DExp <- Matrix::Diagonal(
        n = length(ped_multi@label),
        x = c(0, 0, 0.995, 0.99, 0.49625, 0.49625, 0.4875, 0.49125))
    DExp[1:2, 1:2] <- gamma_multi
    dimnames(DExp) <- list(ped_multi@label, ped_multi@label)

    expect_equal(DEst, DExp, tolerance = 1e-8)
})

# Case 4: Two metafounders crosses between them
test_that("case_4_mvs", {

    ped_cross <- pedigree(
        sire  = c(0, 0, 1, 1, 3, 3, 5, 5),
        dam   = c(0, 0, 1, 2, 4, 4, 6, 4),
        label = 1:8
    )

    gamma_multi <- matrix(
        c(0.01, 0.05,
          0.05, 0.02),
        nrow = 2,
        byrow = TRUE
    )

    DEst <- getD(ped_cross, gamma = gamma_multi, vector = FALSE)
    DExp <- Matrix::Diagonal(
        n = length(ped_cross@label),
        x = c(0.0, 0.0, 0.995, 0.9924999999999999, 0.4925, 0.4925, 0.4925, 0.49))
    DExp[1:2, 1:2] <- gamma_multi
    dimnames(DExp) <- list(ped_cross@label, ped_cross@label)

    expect_equal(DEst, DExp, tolerance = 1e-8)
})


# Case 5: Two metafounders crosses between metafounder and indviduals
test_that("case_5_mvs", {

    ped_cross_ind <- pedigree(
        sire  = c(0, 0, 1, 1, 3, 3, 5, 5),
        dam   = c(0, 0, 1, 2, 4, 2, 6, 4),
        label = 1:8
    )

    gamma_multi <- matrix(
        c(0.01, 0.05,
          0.05, 0.02),
        nrow = 2,
        byrow = TRUE
    )

    DEst <- getD(ped_cross_ind, gamma = gamma_multi, vector = FALSE)
    DExp <- Matrix::Diagonal(
        n = length(ped_cross_ind@label),
        x = c(0.0, 0.0, 0.995, 0.9924999999999999, 0.4925, 0.74375, 0.49, 0.49))
    DExp[1:2, 1:2] <- gamma_multi
    dimnames(DExp) <- list(ped_cross_ind@label, ped_cross_ind@label)

    expect_equal(DEst, DExp, tolerance = 1e-8)
})






