
# Case 1: Without metafounders

test_that("inbreeding - Case 1: Without metafounders", {

    ped <- pedigree(
        sire  = c(NA, NA, 1, 1, 4, 5),
        dam   = c(NA, NA, 2, NA, 3, 2),
        label = 1:6
    )

    FEst <- inbreeding(ped)
    FExp <- c(0.000, 0.000, 0.000, 0.000, 0.125, 0.125)

    expect_equal(FEst, FExp, tolerance = .Machine$double.eps)

})

# Case 2: Single metafounder

test_that("inbreeding - Case 2: Single metafounder", {

    ped_single <- pedigree(
        sire  = c(0, 1, 1, 2, 2, 4, 4),
        dam   = c(0, 1, 1, 3, 3, 5, 3),
        label = 1:7
    )

    gamma_single <- matrix(0.02)

    FEst <- inbreeding(ped_single, gamma = gamma_single)
    FExp <- c(-0.9800, 0.0100, 0.0100, 0.0100, 0.0100, 0.2575, 0.2575)

    expect_equal(FEst, FExp, tolerance = 1e-8)

})

# Case 3: Two metafounders

test_that("inbreeding - Case 3: Two metafounders", {

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

    FEst <- inbreeding(ped_multi, gamma = gamma_multi)
    FExp <- c(-0.990000, -0.980000, 0.005000, 0.010000, 0.025000, 0.025000, 0.264375, 0.265000)

    expect_equal(FEst, FExp, tolerance = 1e-8)
})

# Case 4: Two metafounders crosses between them

test_that("inbreeding - Case 4: Two metafounders crosses between them", {

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

    FEst <- inbreeding(ped_cross, gamma = gamma_multi)
    FExp <- c(-0.99000, -0.98000, 0.00500, 0.02500, 0.01500, 0.01500, 0.26125, 0.26375)

    expect_equal(FEst, FExp, tolerance = 1e-9)
})

# Case 6: Two metafounders crosses between metafounder and indviduals

test_that("inbreeding - Case 5: Two metafounders crosses between metafounder and indvidual", {

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

    FEst <- inbreeding(ped_cross_ind, gamma = gamma_multi)
    FExp <- c(-0.99000, -0.98000, 0.00500, 0.02500, 0.01500, 0.02500, 0.14000, 0.26375)

    expect_equal(FEst, FExp, tolerance = 1e-5)
})



