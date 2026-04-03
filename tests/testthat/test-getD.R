test_that("Mendelian sampling variance - Case 1: Without metafounders ", {

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


test_that("Mendelian sampling variance - Case 2: Single metafounder", {

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


test_that("Mendelian sampling variance - Case 3: Two metafounders", {

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

