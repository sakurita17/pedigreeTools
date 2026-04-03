test_that("Mendelian sampling variance inverse - Case 1: Without metafounders", {

    # Case 1: without gamma
    ped <- pedigree(
        sire  = c(NA, NA, 1, 1, 4, 5),
        dam   = c(NA, NA, 2, NA, 3, 2),
        label = 1:6
    )

    # Lo ajuste adicionando dos decimales a los valores esperados

    # Testing Vector = TRUE
    DInvEst <- getDInv(ped, vector = TRUE)
    DInvExp <- c(1.00000000, 1.00000000, 2.00000000, 1.33333333, 2.00000000, 2.13333333)
    names(DInvExp) <- ped@label
    expect_equal(DInvEst, DInvExp, tolerance = 1e-8, ignore_attr = TRUE)

    # Testing Vector = FALSE
    DInvEst <- getDInv(ped, vector = FALSE)
    DInvExp <- Matrix::Diagonal(
        n = length(ped@label),
        x = c(1.00000000, 1.00000000, 2.00000000, 1.33333333, 2.00000000, 2.13333333))
    dimnames(DInvExp) <- list(ped@label, ped@label)
    expect_equal(DInvEst, DInvExp, tolerance = 1e-8, ignore_attr = TRUE)

})

test_that("Precision Mendelian sampling variance - Case 2: Single metafounder ", {

    ped_single <- pedigree(
        sire  = c(0, 1, 1, 2, 2, 4, 4),
        dam   = c(0, 1, 1, 3, 3, 5, 3),
        label = 1:7
    )

    gamma_single <- matrix(0.02)

    DInvEst <- getDInv(ped_single, gamma = gamma_single, vector = FALSE)
    DInvExp <- Matrix::Diagonal(
        n = length(ped_single@label),
        x = c(0, 1.010101, 1.010101, 2.020202, 2.020202, 2.020202, 2.020202))
    DInvExp[1, 1] <- solve(gamma_single)
    dimnames(DInvExp) <- list(ped_single@label, ped_single@label)

    expect_equal(DInvEst, DInvExp, tolerance = 1e-6)

})


test_that("Precision Mendelian sampling variance - Case 3: Two metafounder", {

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

    DInvEst <- getDInv(ped_multi, gamma = gamma_multi, vector = FALSE)
    DInvExp <- Matrix::Diagonal(
        n = length(ped_multi@label),
        x = c(0, 0, 1.005025, 1.010101, 2.015113, 2.015113, 2.051282, 2.035623))
    DInvExp[1:2, 1:2] <- solve(gamma_multi)
    dimnames(DInvExp) <- list(ped_multi@label, ped_multi@label)

    expect_equal(DInvEst, DInvExp, tolerance = 1e-6)

})

