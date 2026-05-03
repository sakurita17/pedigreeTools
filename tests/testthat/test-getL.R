test_that("getLInv - Case 1: Without metafounders", {
  
  ped <- pedigree(
    sire  = c(NA, NA, 1, 1, 4, 5),
    dam   = c(NA, NA, 2, NA, 3, 2),
    label = 1:6
  )
  
  LEst <- getL(ped)
  LExp <- matrix(data = c(
    1.0000, 0.0000, 0.5000, 0.5000, 0.5000, 0.2500,
    0.0000, 1.0000, 0.5000, 0.0000, 0.2500, 0.6250,
    0.0000, 0.0000, 0.7071, 0.0000, 0.3536, 0.1768,
    0.0000, 0.0000, 0.0000, 0.8660, 0.4330, 0.2165,
    0.0000, 0.0000, 0.0000, 0.0000, 0.7071, 0.3536,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.6847),
    byrow = TRUE, nrow = 6)

  expect_equal(unname(as.matrix(LEst)), unname(as.matrix(LExp)), tolerance = 1e-4)
  expect_equal(unname(as.matrix(LEst)), unname(as.matrix(chol(getA(ped)))), tolerance = 1e-4)
})


test_that("getL - Case 2: Single metafounder", {

    ped_single <- pedigree(
        sire  = c(0, 1, 1, 2, 2, 4, 4),
        dam   = c(0, 1, 1, 3, 3, 5, 3),
        label = 1:7
    )

    gamma_single <- matrix(0.02)

    LEst <- getL(ped_single, gamma = gamma_single)
    A <- matrix(c(
        0.02000, 0.020000, 0.02000, 0.02000, 0.02000, 0.02000, 0.02000,
        0.02000, 1.01000, 0.02000, 0.51500, 0.515000, 0.51500, 0.26750,
        0.02000, 0.02000, 1.01000, 0.515000, 0.51500, 0.51500, 0.76250,
        0.02000, 0.51500, 0.51500, 1.01000, 0.51500, 0.76250, 0.76250,
        0.02000, 0.51500, 0.51500, 0.51500, 1.01000, 0.76250, 0.51500,
        0.02000, 0.51500, 0.51500, 0.76250, 0.76250, 1.25750, 0.63875,
        0.02000, 0.26750, 0.76250, 0.76250, 0.51500, 0.63875, 1.25750
    ), nrow = 7, byrow = TRUE)
    
    dimnames(A) <- list(ped_single@label, ped_single@label)
    expect_equal(unname(as.matrix(crossprod(LEst))), unname(as.matrix(A)), tolerance = 1e-6)
})

test_that("getL - Case 3: Two metafounders", {
    
    ped_multi <- pedigree(
        sire  = c(0, 0, 1, 2, 3, 3, 5, 5),
        dam   = c(0, 0, 1, 2, 4, 4, 6, 4),
        label = 1:8
    )

    gamma_multi <- matrix(
        c(0.1, 0.05,
          0.05, 0.2),
        nrow = 2,
        byrow = TRUE
    )

    LEst <- getL(ped_multi, gamma = gamma_multi)
    A <- matrix(c(
        0.100000, 0.050000, 0.100000, 0.050000, 0.075000, 0.075000, 0.075000, 0.062500,
        0.050000, 0.200000, 0.050000, 0.200000, 0.125000, 0.125000, 0.125000, 0.162500,
        0.100000, 0.050000, 1.050000, 0.050000, 0.550000, 0.550000, 0.550000, 0.300000,
        0.050000, 0.200000, 0.050000, 1.100000, 0.575000, 0.575000, 0.575000, 0.837500,
        0.075000, 0.125000, 0.550000, 0.575000, 1.025000, 0.562500, 0.793750, 0.800000,
        0.075000, 0.125000, 0.550000, 0.575000, 0.562500, 1.025000, 0.793750, 0.568750,
        0.075000, 0.125000, 0.550000, 0.575000, 0.793750, 0.793750, 1.281250, 0.684375,
        0.062500, 0.162500, 0.300000, 0.837500, 0.800000, 0.568750, 0.684375, 1.287500),
        nrow = 8,
        byrow = TRUE
    )

    dimnames(A) <- list(ped_multi@label, ped_multi@label)
    expect_equal(unname(as.matrix(crossprod(LEst))), unname(as.matrix(A)), tolerance = 1e-6)
})

test_that("getL - Case 4: Two metafounders crosses between them", {

    ped_cross <- pedigree(
        sire  = c(0, 0, 1, 1, 3, 3, 5, 5),
        dam   = c(0, 0, 1, 2, 4, 4, 6, 4),
        label = 1:8
    )

    gamma_multi <- matrix(
        c(0.1, 0.05,
          0.05, 0.2),
        nrow = 2,
        byrow = TRUE
    )

    LEst <- getL(ped_cross, gamma = gamma_multi)
    A <- Matrix::Matrix(c(
        0.100000, 0.050000, 0.100000, 0.075000, 0.087500, 0.087500, 0.087500, 0.081250,
        0.050000, 0.200000, 0.050000, 0.125000, 0.087500, 0.087500, 0.087500, 0.106250,
        0.100000, 0.050000, 1.050000, 0.075000, 0.562500, 0.562500, 0.562500, 0.318750,
        0.075000, 0.125000, 0.075000, 1.025000, 0.550000, 0.550000, 0.550000, 0.787500,
        0.087500, 0.087500, 0.562500, 0.550000, 1.037500, 0.556250, 0.796875, 0.793750,
        0.087500, 0.087500, 0.562500, 0.550000, 0.556250, 1.037500, 0.796875, 0.553125,
        0.087500, 0.087500, 0.562500, 0.550000, 0.796875, 0.796875, 1.278120, 0.673437,
        0.081250, 0.106250, 0.318750, 0.787500, 0.793750, 0.553125, 0.673437, 1.275000
    ), nrow = 8, byrow = TRUE)

    dimnames(A) <- list(ped_cross@label, ped_cross@label)
    expect_equal(unname(as.matrix(crossprod(LEst))), unname(as.matrix(A)), tolerance = 1e-6)
})

test_that("getL - Case 5: Two metafounders crosses between metafounder and individual", {

    ped_cross_ind <- pedigree(
        sire  = c(0, 0, 1, 1, 3, 3, 5, 5),
        dam   = c(0, 0, 1, 2, 4, 2, 6, 4),
        label = 1:8
    )

    gamma_multi <- matrix(
        c(0.1, 0.05,
          0.05, 0.2),
        nrow = 2,
        byrow = TRUE
    )

    LEst <- getL(ped_cross_ind, gamma = gamma_multi)
    A <- matrix(c(
        0.100000, 0.050000, 0.100000, 0.075000, 0.087500, 0.075000, 0.081250, 0.081250,
        0.050000, 0.200000, 0.050000, 0.125000, 0.087500, 0.125000, 0.106250, 0.106250,
        0.100000, 0.050000, 1.050000, 0.075000, 0.562500, 0.550000, 0.556250, 0.318750,
        0.075000, 0.125000, 0.075000, 1.025000, 0.550000, 0.100000, 0.325000, 0.787500,
        0.087500, 0.087500, 0.562500, 0.550000, 1.037500, 0.325000, 0.681250, 0.793750,
        0.075000, 0.125000, 0.550000, 0.100000, 0.325000, 1.025000, 0.675000, 0.212500,
        0.081250, 0.106250, 0.556250, 0.325000, 0.681250, 0.675000, 1.162500, 0.503125,
        0.081250, 0.106250, 0.318750, 0.787500, 0.793750, 0.212500, 0.503125, 1.275000),
        nrow = 8,
        byrow = TRUE
    )

    dimnames(A) <- list(ped_cross_ind@label, ped_cross_ind@label)
    expect_equal(unname(as.matrix(crossprod(LEst))), unname(as.matrix(A)), tolerance = 1e-6)
})

# TODO: Test if getL with gamma = chol(A)?