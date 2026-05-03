
test_that("getLInv - Case 1: Without metafounders", {
  
  ped <- pedigree(
    sire  = c(NA, NA, 1, 1, 4, 5),
    dam   = c(NA, NA, 2, NA, 3, 2),
    label = 1:6
  )
  
  LInvEsp <- getLInv(ped)
  LInvExp <- matrix(data = c(
    1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000,  1.0000,  0.0000,  0.0000,  0.0000, 0.0000,
    -0.7071, -0.7071,  1.4142,  0.0000,  0.0000, 0.0000,
    -0.5774,  0.0000,  0.0000,  1.1547,  0.0000, 0.0000,
    0.0000,  0.0000, -0.7071, -0.7071,  1.4142, 0.0000,
    0.0000, -0.7303,  0.0000,  0.0000, -0.7303, 1.4606),
    byrow = TRUE, nrow = 6)
  
  expect_equal(unname(as.matrix(LInvEsp)), unname(as.matrix(LInvExp)), tolerance = 1e-4)
  expect_equal(unname(as.matrix(crossprod(LInvEsp))), unname(as.matrix(getAInv(ped))), tolerance = 1e-6)
})

test_that("getLInv - Case 2: Single metafounder pedigree", {
  
  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )
  
  gamma_single <- matrix(0.02)
  
  LInvEst <- getLInv(ped_single, gamma = gamma_single)
  Ai <- matrix(c(
    52.020200, -1.010100, -1.010100,  0.000000,  0.000000,  0.000000,  0.000000,
    -1.010100,  2.020200,  1.010100, -1.010100, -1.010100,  0.000000,  0.000000,
    -1.010100,  1.010100,  2.525250, -0.505051, -1.010100,  0.000000, -1.010100,
    0.000000, -1.010100, -0.505051,  3.030300,  0.505051, -1.010100, -1.010100,
    0.000000, -1.010100, -1.010100,  0.505051,  2.525250, -1.010100,  0.000000,
    0.000000,  0.000000,  0.000000, -1.010100, -1.010100,  2.020200,  0.000000,
    0.000000,  0.000000, -1.010100, -1.010100,  0.000000,  0.000000,  2.020200
    ), nrow = 7, byrow = TRUE)
    
    dimnames(LInvEst) <- list(ped_single@label, ped_single@label)
    expect_equal(unname(as.matrix(crossprod(LInvEst))), unname(as.matrix(Ai)), tolerance = 1e-6)
})

test_that("getLInv - Case 3: Two metafounders pedigree", {
  ped_multi <- pedigree(
    sire  = c(0, 0, 1, 2, 3, 3, 5, 5),
    dam   = c(0, 0, 1, 2, 4, 4, 6, 4),
    label = 1:8
  )
  
  gamma_multi <- matrix(c(0.1, 0.05, 0.05, 0.2), nrow = 2, byrow = TRUE)
  
  Ai <- matrix(c(
    12.481200, -2.857140, -1.052630,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
    -2.857140,  6.825400,  0.000000, -1.111110,  0.000000,  0.000000,  0.000000,  0.000000,
    -1.052630,  0.000000,  2.133710,  1.081080, -1.081080, -1.081080,  0.000000,  0.000000,
    0.000000, -1.111110,  1.081080,  2.725530, -0.547748, -1.081080,  0.000000, -1.066670,
    0.000000,  0.000000, -1.081080, -0.547748,  3.208320,  0.512821, -1.025640, -1.066670,
    0.000000,  0.000000, -1.081080, -1.081080,  0.512821,  2.674980, -1.025640,  0.000000,
    0.000000,  0.000000,  0.000000,  0.000000, -1.025640, -1.025640,  2.051280,  0.000000,
    0.000000,  0.000000,  0.000000, -1.066670, -1.066670,  0.000000,  0.000000,  2.133330
  ), nrow = 8, byrow = TRUE)

  LInvEst <- getLInv(ped = ped_multi, gamma = gamma_multi)
  dimnames(LInvEst) <- list(ped_multi@label, ped_multi@label)
  expect_equal(unname(as.matrix(crossprod(LInvEst))), unname(as.matrix(Ai)), tolerance = 1e-5)
    
})

test_that("getLInv - Case 4: Two metafounders crosses between them", {
  
  ped_cross <- pedigree(
    sire  = c(0, 0, 1, 1, 3, 3, 5, 5),
    dam   = c(0, 0, 1, 2, 4, 4, 6, 4),
    label = 1:8
  )
  
  gamma_multi <- matrix(c(0.1, 0.05, 0.05, 0.2), nrow = 2, byrow = TRUE)
  
  LInvEst <- getLInv(ped = ped_cross, gamma = gamma_multi)
  Ai <- matrix(c(
    12.751500, -2.586870, -1.052630, -0.540541,  0.000000,  0.000000,  0.000000,  0.000000,
    -2.586870,  5.984560,  0.000000, -0.540541,  0.000000,  0.000000,  0.000000,  0.000000,
    -1.052630,  0.000000,  2.091590,  1.038960, -1.038960, -1.038960,  0.000000,  0.000000,
    -0.540541, -0.540541,  1.038960,  2.636170, -0.522832, -1.038960,  0.000000, -1.032260,
    0.000000,  0.000000, -1.038960, -0.522832,  3.113530,  0.519481, -1.038960, -1.032260,
    0.000000,  0.000000, -1.038960, -1.038960,  0.519481,  2.597400, -1.038960,  0.000000,
    0.000000,  0.000000,  0.000000,  0.000000, -1.038960, -1.038960,  2.077920,  0.000000,
    0.000000,  0.000000,  0.000000, -1.032260, -1.032260,  0.000000,  0.000000,  2.064520
  ), nrow = 8, byrow = TRUE)

  dimnames(LInvEst) <- list(ped_cross@label, ped_cross@label)
  expect_equal(unname(as.matrix(crossprod(LInvEst))), unname(as.matrix(Ai)), tolerance = 1e-5)
})

test_that("getLInv - Case 5: Two metafounders crosses between metafounder and individual", {
  
  ped_cross_ind <- pedigree(
    sire  = c(0, 0, 1, 1, 3, 3, 5, 5),
    dam   = c(0, 0, 1, 2, 4, 2, 6, 4),
    label = 1:8
  )
  
  gamma_multi <- matrix(c(0.1, 0.05, 0.05, 0.2), nrow =
     2, byrow = TRUE)
  
  LInvEst <- getLInv(ped = ped_cross_ind, gamma = gamma_multi)
  Ai <- matrix(c(
    12.751500, -2.586870, -1.052630, -0.540541,  0.000000,  0.000000,  0.000000,  0.000000,
    -2.586870,  6.348190,  0.363636, -0.540541,  0.000000, -0.727273,  0.000000,  0.000000,
    -1.052630,  0.363636,  1.935750,  0.519481, -1.038960, -0.727273,  0.000000,  0.000000,
    -0.540541, -0.540541,  0.519481,  2.116690, -0.522832,  0.000000,  0.000000, -1.032260,
    0.000000,  0.000000, -1.038960, -0.522832,  3.110180,  0.516129, -1.032260, -1.032260,
    0.000000, -0.727273, -0.727273,  0.000000,  0.516129,  1.970670, -1.032260,  0.000000,
    0.000000,  0.000000,  0.000000,  0.000000, -1.032260, -1.032260,  2.064520,  0.000000,
    0.000000,  0.000000,  0.000000, -1.032260, -1.032260,  0.000000,  0.000000,  2.064520
  ), nrow = 8, byrow = TRUE)

  dimnames(LInvEst) <- list(ped_cross_ind@label, ped_cross_ind@label)
  expect_equal(unname(as.matrix(crossprod(LInvEst))), unname(as.matrix(Ai)), tolerance = 1e-5)

})

# TODO: Test if getLInv with gamma = chol(AInv), this is not because something.. 