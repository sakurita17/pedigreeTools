test_that("test inbreeding()", {
  ped <- pedigree(
      sire = c(NA, NA, 1, 1, 4, 5),
      dam = c(NA, NA, 2, NA, 3, 2),
      label = 1:6
  )

  FEst <- inbreeding(ped)
  FExp <- c(0.000, 0.000, 0.000, 0.000, 0.125, 0.125)

  expect_equal(FEst, FExp, tolerance = .Machine$double.eps)
})
