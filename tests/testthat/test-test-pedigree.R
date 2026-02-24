test_that("test-pedigree", {
  ped <- pedigree(
      sire = c(NA, NA, 1, 1, 4, 5),
      dam = c(NA, NA, 2, NA, 3, 2),
      label = 1:6
  )

  F <- inbreeding(ped)


})
