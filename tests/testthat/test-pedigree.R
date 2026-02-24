test_that("ped2DF(), pedigree to data frame", {
  ped <- pedigree(
      sire = c(NA, NA, 1, 1, 4, 5),
      dam = c(NA, NA, 2, NA, 3, 2),
      label = 1:6
  )

  df <- ped2DF(ped)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 6)
  expect_true(all(c("sire", "dam") %in% names(df)))

  expect_error(ped2DF())
})
