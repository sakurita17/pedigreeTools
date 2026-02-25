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

test_that("test Mendelian sampling variance: getD() & getDinv()", {
    ped <- pedigree(
        sire = c(NA, NA, 1, 1, 4, 5),
        dam = c(NA, NA, 2, NA, 3, 2),
        label = 1:6
    )

    D <- getD(ped)
    DInv <- getDInv(ped)

    # Test for correctness
    DExp <- c(1.00, 1.00, 0.50, 0.75, 0.50, 0.46875)
    names(DExp) <- as.character(1:6)
    expect_equal(D, DExp, tolerance = .Machine$double.eps)

    DInvExp <- 1 / DExp
    names(DInvExp) <- as.character(1:6)
    expect_equal(DInv, DInvExp, tolerance = .Machine$double.eps)
})



