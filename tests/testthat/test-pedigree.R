
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

test_that("Inverse gene flow from a pedigree", {

    ped <- pedigree(
        sire = c(NA, NA, 1, 1, 4, 5),
        dam  = c(NA, NA, 2, NA, 3, 2),
        label = 1:6
    )

    TInv <- getTInv(ped)

    TInvExp <- matrix(
        c(1.0,  0.0,  0.0,  0.0,  0.0,  0.0,
          0.0,  1.0,  0.0,  0.0,  0.0,  0.0,
          -0.5, -0.5,  1.0,  0.0,  0.0,  0.0,
          -0.5,  0.0,  0.0,  1.0,  0.0,  0.0,
          0.0,  0.0, -0.5, -0.5,  1.0,  0.0,
          0.0, -0.5,  0.0,  0.0, -0.5,  1.0),
        byrow = TRUE, nrow = 6
    )
    dimnames(TInvExp) <- list(as.character(1:6), as.character(1:6))

    expect_equal(as.matrix(TInv), TInvExp, tolerance = .Machine$double.eps)
    expect_s4_class(TInv, "sparseMatrix")
})


test_that("Gene flow from a pedigree", {
    ped <- pedigree(
        sire = c(NA, NA, 1, 1, 4, 5),
        dam  = c(NA, NA, 2, NA, 3, 2),
        label = 1:6
    )

    T <- getT(ped)

    # Test for correctness
    TExp <- matrix(
        c(1.00, 0.000, 0.00, 0.00, 0.0, 0,
          0.00, 1.000, 0.00, 0.00, 0.0, 0,
          0.50, 0.500, 1.00, 0.00, 0.0, 0,
          0.50, 0.000, 0.00, 1.00, 0.0, 0,
          0.50, 0.250, 0.50, 0.50, 1.0, 0,
          0.25, 0.625, 0.25, 0.25, 0.5, 1),
        byrow = TRUE, nrow = 6
    )
    dimnames(TExp) <- list(as.character(1:6), as.character(1:6))

    expect_equal(as.matrix(T), TExp, tolerance = .Machine$double.eps)
})


test_that("Relationship factor from a pedigree",{
    ped <- pedigree(
        sire = c(NA, NA, 1, 1, 4, 5),
        dam  = c(NA, NA, 2, NA, 3, 2),
        label = 1:6
    )

    L <- getL(ped)

    # Test for correctness (manual expected)
    LExp <- matrix(data = c(
        1.0000, 0.0000, 0.5000, 0.5000, 0.5000, 0.2500,
        0.0000, 1.0000, 0.5000, 0.0000, 0.2500, 0.6250,
        0.0000, 0.0000, 0.7071, 0.0000, 0.3536, 0.1768,
        0.0000, 0.0000, 0.0000, 0.8660, 0.4330, 0.2165,
        0.0000, 0.0000, 0.0000, 0.0000, 0.7071, 0.3536,
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.6847
    ), byrow = TRUE, nrow = 6)

    LExp <- as(Matrix::Matrix(LExp, sparse = TRUE), "dtCMatrix")
    dimnames(LExp) <- list(as.character(1:6), as.character(1:6))

    expect_equal(as.matrix(L), as.matrix(LExp), tolerance = 1e-4)

    # Alternative expected from A
    LExp2 <- chol(getA(ped))
    dimnames(LExp2) <- list(as.character(1:6), as.character(1:6))

    expect_equal(as.matrix(L), as.matrix(LExp2), tolerance = 1e-10)
})



