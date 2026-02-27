
# Testing pedigreeTools functions
# https://github.com/gregorgorjanc/pedigreeTools/blob/master/R/pedigree.R
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

#
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


test_that("Relationship factor from a pedigree; getL() and relfactor()",{
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

    expect_s4_class(L, "sparseMatrix")
})


test_that("Inverse relationship factor from a peidgree; getLInv() and relfactorInv()", {
    ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
                    dam =  c(NA, NA, 2, NA, 3, 2),
                    label = 1:6)
    LInv <- getLInv(ped)

    # Test for correctness
    LInvExp <- matrix(data = c(1.0000,  0.0000,  0.0000,  0.0000,  0.0000, 0.0000,
                               0.0000,  1.0000,  0.0000,  0.0000,  0.0000, 0.0000,
                               -0.7071, -0.7071,  1.4142,  0.0000,  0.0000, 0.0000,
                               -0.5774,  0.0000,  0.0000,  1.1547,  0.0000, 0.0000,
                               0.0000,  0.0000, -0.7071, -0.7071,  1.4142, 0.0000,
                               0.0000, -0.7303,  0.0000,  0.0000, -0.7303, 1.4606),
                      byrow = TRUE, nrow = 6)

    dimnames(LInvExp) <- list(as.character(1:6), as.character(1:6))

    expect_equal(as.matrix(LInv), LInvExp, tolerance = 1e-4)
    expect_s4_class(LInv, "sparseMatrix")

})


test_that("Inverse of the additive relationship matrix; getAinv()", {
    ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
                    dam =  c(NA, NA, 2, NA, 3, 2),
                    label = 1:6)
    AInv <- getAInv(ped)

    # Test for correctness
    AInvExp <- matrix(data = c( 1.833,  0.500, -1.000, -0.667,  0.000,  0.000,
                                0.500,  2.033, -1.000,  0.000,  0.533, -1.067,
                                -1.000, -1.000,  2.500,  0.500, -1.000,  0.000,
                                -0.667,  0.000,  0.500,  1.833, -1.000,  0.000,
                                0.000,  0.533, -1.000, -1.000,  2.533, -1.067,
                                0.000, -1.067,  0.000,  0.000, -1.067,  2.133),
                      byrow = TRUE, nrow = 6)

    dimnames(AInvExp) <- list(as.character(1:6), as.character(1:6))

    expect_equal(as.matrix(AInv), AInvExp, tolerance = 1e-3)
    expect_s4_class(AInv, "sparseMatrix")
    expect_true(isSymmetric(AInv, tol = 1e-12))
})


test_that("Additive relationship matrix; getA()", {
    ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
                    dam =  c(NA, NA, 2, NA, 3, 2),
                    label = 1:6)
    A <- getA(ped)
    # Test for correctness
    AExp <- matrix(data = c(1.0000, 0.0000, 0.5000, 0.5000, 0.5000, 0.2500,
                            0.0000, 1.0000, 0.5000, 0.0000, 0.2500, 0.6250,
                            0.5000, 0.5000, 1.0000, 0.2500, 0.6250, 0.5625,
                            0.5000, 0.0000, 0.2500, 1.0000, 0.6250, 0.3125,
                            0.5000, 0.2500, 0.6250, 0.6250, 1.1250, 0.6875,
                            0.2500, 0.6250, 0.5625, 0.3125, 0.6875, 1.1250),
                   byrow = TRUE, nrow = 6)

    dimnames(AExp) <- list(as.character(1:6), as.character(1:6))

    expect_equal(as.matrix(A), AExp, tolerance = 1e-4)
    expect_s4_class(A, "sparseMatrix")
    expect_true(isSymmetric(A, tol = 1e-12))
})


test_that("Subset of additive relationship matrix getASubset()",{
    ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
                    dam =  c(NA, NA, 2, NA, 3, 2),
                    label = 1:6)
    A <- getA(ped)
    ASubset  <- A[4:6, 4:6]
    ASubset2 <- getASubset(ped, labs = 4:6)
    ASubset3  <- A[6:4, 6:4]
    ASubset4 <- getASubset(ped, labs = 6:4)

    # Test for correctness
    expect_equal(as.matrix(ASubset),  as.matrix(ASubset2), tolerance = .Machine$double.eps)
    expect_equal(as.matrix(ASubset3), as.matrix(ASubset4),tolerance = .Machine$double.eps)

    expect_true(isSymmetric(ASubset2, tol = 1e-12))
    expect_true(isSymmetric(ASubset4, tol = 1e-12))
})


test_that("Counts number of generations of ancestors for one subject. Use recursion; getGenAncestors()", {
    ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
                    dam =  c(NA, NA, 2, NA, 3, 2),
                    label = 1:6)
    ped <- ped2DF(ped)
    ped$id <- row.names(ped)
    ped$generation <- NA

    tmp1 <- getGenAncestors(ped, id = 1)
    tmp2 <- getGenAncestors(ped, id = 4)
    tmp3 <- getGenAncestors(ped, id = 6)

    # Test for correctness
    # tmp1
    expect_equal(tmp1$generation[1], 0)
    expect_true(all(is.na(tmp1$generation[-1])))

    # tmp2
    expect_equal(tmp2$generation[c(1, 4)], c(0, 1))
    expect_true(all(is.na(tmp2$generation[-c(1, 4)])))

    # tmp3
    expect_equal(tmp3$generation, c(0, 0, 1, 1, 2, 3))
})





