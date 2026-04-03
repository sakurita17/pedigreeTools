test_that("Relationship factor from a pedigree - Case 1: Without metafounders ", {

    ped <- pedigree(
        sire  = c(NA, NA, 1, 1, 4, 5),
        dam   = c(NA, NA, 2, NA, 3, 2),
        label = 1:6
    )

    LEst <- getL(ped)
    AChol <- chol(getA(ped))
    LExp <- matrix(data = c(1.0000, 0.0000, 0.5000, 0.5000, 0.5000, 0.2500,
                            0.0000, 1.0000, 0.5000, 0.0000, 0.2500, 0.6250,
                            0.0000, 0.0000, 0.7071, 0.0000, 0.3536, 0.1768,
                            0.0000, 0.0000, 0.0000, 0.8660, 0.4330, 0.2165,
                            0.0000, 0.0000, 0.0000, 0.0000, 0.7071, 0.3536,
                            0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.6847),
                   byrow = TRUE, nrow = 6)
    LExp <- as(Matrix::Matrix(LExp, sparse = TRUE), "dtCMatrix")
    expect_equal(unname(as.matrix(LEst)), unname(as.matrix(LExp)), tolerance = 1e-4)
    expect_equal(unname(as.matrix(LEst)), unname(as.matrix(AChol)), tolerance = 1e-4)
})


