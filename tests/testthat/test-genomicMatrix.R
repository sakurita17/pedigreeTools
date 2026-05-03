genotypicMatrix12 <- matrix(c(
  1, 1, 2, 2,
  1, 2, 1, 1,
  1, 2, 1, 1,
  2, 1, 2, 2,
  0, 2, 2, 1,
  0, 2, 2, 1,
  1, 1, 2, 1,
  2, 1, 1, 2,
  1, 1, 1, 2,
  1, 1, 1, 2,
  1, 1, 2, 2,
  2, 0, 2, 2
), nrow = 12, byrow = TRUE,
dimnames = list(paste0("#", 1:12), paste0("Locus ", 1:4)))

genotypicMatrix12Ids <- rownames(genotypicMatrix12)

genotypicMatrix12G05Exp <- matrix(c(
  1.0000,  0.0000,  0.0000, 1.0000,  0.5000,  0.5000, 0.5000,  0.5000, 0.5000, 0.5000, 1.0000,  1.0000,
  0.0000,  0.5000,  0.5000, 0.0000,  0.5000,  0.5000, 0.0000,  0.0000, 0.0000, 0.0000, 0.0000, -0.5000,
  0.0000,  0.5000,  0.5000, 0.0000,  0.5000,  0.5000, 0.0000,  0.0000, 0.0000, 0.0000, 0.0000, -0.5000,
  1.0000,  0.0000,  0.0000, 1.5000,  0.0000,  0.0000, 0.5000,  1.0000, 0.5000, 0.5000, 1.0000,  1.5000,
  0.5000,  0.5000,  0.5000, 0.0000,  1.5000,  1.5000, 0.5000, -0.5000, 0.0000, 0.0000, 0.5000, -0.5000,
  0.5000,  0.5000,  0.5000, 0.0000,  1.5000,  1.5000, 0.5000, -0.5000, 0.0000, 0.0000, 0.5000, -0.5000,
  0.5000,  0.0000,  0.0000, 0.5000,  0.5000,  0.5000, 0.5000,  0.0000, 0.0000, 0.0000, 0.5000,  0.5000,
  0.5000,  0.0000,  0.0000, 1.0000, -0.5000, -0.5000, 0.0000,  1.0000, 0.5000, 0.5000, 0.5000,  1.0000,
  0.5000,  0.0000,  0.0000, 0.5000,  0.0000,  0.0000, 0.0000,  0.5000, 0.5000, 0.5000, 0.5000,  0.5000,
  0.5000,  0.0000,  0.0000, 0.5000,  0.0000,  0.0000, 0.0000,  0.5000, 0.5000, 0.5000, 0.5000,  0.5000,
  1.0000,  0.0000,  0.0000, 1.0000,  0.5000,  0.5000, 0.5000,  0.5000, 0.5000, 0.5000, 1.0000,  1.0000,
  1.0000, -0.5000, -0.5000, 1.5000, -0.5000, -0.5000, 0.5000,  1.0000, 0.5000, 0.5000, 1.0000,  2.0000
), nrow = 12, byrow = TRUE,
dimnames = list(genotypicMatrix12Ids, genotypicMatrix12Ids))

test_that("getG matches VR1 with 0.5 centering for the 12-individual example", {
  G <- getG(genotypicMatrix12, method = "VR1", freqMode = "0.5", epsilon = 1e-12)
  expect_equal(unname(G), unname(genotypicMatrix12G05Exp), tolerance = 1e-8)
  expect_identical(rownames(G), rownames(genotypicMatrix12))
  expect_identical(colnames(G), rownames(genotypicMatrix12))
})

test_that("getG returns the expected VR1 matrix with observed frequencies", {
  G_exp <- structure(c(
    0.2564, -0.4103, -0.4103, 0.2051, -0.1026, -0.1026,
    0.0000, -0.0513, 0.0000, 0.0000, 0.2564, 0.3590,
   -0.4103, 0.7692, 0.7692, -0.4615, 0.4615, 0.4615,
   -0.0513, -0.1026, -0.0513, -0.0513, -0.4103, -0.9231,
   -0.4103, 0.7692, 0.7692, -0.4615, 0.4615, 0.4615,
   -0.0513, -0.1026, -0.0513, -0.0513, -0.4103, -0.9231,
    0.2051, -0.4615, -0.4615, 0.7692, -0.7692, -0.7692,
   -0.0513, 0.5128, -0.0513, -0.0513, 0.2051, 0.9231,
   -0.1026, 0.4615, 0.4615, -0.7692, 1.3846, 1.3846,
    0.2564, -1.0256, -0.3590, -0.3590, -0.1026, -1.2308,
   -0.1026, 0.4615, 0.4615, -0.7692, 1.3846, 1.3846,
    0.2564, -1.0256, -0.3590, -0.3590, -0.1026, -1.2308,
    0.0000, -0.0513, -0.0513, -0.0513, 0.2564, 0.2564,
    0.3590, -0.3077, -0.2564, -0.2564, 0.0000, 0.1026,
   -0.0513, -0.1026, -0.1026, 0.5128, -1.0256, -1.0256,
   -0.3077, 0.8718, 0.3077, 0.3077, -0.0513, 0.6667,
    0.0000, -0.0513, -0.0513, -0.0513, -0.3590, -0.3590,
   -0.2564, 0.3077, 0.3590, 0.3590, 0.0000, 0.1026,
    0.0000, -0.0513, -0.0513, -0.0513, -0.3590, -0.3590,
   -0.2564, 0.3077, 0.3590, 0.3590, 0.0000, 0.1026,
    0.2564, -0.4103, -0.4103, 0.2051, -0.1026, -0.1026,
    0.0000, -0.0513, 0.0000, 0.0000, 0.2564, 0.3590,
    0.3590, -0.9231, -0.9231, 0.9231, -1.2308, -1.2308,
    0.1026, 0.6667, 0.1026, 0.1026, 0.3590, 1.6923
  ), dim = c(12L, 12L),
  dimnames = list(genotypicMatrix12Ids, genotypicMatrix12Ids))

  G <- getG(genotypicMatrix12, method = "VR1", freqMode = "observed", epsilon = 1e-12)

  expect_equal(unname(G), unname(G_exp), tolerance = 1e-4)
  expect_identical(rownames(G), genotypicMatrix12Ids)
  expect_identical(colnames(G), genotypicMatrix12Ids)
})

test_that("getG returns the expected VR2 matrix with observed frequencies", {
  G_exp <- structure(c(
    0.3000, -0.4649, -0.4649, 0.2580, -0.1072, -0.1072,
   -0.0158, -0.0578, -0.0158, -0.0158, 0.3000, 0.3914,
   -0.4649, 0.8193, 0.8193, -0.5069, 0.4191, 0.4191,
   -0.0228, -0.0648, -0.0228, -0.0228, -0.4649, -0.9069,
   -0.4649, 0.8193, 0.8193, -0.5069, 0.4191, 0.4191,
   -0.0228, -0.0648, -0.0228, -0.0228, -0.4649, -0.9069,
    0.2580, -0.5069, -0.5069, 0.7196, -0.6526, -0.6526,
   -0.0578, 0.4038, -0.0578, -0.0578, 0.2580, 0.8529,
   -0.1072, 0.4191, 0.4191, -0.6526, 1.2804, 1.2804,
    0.3349, -0.9684, -0.4230, -0.4230, -0.1072, -1.0526,
   -0.1072, 0.4191, 0.4191, -0.6526, 1.2804, 1.2804,
    0.3349, -0.9684, -0.4230, -0.4230, -0.1072, -1.0526,
   -0.0158, -0.0228, -0.0228, -0.0578, 0.3349, 0.3349,
    0.4263, -0.3735, -0.3316, -0.3316, -0.0158, 0.0756,
   -0.0578, -0.0648, -0.0648, 0.4038, -0.9684, -0.9684,
   -0.3735, 0.8459, 0.3843, 0.3843, -0.0578, 0.5371,
   -0.0158, -0.0228, -0.0228, -0.0578, -0.4230, -0.4230,
   -0.3316, 0.3843, 0.4263, 0.4263, -0.0158, 0.0756,
   -0.0158, -0.0228, -0.0228, -0.0578, -0.4230, -0.4230,
   -0.3316, 0.3843, 0.4263, 0.4263, -0.0158, 0.0756,
    0.3000, -0.4649, -0.4649, 0.2580, -0.1072, -0.1072,
   -0.0158, -0.0578, -0.0158, -0.0158, 0.3000, 0.3914,
    0.3914, -0.9069, -0.9069, 0.8529, -1.0526, -1.0526,
    0.0756, 0.5371, 0.0756, 0.0756, 0.3914, 1.5196
  ), dim = c(12L, 12L),
  dimnames = list(genotypicMatrix12Ids, genotypicMatrix12Ids))

  G <- getG(genotypicMatrix12, method = "VR2", freqMode = "observed", epsilon = 1e-12)

  expect_equal(unname(G), unname(G_exp), tolerance = 1e-4)
  expect_identical(rownames(G), genotypicMatrix12Ids)
  expect_identical(colnames(G), genotypicMatrix12Ids)
})

test_that("getG matches VR2 with 0.5 centering for the 12-individual example", {
  G <- getG(genotypicMatrix12, method = "VR2", freqMode = "0.5", epsilon = 1e-12)

  expect_equal(unname(G), unname(genotypicMatrix12G05Exp), tolerance = 1e-8)
  expect_identical(rownames(G), rownames(genotypicMatrix12))
  expect_identical(colnames(G), rownames(genotypicMatrix12))
})

test_that("getG supports big.matrix-like genotype inputs", {
  big_like <- genotypicMatrix12
  class(big_like) <- c("big.matrix", class(big_like))

  G_matrix <- getG(genotypicMatrix12, method = "VR1", freqMode = "observed", epsilon = 1e-12)
  G_big <- getG(big_like, method = "VR1", freqMode = "observed", epsilon = 1e-12)

  expect_equal(unname(G_big), unname(G_matrix), tolerance = 1e-8)
  expect_identical(rownames(G_big), rownames(genotypicMatrix12))
  expect_identical(colnames(G_big), rownames(genotypicMatrix12))
})
