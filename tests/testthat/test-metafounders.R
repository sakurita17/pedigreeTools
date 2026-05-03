metMatch <- c(
  "#1" = "M1", "#2" = "M1", "#3" = "M2", "#4" = "M2",
  "#5" = "M1", "#6" = "M2", "#7" = "M2", "#8" = "M2",
  "#9" = "M1", "#10" = "M1", "#11" = "M1", "#12" = "M1"
)


test_that("masuda_implicit is linear in q", {
  masuda_implicit <- pedigreeTools:::masuda_implicit

  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  gamma_single <- matrix(0.02)
  idx.geno <- as.character(4:7)
  q1 <- c(1, -1, 0.5, 2)
  q2 <- c(-0.25, 0.75, 1, -1.5)

  Li <- getLInv(ped_single, gamma = gamma_single)
  m1 <- masuda_implicit(q1, Li, idx.geno)
  m2 <- masuda_implicit(q2, Li, idx.geno)
  m12 <- masuda_implicit(q1 + q2, Li, idx.geno)

  expect_equal(
    unname(as.matrix(m12)),
    unname(as.matrix(m1 + m2)),
    tolerance = 1e-8
  )
})

test_that("masuda_implicit accepts marker matrices in block form", {
  masuda_implicit <- pedigreeTools:::masuda_implicit

  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  gamma_single <- matrix(0.02)
  idx.geno <- as.character(4:7)
  Z <- matrix(c(
    -1,  0,  1,
     0,  0, -1,
     1, -1,  0,
     0,  1,  0
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  Li <- getLInv(ped_single, gamma = gamma_single)
  block <- masuda_implicit(Z, Li, idx.geno)
  by_col <- lapply(seq_len(ncol(Z)), function(j) {
    masuda_implicit(Z[, j], Li, idx.geno)
  })
  expected <- do.call(cbind, by_col)

  expect_equal(
    unname(as.matrix(block)),
    unname(as.matrix(expected)),
    tolerance = 1e-8
  )
})

test_that("logdet_block_A22_from_Li matches the block determinant identity", {
  logdet_block_A22_from_Li <- pedigreeTools:::logdet_block_A22_from_Li

  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  gamma_single <- matrix(0.02)
  idx.geno <- as.character(4:7)
  Li <- getLInv(ped_single, gamma = gamma_single)
  Ainv <- Matrix::crossprod(Li)
  idx.non.geno <- setdiff(rownames(Ainv), idx.geno)

  det_Ainv <- determinant(Ainv, logarithm = TRUE)
  det_A11 <- determinant(Ainv[idx.non.geno, idx.non.geno, drop = FALSE], logarithm = TRUE)
  expected <- -as.numeric(det_Ainv$modulus) + as.numeric(det_A11$modulus)

  expect_equal(
    logdet_block_A22_from_Li(Li, idx.geno),
    expected,
    tolerance = 1e-8
  )
})

test_that("logLikeGamma agrees by marker and by individual", {
  logLikeGamma <- pedigreeTools:::logLikeGamma

  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  gamma_single <- matrix(0.02)
  idx.geno <- as.character(4:7)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  ll_marker <- logLikeGamma(ped_single, gamma_single, idx.geno, genotypicMatrix, method = "byMarker")
  ll_individual <- logLikeGamma(ped_single, gamma_single, idx.geno, genotypicMatrix, method = "byIndividual")

  expect_equal(ll_marker, ll_individual, tolerance = 1e-8)
})

test_that("single gamma components reconstruct the Legarra likelihood", {
  components_fn <- function(ped, gamma_value, genotypedIds, Z, method) {
    pedigreeTools:::legarra_gamma_components(
      ped = ped,
      gamma = matrix(gamma_value, nrow = 1L),
      genotypedIds = genotypedIds,
      Z = Z,
      method = method
    )
  }
  prepare_legarra_data <- pedigreeTools:::prepare_legarra_data
  logLikeGamma <- pedigreeTools:::logLikeGamma
  masuda_implicit <- pedigreeTools:::masuda_implicit

  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  idx.geno <- as.character(4:7)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  prepared <- prepare_legarra_data(genotypicMatrix, ped_single, genotypedIds = idx.geno)
  by_marker <- components_fn(
    ped = ped_single,
    gamma_value = 0.02,
    genotypedIds = idx.geno,
    Z = prepared$Z,
    method = "byMarker"
  )
  by_individual <- components_fn(
    ped = ped_single,
    gamma_value = 0.02,
    genotypedIds = idx.geno,
    Z = prepared$Z,
    method = "byIndividual"
  )

  expect_equal(
    by_marker$logLik,
    logLikeGamma(ped_single, matrix(0.02), idx.geno, genotypicMatrix, method = "byMarker"),
    tolerance = 1e-8
  )
  expect_equal(
    by_individual$logLik,
    logLikeGamma(ped_single, matrix(0.02), idx.geno, genotypicMatrix, method = "byIndividual"),
    tolerance = 1e-8
  )
  expect_equal(
    by_marker$logLik,
    by_marker$determinant_term + by_marker$quadratic_term,
    tolerance = 1e-8
  )
  expect_equal(
    by_individual$logLik,
    by_individual$determinant_term + by_individual$quadratic_term,
    tolerance = 1e-8
  )

  Li <- getLInv(ped_single, gamma = matrix(0.02))
  A22inv_Z <- masuda_implicit(prepared$Z, Li, idx.geno)
  quad_by_marker <- colSums(prepared$Z * A22inv_Z)

  expect_equal(
    by_individual$trace_term,
    (2 / ncol(prepared$Z)) * sum(quad_by_marker),
    tolerance = 1e-8
  )
})

test_that("logLikeGamma rejects non-positive-definite gamma matrices", {
  logLikeGamma <- pedigreeTools:::logLikeGamma

  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  idx.geno <- as.character(4:7)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  expect_error(
    logLikeGamma(
      ped = ped_single,
      gamma = matrix(c(1, 2, 2, 1), nrow = 2),
      genotypedIds = idx.geno,
      genotypicMatrix = genotypicMatrix
    ),
    "gamma must be a symmetric positive-definite matrix"
  )
})

test_that("getGamma auto-detects centered -1, 0, 1 coding for Legarra", {
  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  idx.geno <- as.character(4:7)
  centered_matrix <- matrix(c(
    -1,  0,  1,
     0,  0, -1,
     1, -1,  0,
     0,  1,  0
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))
  uncentered_matrix <- centered_matrix + 1
  fit_centered <- getGamma(
    ped = ped_single,
    genotypicMatrix = centered_matrix,
    groups = 1,
    method = "legarra"
  )
  fit_uncentered <- getGamma(
    ped = ped_single,
    genotypicMatrix = uncentered_matrix,
    groups = 1,
    method = "legarra"
  )

  expect_equal(fit_centered$gamma, fit_uncentered$gamma, tolerance = 1e-8)
  expect_equal(fit_centered$logLik, fit_uncentered$logLik, tolerance = 1e-8)
})

test_that("getGamma auto-detects centered -1, 0, 1 coding for Garcia-Baccino", {
  genotypicMatrix <- matrix(c(
    0, 2,
    2, 2,
    0, 0,
    2, 0
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("i", 1:4), paste0("m", 1:2)))
  centered_matrix <- genotypicMatrix - 1
  metMatch <- c(i1 = "MF1", i2 = "MF1", i3 = "MF2", i4 = "MF2")

  fit_uncentered <- getGamma(genotypicMatrix = genotypicMatrix, metMatch = metMatch, method = "garcia-baccino")
  fit_centered <- getGamma(genotypicMatrix = centered_matrix, metMatch = metMatch, method = "garcia-baccino")

  expect_equal(fit_centered$gamma, fit_uncentered$gamma, tolerance = 1e-8)
  expect_equal(fit_centered$pCentered, fit_uncentered$pCentered, tolerance = 1e-8)
})

test_that("Garcia-Baccino supports big.matrix-like genotype inputs", {
  genotypicMatrix <- matrix(c(
    0, 2,
    2, 2,
    0, 0,
    2, 0
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("i", 1:4), paste0("m", 1:2)))
  big_like <- genotypicMatrix
  class(big_like) <- c("big.matrix", class(big_like))
  metMatch <- c(i1 = "MF1", i2 = "MF1", i3 = "MF2", i4 = "MF2")

  fit_matrix <- getGammaGarciaBaccino(genotypicMatrix = genotypicMatrix, metMatch = metMatch)
  fit_big <- getGammaGarciaBaccino(genotypicMatrix = big_like, metMatch = metMatch)

  expect_equal(unname(fit_big$gamma), unname(fit_matrix$gamma), tolerance = 1e-8)
  expect_equal(unname(fit_big$pCentered), unname(fit_matrix$pCentered), tolerance = 1e-8)
})

test_that("pseudo-EM update returns a symmetric positive-definite gamma", {
  update_fn <- pedigreeTools:::pseudo_em_gamma_update
  prepare_legarra_data <- pedigreeTools:::prepare_legarra_data

  ped_multi <- pedigree(
    sire  = c(0, 0, 1, 2, 3, 3, 5, 5),
    dam   = c(0, 0, 1, 2, 4, 4, 6, 4),
    label = 1:8
  )

  idx.geno <- as.character(5:8)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  prepared <- prepare_legarra_data(genotypicMatrix, ped_multi, genotypedIds = idx.geno)
  start_gamma <- diag(0.1, 2)
  dimnames(start_gamma) <- list(as.character(1:2), as.character(1:2))
  update <- update_fn(
    ped = ped_multi,
    gamma = start_gamma,
    genotypedIds = idx.geno,
    Z = prepared$Z
  )

  expect_true(isSymmetric(update$gamma_next))
  expect_true(isTRUE(all(eigen(update$gamma_next, symmetric = TRUE, only.values = TRUE)$values > 0)))
})

test_that("Legarra supports big.matrix-like genotype inputs", {
  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  idx.geno <- as.character(4:7)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))
  big_like <- genotypicMatrix
  class(big_like) <- c("big.matrix", class(big_like))

  fit_matrix <- getGammaLegarra(
    ped = ped_single,
    genotypicMatrix = genotypicMatrix,
    groups = 1,
    genotypedIds = idx.geno,
    idSource = "genotyped"
  )
  fit_big <- getGammaLegarra(
    ped = ped_single,
    genotypicMatrix = big_like,
    groups = 1,
    genotypedIds = idx.geno,
    idSource = "genotyped"
  )

  expect_equal(fit_big$gamma, fit_matrix$gamma, tolerance = 1e-8)
  expect_equal(fit_big$logLik, fit_matrix$logLik, tolerance = 1e-8)
})

test_that("getGamma returns a finite estimate for one metafounder", {
  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  idx.geno <- as.character(4:7)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  fit_marker <- getGamma(
    ped = ped_single,
    genotypicMatrix = genotypicMatrix,
    groups = 1,
    method = "legarra",
    likelihood = "byMarker"
  )

  fit_individual <- getGamma(
    ped = ped_single,
    genotypicMatrix = genotypicMatrix,
    groups = 1,
    method = "legarra",
    likelihood = "byIndividual"
  )

  expect_true(is.finite(fit_marker$gamma[1, 1]))
  expect_true(is.finite(fit_marker$logLik))
  expect_true(fit_marker$gamma[1, 1] > 0 && fit_marker$gamma[1, 1] <= 1)
  expect_equal(fit_marker$gamma, fit_individual$gamma, tolerance = 1e-6)
  expect_equal(fit_marker$genotyped, idx.geno)
  expect_equal(fit_marker$id_source, "rownames")
})

test_that("getGammaLegarra matches the compatibility wrapper", {
  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  idx.geno <- as.character(4:7)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  fit_specific <- getGammaLegarra(
    ped = ped_single,
    genotypicMatrix = genotypicMatrix,
    groups = 1,
    likelihood = "byMarker"
  )
  fit_wrapper <- getGamma(
    ped = ped_single,
    genotypicMatrix = genotypicMatrix,
    groups = 1,
    method = "legarra",
    likelihood = "byMarker"
  )

  expect_equal(fit_specific$gamma, fit_wrapper$gamma, tolerance = 1e-8)
  expect_equal(fit_specific$logLik, fit_wrapper$logLik, tolerance = 1e-8)
  expect_equal(fit_specific$optimizer, fit_wrapper$optimizer)
})

test_that("getGamma keeps the original closed-form solver for one metafounder", {
  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  idx.geno <- as.character(4:7)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  fit_default <- getGamma(
    ped = ped_single,
    genotypicMatrix = genotypicMatrix,
    groups = 1,
    method = "legarra"
  )
  fit_repeat <- getGamma(
    ped = ped_single,
    genotypicMatrix = genotypicMatrix,
    groups = 1,
    method = "legarra"
  )

  expect_equal(fit_default$gamma, fit_repeat$gamma, tolerance = 1e-8)
  expect_equal(fit_default$logLik, fit_repeat$logLik, tolerance = 1e-8)
  expect_equal(fit_default$optimizer, "closed-form")
  expect_equal(fit_repeat$optimizer, "closed-form")
})

test_that("getGamma matches the closed-form single-metafounder solver", {
  single_solver <- pedigreeTools:::legarra_single_exact_ml
  single_selector <- pedigreeTools:::select_legarra_single_gamma
  prepare_legarra_data <- pedigreeTools:::prepare_legarra_data
  logLikeGamma <- pedigreeTools:::logLikeGamma

  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  idx.geno <- as.character(4:7)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  fit <- getGamma(
    ped = ped_single,
    genotypicMatrix = genotypicMatrix,
    groups = 1,
    method = "legarra",
    likelihood = "byMarker"
  )

  prepared <- prepare_legarra_data(genotypicMatrix, ped_single, genotypedIds = idx.geno)
  roots <- single_solver(
    ped = ped_single,
    genotypedIds = idx.geno,
    Z = prepared$Z
  )
  reference <- single_selector(
    ped = ped_single,
    genotypedIds = idx.geno,
    Z = prepared$Z
  )

  expect_s3_class(roots, "legarra_single_gamma_candidates")
  expect_equal(fit$gamma[1, 1], reference$gamma, tolerance = 1e-12)
  expect_equal(fit$logLik, reference$logLik, tolerance = 1e-10)
  expect_equal(names(roots), c("real_roots", "admissible_roots", "candidates"))
  expect_false(any(c("gamma", "logLik") %in% names(roots)))
  expect_true(all(roots$admissible_roots >= 1e-6 & roots$admissible_roots <= 1))
  expect_equal(
    names(roots$candidates),
    c("gamma", "logLik", "is_root", "is_boundary")
  )
  expect_true(isTRUE(all(diff(roots$candidates$logLik) <= 0)))
  expect_equal(roots$candidates$gamma[1], reference$gamma, tolerance = 1e-12)
  expect_equal(roots$candidates$logLik[1], reference$logLik, tolerance = 1e-10)
  expect_match(capture_output(print(roots)), "Recommended gamma:")
  expect_equal(
    fit$logLik,
    logLikeGamma(
      ped = ped_single,
      gamma = fit$gamma,
      genotypedIds = idx.geno,
      genotypicMatrix = genotypicMatrix,
      method = "byMarker"
    ),
    tolerance = 1e-8
  )
})

test_that("getGamma can infer genotyped IDs from rownames(genotypicMatrix)", {
  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  idx.geno <- as.character(4:7)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  fit <- getGamma(
    ped = ped_single,
    genotypicMatrix = genotypicMatrix,
    groups = 1,
    method = "legarra"
  )

  expect_equal(fit$genotyped, idx.geno)
  expect_equal(fit$id_source, "rownames")
})

test_that("getGamma can use genotypedIds by row order when requested", {
  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  idx.geno <- as.character(4:7)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("x", 1:4), paste0("m", 1:3)))

  fit <- getGamma(
    ped = ped_single,
    genotypicMatrix = genotypicMatrix,
    groups = 1,
    genotypedIds = idx.geno,
    idSource = "genotyped",
    method = "legarra"
  )

  expect_equal(fit$genotyped, idx.geno)
  expect_equal(fit$id_source, "genotyped")
  expect_true(is.finite(fit$logLik))
})

test_that("getGamma errors if Legarra IDs cannot be resolved", {
  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("x", 1:4), paste0("m", 1:3)))

  expect_error(
    getGamma(ped = ped_single, genotypicMatrix = genotypicMatrix, method = "legarra"),
    "provide rownames\\(genotypicMatrix\\) matching ped@label or supply genotypedIds"
  )
})

test_that("getGamma requires groups for Legarra", {
  ped_single <- pedigree(
    sire  = c(0, 1, 1, 2, 2, 4, 4),
    dam   = c(0, 1, 1, 3, 3, 5, 3),
    label = 1:7
  )

  idx.geno <- as.character(4:7)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  expect_error(
    getGamma(ped = ped_single, genotypicMatrix = genotypicMatrix, method = "legarra"),
    "groups must be a single positive integer for method = 'legarra'"
  )
})

test_that("getGamma validates pseudo-EM control entries", {
  ped_multi <- pedigree(
    sire  = c(0, 0, 1, 2, 3, 3, 5, 5),
    dam   = c(0, 0, 1, 2, 4, 4, 6, 4),
    label = 1:8
  )

  idx.geno <- as.character(5:8)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  expect_error(
    getGamma(
      ped = ped_multi,
      genotypicMatrix = genotypicMatrix,
      groups = 2,
      method = "legarra",
      control = list(100)
    ),
    "control must be a named list"
  )

  expect_error(
    getGamma(
      ped = ped_multi,
      genotypicMatrix = genotypicMatrix,
      groups = 2,
      method = "legarra",
      control = list(foo = 100)
    ),
    "Unknown control entries"
  )
})

test_that("getGamma supports pseudo-EM optimization with two metafounders", {
  logLikeGamma <- pedigreeTools:::logLikeGamma

  ped_multi <- pedigree(
    sire  = c(0, 0, 1, 2, 3, 3, 5, 5),
    dam   = c(0, 0, 1, 2, 4, 4, 6, 4),
    label = 1:8
  )

  idx.geno <- as.character(5:8)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))
  fit <- getGamma(
    ped = ped_multi,
    genotypicMatrix = genotypicMatrix,
    groups = 2,
    method = "legarra",
    likelihood = "byIndividual",
    control = list(
      maxit = 200,
      tol = 1e-10,
      keep_history = TRUE,
      compute_logLik = TRUE
    )
  )

  expected_ll <- logLikeGamma(
    ped = ped_multi,
    gamma = fit$gamma,
    genotypedIds = idx.geno,
    genotypicMatrix = genotypicMatrix,
    method = "byIndividual"
  )

  expect_equal(dim(fit$gamma), c(2, 2))
  expect_true(all(is.finite(fit$gamma)))
  expect_true(isTRUE(all(eigen(fit$gamma, symmetric = TRUE, only.values = TRUE)$values > 0)))
  expect_true(is.finite(fit$logLik))
  expect_equal(fit$logLik, expected_ll, tolerance = 1e-8)
  expect_equal(fit$optimizer, "pseudo-em")
  expect_true(nrow(fit$history) >= 1L)
  expect_true("logLik" %in% names(fit$history))
})

test_that("getGamma uses pseudo-EM by default for two metafounders", {
  ped_multi <- pedigree(
    sire  = c(0, 0, 1, 2, 3, 3, 5, 5),
    dam   = c(0, 0, 1, 2, 4, 4, 6, 4),
    label = 1:8
  )

  idx.geno <- as.character(5:8)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  fit <- getGamma(
    ped = ped_multi,
    genotypicMatrix = genotypicMatrix,
    groups = 2,
    method = "legarra",
    likelihood = "byIndividual",
    control = list(maxit = 50, tol = 1e-8)
  )

  expect_equal(dim(fit$gamma), c(2, 2))
  expect_true(all(is.finite(fit$gamma)))
  expect_true(isTRUE(all(eigen(fit$gamma, symmetric = TRUE, only.values = TRUE)$values > 0)))
  expect_equal(fit$optimizer, "pseudo-em")
  expect_null(fit$logLik)
})

test_that("getGamma uses the internal default multi-metafounder start", {
  ped_multi <- pedigree(
    sire  = c(0, 0, 1, 2, 3, 3, 5, 5),
    dam   = c(0, 0, 1, 2, 4, 4, 6, 4),
    label = 1:8
  )

  idx.geno <- as.character(5:8)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  fit_default <- getGamma(
    ped = ped_multi,
    genotypicMatrix = genotypicMatrix,
    groups = 2,
    method = "legarra",
    likelihood = "byIndividual",
    control = list(maxit = 50, tol = 1e-8)
  )

  fit_repeat <- getGamma(
    ped = ped_multi,
    genotypicMatrix = genotypicMatrix,
    groups = 2,
    method = "legarra",
    likelihood = "byIndividual",
    control = list(maxit = 50, tol = 1e-8)
  )

  expect_equal(fit_default$gamma, fit_repeat$gamma, tolerance = 1e-8)
  expect_equal(fit_default$criterion, fit_repeat$criterion, tolerance = 1e-12)
  expect_equal(fit_default$convergence, fit_repeat$convergence)
})

test_that("getGammaLegarra supports the multi-metafounder pseudo-EM API directly", {
  ped_multi <- pedigree(
    sire  = c(0, 0, 1, 2, 3, 3, 5, 5),
    dam   = c(0, 0, 1, 2, 4, 4, 6, 4),
    label = 1:8
  )

  idx.geno <- as.character(5:8)
  genotypicMatrix <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE,
  dimnames = list(idx.geno, paste0("m", 1:3)))

  fit_specific <- getGammaLegarra(
    ped = ped_multi,
    genotypicMatrix = genotypicMatrix,
    groups = 2,
    likelihood = "byIndividual",
    control = list(maxit = 50, tol = 1e-8)
  )
  fit_wrapper <- getGamma(
    ped = ped_multi,
    genotypicMatrix = genotypicMatrix,
    groups = 2,
    method = "legarra",
    likelihood = "byIndividual",
    control = list(maxit = 50, tol = 1e-8)
  )

  expect_equal(fit_specific$gamma, fit_wrapper$gamma, tolerance = 1e-8)
  expect_equal(fit_specific$criterion, fit_wrapper$criterion, tolerance = 1e-12)
  expect_equal(fit_specific$convergence, fit_wrapper$convergence)
})

test_that("getGamma supports Garcia-Baccino from grouped genotypes", {
  genotypicMatrix <- matrix(c(
    0, 2,
    2, 2,
    0, 0,
    2, 0
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("i", 1:4), paste0("m", 1:2)))
  metMatch <- c(i1 = "MF1", i2 = "MF1", i3 = "MF2", i4 = "MF2")

  fit <- getGamma(genotypicMatrix = genotypicMatrix, metMatch = metMatch, method = "garcia-baccino")

  p_exp <- rbind(
    MF1 = c(0.5, 1.0),
    MF2 = c(0.5, 0.0)
  )
  p_centered_exp <- p_exp - 0.5
  gamma_exp <- (8 / 2) * tcrossprod(p_centered_exp)

  expect_equal(unname(fit$pCentered), unname(p_centered_exp), tolerance = 1e-8)
  expect_equal(unname(fit$alleleFrequencies), unname(p_exp), tolerance = 1e-8)
  expect_equal(unname(fit$gamma), unname(gamma_exp), tolerance = 1e-8)
  expect_equal(unname(fit$metafounders), c("MF1", "MF2"))
  expect_equal(fit$estimator, "naive")
})

test_that("getGammaGarciaBaccino matches the compatibility wrapper", {
  genotypicMatrix <- matrix(c(
    0, 2,
    2, 2,
    0, 0,
    2, 0
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("i", 1:4), paste0("m", 1:2)))
  metMatch <- c(i1 = "MF1", i2 = "MF1", i3 = "MF2", i4 = "MF2")

  fit_specific <- getGammaGarciaBaccino(
    genotypicMatrix = genotypicMatrix,
    metMatch = metMatch
  )
  fit_wrapper <- getGamma(
    genotypicMatrix = genotypicMatrix,
    metMatch = metMatch,
    method = "garcia-baccino"
  )

  expect_equal(fit_specific$gamma, fit_wrapper$gamma, tolerance = 1e-8)
  expect_equal(fit_specific$pCentered, fit_wrapper$pCentered, tolerance = 1e-8)
  expect_equal(fit_specific$alleleFrequencies, fit_wrapper$alleleFrequencies, tolerance = 1e-8)
})

test_that("getGamma requires genotypicMatrix for Garcia-Baccino", {
  expect_error(
    getGamma(method = "garcia-baccino"),
    "genotypicMatrix must be provided for method = 'garcia-baccino'"
  )
})

test_that("Garcia-Baccino validates metMatch IDs against rownames(genotypicMatrix)", {
  genotypicMatrix <- matrix(c(
    0, 2,
    2, 2,
    0, 0,
    2, 0
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("i", 1:4), paste0("m", 1:2)))

  expect_error(
    getGammaGarciaBaccino(
      genotypicMatrix = genotypicMatrix,
      metMatch = c(i1 = "MF1", i2 = "MF1", i3 = "MF2", i5 = "MF2")
    ),
    "All IDs in metMatch must be present in rownames\\(genotypicMatrix\\)"
  )

  expect_error(
    getGammaGarciaBaccino(
      genotypicMatrix = genotypicMatrix,
      metMatch = c(i1 = "MF1", i2 = "MF1", i3 = "MF2")
    ),
    "metMatch must include one metafounder assignment for each row of genotypicMatrix"
  )
})

test_that("getGamma keeps groups as a legacy alias for Garcia-Baccino metMatch", {
  genotypicMatrix <- matrix(c(
    0, 2,
    2, 2,
    0, 0,
    2, 0
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("i", 1:4), paste0("m", 1:2)))
  groups <- c(i1 = "MF1", i2 = "MF1", i3 = "MF2", i4 = "MF2")

  fit_legacy <- getGamma(
    genotypicMatrix = genotypicMatrix,
    groups = groups,
    method = "garcia-baccino"
  )
  fit_explicit <- getGamma(
    genotypicMatrix = genotypicMatrix,
    metMatch = groups,
    method = "garcia-baccino"
  )

  expect_equal(fit_legacy$gamma, fit_explicit$gamma, tolerance = 1e-8)
  expect_equal(fit_legacy$pCentered, fit_explicit$pCentered, tolerance = 1e-8)
})

test_that("Garcia-Baccino keeps backward-compatible aliases for renamed outputs", {
  genotypicMatrix <- matrix(c(
    0, 2,
    2, 2,
    0, 0,
    2, 0
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("i", 1:4), paste0("m", 1:2)))
  metMatch <- c(i1 = "MF1", i2 = "MF1", i3 = "MF2", i4 = "MF2")

  fit <- getGammaGarciaBaccino(
    genotypicMatrix = genotypicMatrix,
    metMatch = metMatch
  )

  expect_equal(fit$alleleFrequencies, fit$frequencies, tolerance = 1e-8)
  expect_equal(fit$pCentered, fit$alleFrequency, tolerance = 1e-8)
})
