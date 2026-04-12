# Gene flow from a pedigree with and without metafounders

test_that("gene_flow", {

  ped <- pedigree(
    sire  = c(NA, NA, 1, 1, 4, 5),
    dam   = c(NA, NA, 2, NA, 3, 2),
    label = 1:6
  )
  
  TEsp <- getT(ped)
  TExp <- matrix(data = 
    c(
    1.00, 0.000, 0.00, 0.00, 0.0, 0,
    0.00, 1.000, 0.00, 0.00, 0.0, 0,
    0.50, 0.500, 1.00, 0.00, 0.0, 0,
    0.50, 0.000, 0.00, 1.00, 0.0, 0,
    0.50, 0.250, 0.50, 0.50, 1.0, 0,
    0.25, 0.625, 0.25, 0.25, 0.5, 1
  ), byrow = TRUE, nrow = 6)
  
  TExp <- as(TExp, "dtCMatrix")
  dimnames(TExp) <- list(ped@label, ped@label)
  
  
  expect_equal(TEsp, TExp, tolerance = .Machine$double.eps)

  # with metafounders
  ped_multi <- pedigree(
    sire  = c(0, 0, 1, 2, 3, 3, 5, 5),
    dam   = c(0, 0, 1, 2, 4, 4, 6, 4),
    label = 1:8
  )

  TEsp <- getT(ped_multi)
  TExp <- matrix(data = 
    c(
    1.00, 0.000, 0.00, 0.00, 0.0, 0.0, 0.0, 0.0,
    0.00, 1.000, 0.00, 0.00, 0.0, 0.0, 0.0, 0.0,
    1.00, 0.000, 1.00, 0.00, 0.0, 0.0, 0.0, 0.0,
    0.00, 1.000, 0.00, 1.00, 0.0, 0.0, 0.0, 0.0,
    0.50, 0.500, 0.50, 0.50, 1.0, 0.0, 0.0, 0.0,
    0.50, 0.500, 0.50, 0.50, 0.0, 1.0, 0.0, 0.0,
    0.50, 0.500, 0.50, 0.50, 0.5, 0.5, 1.0, 0.0,
    0.25, 0.750, 0.25, 0.75, 0.5, 0.0, 0.0, 1.0
  ), byrow = TRUE, nrow = 8)

  TExp <- as(TExp, "dtCMatrix")
  dimnames(TExp) <- list(ped_multi@label, ped_multi@label)
  
  expect_equal(TEsp, TExp, tolerance = .Machine$double.eps)

})



# Gene flow from a pedigree with and without metafounders
test_that("gene_flow_inv", {

  ped <- pedigree(
    sire  = c(NA, NA, 1, 1, 4, 5),
    dam   = c(NA, NA, 2, NA, 3, 2),
    label = 1:6
  )

  TInvEsp <- getTInv(ped)
  TInvExp <- matrix(
    c(
      1.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  1.0,  0.0,  0.0,  0.0,  0.0,
     -0.5, -0.5,  1.0,  0.0,  0.0,  0.0,
     -0.5,  0.0,  0.0,  1.0,  0.0,  0.0,
      0.0,  0.0, -0.5, -0.5,  1.0,  0.0,
      0.0, -0.5,  0.0,  0.0, -0.5,  1.0
    ),
    byrow = TRUE,
    nrow = 6
  )
  dimnames(TInvExp) <- list(as.character(ped@label), as.character(ped@label))

  expect_equal(
    as.matrix(TInvEsp),
    TInvExp,
    tolerance = .Machine$double.eps
  )

  ped_multi <- pedigree(
    sire  = c(0, 0, 1, 2, 3, 3, 5, 5),
    dam   = c(0, 0, 1, 2, 4, 4, 6, 4),
    label = 1:8
  )

  TInvEsp <- getTInv(ped_multi)
  TInvExp <- matrix(
    c(
      1.00,  0.00,  0.00,  0.00,  0.0,  0.0,  0.0,  0.0,
      0.00,  1.00,  0.00,  0.00,  0.0,  0.0,  0.0,  0.0,
     -1.00,  0.00,  1.00,  0.00,  0.0,  0.0,  0.0,  0.0,
      0.00, -1.00,  0.00,  1.00,  0.0,  0.0,  0.0,  0.0,
      0.00,  0.00, -0.50, -0.50,  1.0,  0.0,  0.0,  0.0,
      0.00,  0.00, -0.50, -0.50,  0.0,  1.0,  0.0,  0.0,
      0.00,  0.00,  0.00,  0.00, -0.5, -0.5,  1.0,  0.0,
      0.00,  0.00,  0.00, -0.50, -0.5,  0.0,  0.0,  1.0
    ),
    byrow = TRUE,
    nrow = 8
  )
  dimnames(TInvExp) <- list(as.character(ped_multi@label), as.character(ped_multi@label))

  expect_equal(
    as.matrix(TInvEsp),
    TInvExp,
    tolerance = .Machine$double.eps
  )
})
