rosa <- function(Li, z) {
    Matrix::crossprod(Li, Li %*% z)
}

is_big_genotype_matrix <- function(x) {
    inherits(x, "big.matrix")
}

prepare_genotype_matrix <- function(genotypicMatrix, arg = "genotypicMatrix") {
    is_bigmat <- is_big_genotype_matrix(genotypicMatrix)
    if (!is.matrix(genotypicMatrix) && !is_bigmat) {
        stop(arg, " must be a matrix or big.matrix.")
    }

    if (!is_bigmat) {
        genotypicMatrix <- as.matrix(genotypicMatrix)
        storage.mode(genotypicMatrix) <- "double"
    }

    markerCount <- ncol(genotypicMatrix)
    if (markerCount < 1L) {
        stop(arg, " must contain at least one marker column.")
    }

    structure(
        list(
            source = genotypicMatrix,
            is_bigmat = is_bigmat,
            individualCount = nrow(genotypicMatrix),
            markerCount = markerCount,
            block_size = if (is_bigmat) min(1024L, markerCount) else markerCount
        ),
        class = "prepared_genotype_matrix"
    )
}

genotype_block_apply <- function(genotypes, FUN) {
    if (!inherits(genotypes, "prepared_genotype_matrix")) {
        stop("genotypes must come from prepare_genotype_matrix().")
    }

    markerCount <- genotypes$markerCount
    if (markerCount < 1L) {
        return(invisible(NULL))
    }

    if (!genotypes$is_bigmat) {
        FUN(genotypes$source, seq_len(markerCount))
        return(invisible(NULL))
    }

    for (start in seq.int(1L, markerCount, by = genotypes$block_size)) {
        stop_col <- min(start + genotypes$block_size - 1L, markerCount)
        cols <- start:stop_col
        block <- as.matrix(genotypes$source[, cols, drop = FALSE])
        storage.mode(block) <- "double"
        FUN(block, cols)
    }

    invisible(NULL)
}


validate_gamma_matrix <- function(gamma, n_groups = NULL, arg = "gamma") {
    if (!is.matrix(gamma) || !is.numeric(gamma)) {
        stop(arg, " must be a numeric matrix.")
    }

    if (nrow(gamma) != ncol(gamma)) {
        stop(arg, " must be a square matrix.")
    }

    if (anyNA(gamma) || any(!is.finite(gamma))) {
        stop(arg, " must contain only finite numeric values.")
    }

    if (!is.null(n_groups) && nrow(gamma) != n_groups) {
        stop(arg, " must be a square matrix with dimension ", n_groups, " x ", n_groups, ".")
    }

    if (!isSymmetric(gamma)) {
        stop(arg, " must be a symmetric matrix.")
    }

    chol_ok <- tryCatch({
        chol(gamma)
        TRUE
    }, error = function(e) FALSE)

    if (!chol_ok) {
        stop(arg, " must be a symmetric positive-definite matrix.")
    }

    gamma
}
standardize_ids <- function(ids, all_ids, arg = "genotyped") {
    ids <- as.character(ids)

    if (length(ids) < 1L) {
        stop(arg, " must contain at least one ID.")
    }

    if (!all(ids %in% all_ids)) {
        missing_ids <- ids[!ids %in% all_ids]
        stop("Some ", arg, " IDs are not present: ", paste(missing_ids, collapse = ", "))
    }

    ids
}

row_ids <- function(x, arg = "x", required = TRUE) {
    ids <- rownames(x)

    if (is.null(ids)) {
        if (required) {
            stop("rownames(", arg, ") must be provided.")
        }
        return(NULL)
    }

    ids <- as.character(ids)
    if (anyNA(ids) || any(ids == "")) {
        stop("rownames(", arg, ") must not contain NA or empty values.")
    }
    if (anyDuplicated(ids)) {
        stop("rownames(", arg, ") must be unique.")
    }

    ids
}

align_ids <- function(x, ids, arg = "x", ids_arg = "ids",
                      item_label = "value", target_label = ids_arg) {
    if (is.null(x)) {
        stop(arg, " must be provided.")
    }
    if (length(x) < 1L) {
        stop(arg, " must contain at least one ", item_label, ".")
    }

    ids <- as.character(ids)
    if (length(ids) < 1L) {
        stop(ids_arg, " must contain at least one ID.")
    }
    if (anyNA(ids) || any(ids == "")) {
        stop(ids_arg, " must not contain NA or empty values.")
    }
    if (anyDuplicated(ids)) {
        stop(ids_arg, " must be unique.")
    }

    x_ids <- names(x)
    if (is.null(x_ids) || anyNA(x_ids) || any(x_ids == "")) {
        stop(arg, " must be a named vector whose names match ", ids_arg, ".")
    }
    if (anyDuplicated(x_ids)) {
        stop("names(", arg, ") must be unique.")
    }

    ids_not_in_target <- setdiff(x_ids, ids)
    if (length(ids_not_in_target) > 0L) {
        stop(
            "All IDs in ", arg, " must be present in ", ids_arg, ". Unknown IDs: ",
            paste(ids_not_in_target, collapse = ", "),
            "."
        )
    }

    ids_without_value <- setdiff(ids, x_ids)
    if (length(ids_without_value) > 0L) {
        stop(
            arg,
            " must include one ", item_label, " for each ", target_label, ". Missing IDs: ",
            paste(ids_without_value, collapse = ", "),
            "."
        )
    }

    x <- x[ids]
    x <- as.character(x)

    if (anyNA(x) || any(x == "")) {
        stop(arg, " must not contain NA or empty ", paste0(item_label, "s"), ".")
    }

    x
}

infer_genotype_matrix_coding <- function(genotypicMatrix, arg = "genotypicMatrix") {
    values <- genotypicMatrix[!is.na(genotypicMatrix)]

    if (any(!(values %in% c(-1, 0, 1, 2)))) {
        stop(arg, " must be coded as 0, 1, 2 or centered -1, 0, 1, with optional NA values.")
    }

    has_neg1 <- any(values == -1)
    has_two <- any(values == 2)

    if (has_neg1 && has_two) {
        stop(arg, " cannot mix centered (-1, 0, 1) and uncentered (0, 1, 2) codings.")
    }

    if (has_neg1) {
        return("centered")
    }

    "012"
} # Test si ingresa una matriz centrada -1, 0 , 1 o no centrada 0, 1, 2

resolve_genotype_matrix_coding <- function(genotypicMatrix, arg = "genotypicMatrix") {
    if (!inherits(genotypicMatrix, "prepared_genotype_matrix")) {
        stop("genotypicMatrix must come from prepare_genotype_matrix().")
    }

    if (!genotypicMatrix$is_bigmat) {
        return(infer_genotype_matrix_coding(genotypicMatrix$source, arg = arg))
    }

    has_neg1 <- FALSE
    has_two <- FALSE

    genotype_block_apply(genotypicMatrix, function(block, cols) {
        values <- block[!is.na(block)]

        if (any(!(values %in% c(-1, 0, 1, 2)))) {
            stop(arg, " must be coded as 0, 1, 2 or centered -1, 0, 1, with optional NA values.")
        }

        if (any(values == -1)) {
            has_neg1 <<- TRUE
        }
        if (any(values == 2)) {
            has_two <<- TRUE
        }

        if (has_neg1 && has_two) {
            stop(arg, " cannot mix centered (-1, 0, 1) and uncentered (0, 1, 2) codings.")
        }
    })

    if (has_neg1) {
        return("centered")
    }

    "012"
}

center_genotype_block <- function(block, coding) {
    storage.mode(block) <- "double"

    if (coding == "centered") {
        block[is.na(block)] <- 0
        return(block)
    }

    block[is.na(block)] <- 1
    block - 1
}

z_marker_count <- function(Z) {
    if (inherits(Z, "legarra_Z_bigmemory")) {
        return(Z$source$markerCount)
    }

    ncol(Z)
}

z_block_apply <- function(Z, FUN) {
    if (inherits(Z, "legarra_Z_bigmemory")) {
        genotype_block_apply(Z$source, function(block, cols) {
            Z_block <- center_genotype_block(block, coding = Z$coding)
            rownames(Z_block) <- Z$genotyped
            FUN(Z_block, cols)
        })
        return(invisible(NULL))
    }

    FUN(Z, seq_len(ncol(Z)))
    invisible(NULL)
}

# Garcia-Baccino method --------------------------------------------------------

garcia_baccino_naive <- function(genotypicMatrix, metMatch) {
    genotypes <- prepare_genotype_matrix(genotypicMatrix)
    genotypicMatrix <- genotypes$source

    if (genotypes$individualCount < 1L) {
        stop("genotypicMatrix must contain at least one individual row.")
    }

    m <- genotypes$markerCount

    matrix_ids <- row_ids(genotypicMatrix, arg = "genotypicMatrix")
    metMatch <- align_ids(
        metMatch,
        ids = matrix_ids,
        arg = "metMatch",
        ids_arg = "rownames(genotypicMatrix)",
        item_label = "metafounder assignment",
        target_label = "row of genotypicMatrix"
    )
    metafounders <- unique(metMatch)

    coding <- resolve_genotype_matrix_coding(genotypes)

    # For each metafounder, estimate marker allele frequencies by averaging the
    # genotypes of the individuals assigned to that metafounder.
    alleleFrequencies <- t(vapply(metafounders, function(metafounder) {
        idx <- which(metMatch == metafounder)
        avg <- colMeans(genotypicMatrix[idx, , drop = FALSE], na.rm = TRUE)

        if (coding == "centered") {
            return((avg + 1) / 2)
        }

        avg / 2
    }, numeric(m)))

    if (any(!is.finite(alleleFrequencies))) {
        stop("Non-finite metafounder frequencies detected. Check for empty groups or markers with only missing values.")
    }

    rownames(alleleFrequencies) <- metafounders
    colnames(alleleFrequencies) <- colnames(genotypicMatrix)

    # Garcia-Baccino works with allele frequencies centered at the conceptual
    # base p = 0.5. This quantity is analogous to the metafounder "Z" term,
    # but on the allele-frequency scale: pCentered = p - 0.5 = (2p - 1) / 2.
    pCentered <- alleleFrequencies - 0.5
    Gamma <- (8 / m) * tcrossprod(pCentered)

    list(
        gamma = Gamma,
        alleleFrequencies = alleleFrequencies,
        pCentered = pCentered,
        frequencies = alleleFrequencies,
        alleFrequency = pCentered,
        metafounders = metafounders,
        estimator = "naive"
    )
}
prepare_legarra_data <- function(genotypicMatrix, ped, genotypedIds = NULL,
                                 idSource = c("auto", "rownames", "genotyped")) {
    stopifnot(is(ped, "pedigree"))
    idSource <- match.arg(idSource)
    genotypes <- prepare_genotype_matrix(genotypicMatrix)
    genotypicMatrix <- genotypes$source

    ped_ids <- ped@label
    matrixIds <- row_ids(genotypicMatrix, arg = "genotypicMatrix", required = FALSE)
    hasMatrixIds <- !is.null(matrixIds)

    validMatrixIds <- hasMatrixIds && all(matrixIds %in% ped_ids)

    if (idSource == "rownames") {
        if (!validMatrixIds) {
            stop("rownames(genotypicMatrix) must match ped@label when idSource = 'rownames'.")
        }
        genotypedIds <- matrixIds
        usedSource <- "rownames"
    } else if (idSource == "genotyped") {
        if (is.null(genotypedIds)) {
            stop("genotypedIds must be supplied when idSource = 'genotyped'.")
        }
        genotypedIds <- standardize_ids(genotypedIds, ped_ids, arg = "genotypedIds")
        if (nrow(genotypicMatrix) != length(genotypedIds)) {
            stop("nrow(genotypicMatrix) must equal length(genotypedIds) when idSource = 'genotyped'.")
        }
        usedSource <- "genotyped"
    } else if (validMatrixIds) {
        genotypedIds <- matrixIds
        usedSource <- "rownames"
    } else {
        if (is.null(genotypedIds)) {
            stop(
                "For method = 'legarra', provide rownames(genotypicMatrix) matching ped@label ",
                "or supply genotypedIds and set idSource = 'genotyped'."
            )
        }
        genotypedIds <- standardize_ids(genotypedIds, ped_ids, arg = "genotypedIds")
        if (nrow(genotypicMatrix) != length(genotypedIds)) {
            stop("nrow(genotypicMatrix) must equal length(genotypedIds) when rownames(genotypicMatrix) are not used.")
        }
        usedSource <- "genotyped"
    }

    coding <- resolve_genotype_matrix_coding(genotypes)
    if (genotypes$is_bigmat) {
        Z <- structure(
            list(
                source = genotypes,
                coding = coding,
                genotyped = genotypedIds
            ),
            class = "legarra_Z_bigmemory"
        )
    } else if (coding == "centered") {
        genotypicMatrix[is.na(genotypicMatrix)] <- 0
        Z <- genotypicMatrix
    } else {
        genotypicMatrix[is.na(genotypicMatrix)] <- 1
        Z <- genotypicMatrix - 1
    }

    if (!genotypes$is_bigmat) {
        rownames(Z) <- genotypedIds
    }

    list(
        Z = Z,
        genotyped = genotypedIds,
        id_source = usedSource
    )
}

logdet_block_A22_from_Li <- function(Li, genotypedIds) {
    if (is.null(rownames(Li)) || is.null(colnames(Li))) {
        stop("Li must have rownames and colnames corresponding to individual IDs.")
    }

    all_ids <- rownames(Li)
    genotypedIds <- standardize_ids(genotypedIds, all_ids, arg = "genotypedIds")
    idx.non.geno <- setdiff(all_ids, genotypedIds)

    Li_diag <- as.numeric(Matrix::diag(Li))
    if (any(!is.finite(Li_diag)) || any(Li_diag <= 0)) {
        stop("Li must have a finite positive diagonal.")
    }
    logdet_Ainv <- 2 * sum(log(Li_diag))

    if (length(idx.non.geno) == 0L) {
        return(-logdet_Ainv)
    }

    Li11 <- Li[, idx.non.geno, drop = FALSE]
    A11 <- Matrix::crossprod(Li11)
    det_A11 <- determinant(A11, logarithm = TRUE)
    if (det_A11$sign <= 0) {
        stop("The non-genotyped precision block is not positive definite.")
    }

    -logdet_Ainv + as.numeric(det_A11$modulus)
}

masuda_implicit <- function(q, Li, genotypedIds) {
    if (is.null(rownames(Li)) || is.null(colnames(Li))) {
        stop("Li must have rownames and colnames corresponding to individual IDs.")
    }

    all_ids <- rownames(Li)
    genotypedIds <- standardize_ids(genotypedIds, all_ids, arg = "genotypedIds")
    idx.non.geno <- setdiff(all_ids, genotypedIds)

    if (is.null(dim(q))) {
        q <- matrix(as.vector(q), ncol = 1)
    } else {
        q <- as.matrix(q)
        storage.mode(q) <- "double"
    }

    if (nrow(q) != length(genotypedIds)) {
        stop("q must have length(genotypedIds) rows.")
    }

    subMatrix <- matrix(0, nrow = length(all_ids), ncol = ncol(q),
                        dimnames = list(all_ids, colnames(q)))

    subMatrix[genotypedIds, ] <- q
    y_full <- rosa(Li, subMatrix)
    y <- y_full[genotypedIds, , drop = FALSE]

    if (length(idx.non.geno) == 0L) {
        return(y)
    }

    v <- y_full[idx.non.geno, , drop = FALSE]

    Li11 <- Li[, idx.non.geno, drop = FALSE]
    A11 <- Matrix::crossprod(Li11)
    w <- as.matrix(solve(A11, v))

    subMatrix[,] <- 0
    subMatrix[idx.non.geno, ] <- w
    x_full <- rosa(Li, subMatrix)
    x <- x_full[genotypedIds, , drop = FALSE]

    y - x
}

legarra_gamma_components <- function(ped, gamma, genotypedIds, Z,
                                     method = c("byMarker", "byIndividual")) {
    stopifnot(is(ped, "pedigree"))
    method <- match.arg(method)

    gamma <- validate_gamma_matrix(gamma, arg = "gamma")
    all_ids <- ped@label
    genotypedIds <- standardize_ids(genotypedIds, all_ids, arg = "genotypedIds")
    markerCount <- z_marker_count(Z)

    Li <- getLInv(ped, gamma = gamma)
    logdetA22 <- logdet_block_A22_from_Li(Li, genotypedIds)
    determinant_term <- -(markerCount / 2) * logdetA22

    quad_sum <- 0
    z_block_apply(Z, function(Z_block, cols) {
        A22inv_Z_block <- masuda_implicit(Z_block, Li, genotypedIds)
        quad_sum <<- quad_sum + sum(Z_block * A22inv_Z_block)
    })

    if (method == "byMarker") {
        trace_term <- NULL
        quadratic_term <- -quad_sum
    } else {
        trace_term <- (2 / markerCount) * quad_sum
        quadratic_term <- -(markerCount / 2) * trace_term
    }

    list(
        determinant_term = determinant_term,
        trace_term = trace_term,
        quadratic_term = quadratic_term,
        logLik = determinant_term + quadratic_term
    )
}

legarra_single_exact_ml_components <- function(ped, genotypedIds, Z,
                                               zero_gamma = 1e-12) {
    stopifnot(is(ped, "pedigree"))
    genotypedIds <- standardize_ids(genotypedIds, ped@label, arg = "genotypedIds")

    if (!is.numeric(zero_gamma) || length(zero_gamma) != 1L || !is.finite(zero_gamma) ||
        zero_gamma <= 0) {
        stop("zero_gamma must be a single positive number.")
    }

    n <- length(genotypedIds)
    k <- z_marker_count(Z)
    Li0 <- getLInv(ped, gamma = matrix(zero_gamma, nrow = 1L))
    ones_vec <- matrix(1, nrow = n, ncol = 1)
    A22inv_1 <- masuda_implicit(ones_vec, Li0, genotypedIds)
    a <- sum(A22inv_1)
    qprime <- as.numeric(A22inv_1[, 1])
    logdetA22 <- logdet_block_A22_from_Li(Li0, genotypedIds)

    quad_sum <- 0
    t_sum_sq <- 0
    z_block_apply(Z, function(Z_block, cols) {
        A22inv_Z_block <- masuda_implicit(Z_block, Li0, genotypedIds)
        quad_sum <<- quad_sum + sum(Z_block * A22inv_Z_block)
        t_block <- as.numeric(crossprod(qprime, Z_block))
        t_sum_sq <<- t_sum_sq + sum(t_block^2)
    })

    b <- (2 / k) * quad_sum
    c <- (2 / k) * t_sum_sq

    e3 <- -n * (-1/2 + a)^2 / 2
    e2 <- n * (((-3/2 + a) + a - b * (-1/2 + a) + c) * (-1/2 + a))
    e1 <- (n - 1) * (-3/2 + 2 * a) - (-1 + 2 * a) * (-3/2 + a + b)
    e0 <- n - 2 * a - b + 2 * c

    list(
        n = n,
        k = k,
        a = a,
        b = b,
        c = c,
        logdetA22 = logdetA22,
        coefficients = c(e0, e1, e2, e3)
    )
}

legarra_single_real_roots <- function(fit, lower = 1e-6, upper = 1) {
    roots <- polyroot(fit$coefficients)
    real_roots <- sort(Re(roots[abs(Im(roots)) < 1e-8]))
    admissible_roots <- real_roots[real_roots >= lower & real_roots <= upper]

    list(
        real_roots = real_roots,
        admissible_roots = admissible_roots
    )
}

legarra_single_candidates <- function(fit, roots, lower = 1e-6, upper = 1) {
    gamma_values <- sort(unique(c(lower, roots$admissible_roots, upper)))
    logLik_values <- vapply(gamma_values, legarra_single_logLik, numeric(1), fit = fit)

    candidates <- data.frame(
        gamma = gamma_values,
        logLik = logLik_values,
        is_root = gamma_values %in% roots$admissible_roots,
        is_boundary = gamma_values %in% c(lower, upper)
    )
    candidates <- candidates[order(candidates$logLik, decreasing = TRUE), , drop = FALSE]
    rownames(candidates) <- NULL
    candidates
}

print.legarra_single_gamma_candidates <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    recommended <- x$candidates[1L, , drop = FALSE]

    cat("Choose the gamma with the highest logLik.\n")
    cat("Recommended gamma:", format(recommended$gamma, digits = digits), "\n")
    cat("Recommended logLik:", format(recommended$logLik, digits = digits), "\n")
    cat("\nCandidates ordered by logLik:\n")
    print(x$candidates, row.names = FALSE, digits = digits, ...)

    invisible(x)
}

legarra_single_exact_ml <- function(ped, genotypedIds, Z,
                                    lower = 1e-6, upper = 1,
                                    zero_gamma = 1e-12) {
    if (!is.numeric(lower) || length(lower) != 1L || !is.finite(lower) ||
        lower < 0 || lower >= upper) {
        stop("lower must be a single non-negative number smaller than upper.")
    }
    if (!is.numeric(upper) || length(upper) != 1L || !is.finite(upper) ||
        upper <= 0 || upper <= lower) {
        stop("upper must be a single positive number greater than lower.")
    }

    fit <- legarra_single_exact_ml_components(
        ped = ped,
        genotypedIds = genotypedIds,
        Z = Z,
        zero_gamma = zero_gamma
    )
    roots <- legarra_single_real_roots(
        fit = fit,
        lower = lower,
        upper = upper
    )
    candidates <- legarra_single_candidates(
        fit = fit,
        roots = roots,
        lower = lower,
        upper = upper
    )

    structure(
        list(
            real_roots = roots$real_roots,
            admissible_roots = roots$admissible_roots,
            candidates = candidates
        ),
        class = "legarra_single_gamma_candidates"
    )
}

legarra_single_logLik <- function(gamma_value, fit) {
    denom1 <- 1 - gamma_value / 2
    denom2 <- 1 - gamma_value / 2 + gamma_value * fit$a

    -(fit$k / 2) * fit$logdetA22 -
        (fit$k / 2) * (fit$n - 1) * log(denom1) -
        (fit$k / 2) * log(denom2) -
        (fit$k / 2) * fit$b / denom1 +
        (fit$k / 2) * fit$c * gamma_value / (denom1 * denom2)
}

select_legarra_single_gamma <- function(ped, genotypedIds, Z,
                                        lower = 1e-6, upper = 1,
                                        zero_gamma = 1e-12) {
    fit <- legarra_single_exact_ml_components(
        ped = ped,
        genotypedIds = genotypedIds,
        Z = Z,
        zero_gamma = zero_gamma
    )
    roots <- legarra_single_real_roots(
        fit = fit,
        lower = lower,
        upper = upper
    )
    candidates <- legarra_single_candidates(
        fit = fit,
        roots = roots,
        lower = lower,
        upper = upper
    )

    list(
        gamma = candidates$gamma[1],
        logLik = candidates$logLik[1],
        real_roots = roots$real_roots,
        admissible_roots = roots$admissible_roots,
        candidates = candidates
    )
}

pseudo_em_gamma_update <- function(ped, gamma, genotypedIds, Z) {
    stopifnot(is(ped, "pedigree"))
    gamma <- validate_gamma_matrix(gamma, arg = "gamma")
    genotypedIds <- standardize_ids(genotypedIds, ped@label, arg = "genotypedIds")
    metafounderNames <- rownames(gamma)
    if (is.null(metafounderNames)) {
        metafounderNames <- ped@label[seq_len(nrow(gamma))]
    }

    Tmat <- getT(ped)
    Q2 <- as.matrix(Tmat[genotypedIds, metafounderNames, drop = FALSE])
    storage.mode(Q2) <- "double"

    if (nrow(Q2) != length(genotypedIds) || ncol(Q2) != nrow(gamma)) {
        stop("Q2 must have length(genotypedIds) rows and nrow(gamma) columns.")
    }

    Li <- getLInv(ped, gamma = gamma)
    markerCount <- z_marker_count(Z)
    Qgamma <- Q2 %*% gamma
    mf_precision <- as.matrix(Matrix::crossprod(Li[, metafounderNames, drop = FALSE]))
    base_term <- solve(mf_precision)
    sum_Fhat_crossprod <- matrix(0, nrow = nrow(gamma), ncol = ncol(gamma))
    z_block_apply(Z, function(Z_block, cols) {
        A22inv_Z_block <- masuda_implicit(Z_block, Li, genotypedIds)
        Fhat_block <- crossprod(Qgamma, A22inv_Z_block)
        sum_Fhat_crossprod <<- sum_Fhat_crossprod + tcrossprod(Fhat_block)
    })
    gamma_next <- base_term + (2 / markerCount) * sum_Fhat_crossprod
    gamma_next <- as.matrix(Matrix::forceSymmetric(gamma_next, uplo = "U"))

    dimnames(base_term) <- dimnames(gamma)
    dimnames(mf_precision) <- dimnames(gamma)
    dimnames(gamma_next) <- dimnames(gamma)

    list(
        gamma_next = gamma_next
    )
}

logLikeGamma <- function(ped, gamma, genotypedIds, genotypicMatrix,
                         method = c("byMarker", "byIndividual")) {
    stopifnot(is(ped, "pedigree"))
    method <- match.arg(method)

    legarra <- prepare_legarra_data(
        genotypicMatrix = genotypicMatrix,
        ped = ped,
        genotypedIds = genotypedIds,
        idSource = "genotyped"
    )

    legarra_gamma_components(
        ped = ped,
        gamma = gamma,
        genotypedIds = legarra$genotyped,
        Z = legarra$Z,
        method = method
    )$logLik
}
