#' Estimate a Metafounder Gamma Matrix by the Legarra Method
#'
#' Estimate \code{gamma} for one or more metafounders from pedigree and marker
#' data following the likelihood-based implementation of Legarra.
#'
#' @param ped A \code{\link{pedigree}} object.
#' @param genotypicMatrix Genotype matrix with individuals in rows and markers
#'   in columns. The standard coding is \code{0/1/2}. Centered coding
#'   \code{-1/0/1} is also accepted when it can be identified unambiguously.
#'   If a matrix contains only \code{0/1} values, it is interpreted as
#'   standard \code{0/1/2} input.
#' @param groups Single positive integer giving the number of metafounders.
#' @param genotypedIds Character or integer IDs for the genotyped individuals.
#'   This is only required when the rows of \code{genotypicMatrix} are not
#'   already identified by pedigree labels, or when
#'   \code{idSource = "genotyped"} is used.
#' @param likelihood One of \code{"byMarker"} or \code{"byIndividual"} for the
#'   exact Legarra log-likelihood. For \code{groups = 1}, it is the criterion
#'   maximized by the scalar optimization. For \code{groups > 1}, the
#'   pseudo-EM update does not depend on this choice. The exact log-likelihood
#'   is only computed if requested by the user, either through
#'   \code{logLikeGamma()} or via \code{control$compute_logLik = TRUE}.
#' @param idSource Controls where the IDs for the rows of
#'   \code{genotypicMatrix} come from. Use \code{"rownames"} to read them from
#'   \code{rownames(genotypicMatrix)}, use \code{"genotyped"} to read them from
#'   \code{genotypedIds} in row order, and use \code{"auto"} to prefer
#'   \code{rownames(genotypicMatrix)} when they match \code{ped@label} and
#'   otherwise fall back to \code{genotypedIds}. If both are available and
#'   \code{idSource = "auto"}, the row names take precedence.
#' @param control Optional named list used only when \code{groups > 1}.
#'   Supported entries are \code{maxit}, \code{tol}, \code{keep_history},
#'   and \code{compute_logLik}. Defaults are
#'   \code{list(maxit = 200, tol = 1e-8, keep_history = FALSE, compute_logLik = FALSE)}.
#'   For example,
#'   \code{control = list(maxit = 500, tol = 1e-10, compute_logLik = TRUE)}.
#' @return A list containing the estimated \code{gamma} matrix and
#'   method-specific metadata.
#' @details If \code{groups = 1}, \code{gamma} is a \eqn{1 \times 1} matrix and
#'   the exact single-metafounder cubic equation from Legarra is used. If
#'   \code{groups > 1}, the estimate is obtained by the pseudo-EM update
#'   described by Legarra (2024), while the exact log-likelihood remains
#'   available through \code{logLikeGamma()} for verification. The
#'   \code{likelihood} argument chooses only how that exact log-likelihood is
#'   evaluated: by marker or by individual. The pseudo-EM iterations themselves
#'   use only their own update and convergence criterion; exact log-likelihood
#'   checks are optional because they can be costly for large matrices.
#'   Genotypes coded as \code{0/1/2} are internally recoded to \code{Z = M - 1};
#'   centered input coded as \code{-1/0/1} is used directly when detected
#'   unambiguously.
#' @export
getGammaLegarra <- function(ped, genotypicMatrix, groups, genotypedIds = NULL,
                            likelihood = c("byMarker", "byIndividual"),
                            idSource = c("auto", "rownames", "genotyped"),
                            control = NULL) {
    if (is.null(ped)) {
        stop("ped must be provided for getGammaLegarra().")
    }
    if (is.null(genotypicMatrix)) {
        stop("genotypicMatrix must be provided for method = 'legarra'.")
    }

    stopifnot(is(ped, "pedigree"))

    likelihood <- match.arg(likelihood)
    idSource <- match.arg(idSource)

    legarra <- prepare_legarra_data(
        genotypicMatrix = genotypicMatrix,
        ped = ped,
        genotypedIds = genotypedIds,
        idSource = idSource
    )
    genotypedIds <- legarra$genotyped
    Z <- legarra$Z

    if (is.null(groups) || !is.numeric(groups) || length(groups) != 1L ||
        is.na(groups) || groups < 1) {
        stop("groups must be a single positive integer for method = 'legarra'.")
    }
    n_groups <- as.integer(groups)
    metafounderNames <- ped@label[seq_len(n_groups)]

    if (n_groups == 1L) {
        interval <- c(1e-6, 1)
        single_fit <- select_legarra_single_gamma(
            ped = ped,
            genotypedIds = genotypedIds,
            Z = Z,
            lower = interval[1],
            upper = interval[2]
        )
        gamma <- matrix(single_fit$gamma, nrow = 1L, dimnames = list(metafounderNames, metafounderNames))
        return(list(
            gamma = gamma,
            logLik = single_fit$logLik,
            genotyped = genotypedIds,
            id_source = legarra$id_source,
            groups = 1L,
            method = "legarra",
            likelihood = likelihood,
            optimizer = "closed-form",
            interval = interval
        ))
    }

    defaults <- list(
        maxit = 200L,
        tol = 1e-8,
        keep_history = FALSE,
        compute_logLik = FALSE
    )
    if (is.null(control)) {
        control <- defaults
    } else {
        if (!is.list(control)) {
            stop("control must be NULL or a named list for multi-metafounder pseudo-EM.")
        }
        if (length(control) == 0L) {
            control <- defaults
        } else {
            nm <- names(control)
            if (is.null(nm) || anyNA(nm) || any(!nzchar(nm))) {
                stop("control must be a named list for multi-metafounder pseudo-EM.")
            }

            unknown <- setdiff(nm, names(defaults))
            if (length(unknown) > 0L) {
                stop(
                    "Unknown control entries for multi-metafounder pseudo-EM: ",
                    paste(sort(unknown), collapse = ", "),
                    "."
                )
            }

            defaults[nm] <- control
            control <- defaults
        }
    }

    maxit <- as.integer(control$maxit)
    tol <- as.numeric(control$tol)
    keep_history <- isTRUE(control$keep_history)
    compute_logLik <- isTRUE(control$compute_logLik)
    if (!is.finite(maxit) || maxit < 1L) {
        stop("control$maxit must be a positive integer for multi-metafounder pseudo-EM.")
    }
    if (!is.finite(tol) || tol <= 0) {
        stop("control$tol must be a single positive number for multi-metafounder pseudo-EM.")
    }

    start_gamma <- diag(0.1, n_groups)
    gamma <- validate_gamma_matrix(start_gamma, n_groups = n_groups, arg = "start")
    dimnames(gamma) <- list(metafounderNames, metafounderNames)

    chol_prev <- chol(gamma)
    criterion <- Inf
    converged <- FALSE
    history <- if (keep_history) vector("list", maxit) else NULL
    iter <- 0L

    for (step in seq_len(maxit)) {
        update <- pseudo_em_gamma_update(
            ped = ped,
            gamma = gamma,
            genotypedIds = genotypedIds,
            Z = Z
        )

        gamma_next <- update$gamma_next
        dimnames(gamma_next) <- list(metafounderNames, metafounderNames)
        chol_next <- tryCatch(chol(gamma_next), error = function(e) NULL)

        if (is.null(chol_next)) {
            stop("Pseudo-EM update produced a non-positive-definite gamma matrix.")
        }

        criterion <- sum((chol_next - chol_prev)^2) / sum(chol_prev^2)
        gamma <- gamma_next
        chol_prev <- chol_next
        iter <- step

        if (keep_history) {
            history_row <- list(iteration = step, criterion = criterion)
            if (compute_logLik) {
                history_row$logLik <- legarra_gamma_components(
                    ped = ped,
                    gamma = gamma,
                    genotypedIds = genotypedIds,
                    Z = Z,
                    method = likelihood
                )$logLik
            }
            history[[step]] <- history_row
        }

        if (criterion < tol) {
            converged <- TRUE
            break
        }
    }

    final_logLik <- NULL
    if (compute_logLik) {
        final_logLik <- legarra_gamma_components(
            ped = ped,
            gamma = gamma,
            genotypedIds = genotypedIds,
            Z = Z,
            method = likelihood
        )$logLik
    }

    out <- list(
        gamma = gamma,
        logLik = final_logLik,
        genotyped = genotypedIds,
        id_source = legarra$id_source,
        groups = n_groups,
        method = "legarra",
        likelihood = likelihood,
        optimizer = "pseudo-em",
        convergence = if (converged) 0L else 1L,
        counts = c(iterations = iter),
        criterion = criterion
    )

    if (keep_history && iter > 0L) {
        out$history <- do.call(rbind, lapply(history[seq_len(iter)], as.data.frame))
    }

    out
}


#' Estimate a Metafounder Gamma Matrix by the Garcia-Baccino Method
#'
#' Estimate \code{gamma} directly from grouped genotype frequencies following
#' the Garcia-Baccino allele-frequency formulation.
#'
#' @param genotypicMatrix Genotype matrix with individuals in rows and markers
#'   in columns. The standard coding is \code{0/1/2}. Centered coding
#'   \code{-1/0/1} is also accepted when it can be identified unambiguously.
#'   If a matrix contains only \code{0/1} values, it is interpreted as
#'   standard \code{0/1/2} input.
#' @param metMatch Named vector assigning each row of \code{genotypicMatrix} to
#'   a metafounder. Its names must match
#'   \code{rownames(genotypicMatrix)}.
#' @return A list containing the estimated \code{gamma} matrix, raw allele
#'   frequencies \code{alleleFrequencies}, centered frequencies
#'   \code{pCentered = alleleFrequencies - 0.5}, resolved metafounder labels
#'   \code{metafounders}, and the estimator label \code{estimator = "naive"}.
#'   For backward compatibility, the aliases \code{frequencies} and
#'   \code{alleFrequency} are also returned.
#' @details \code{getGammaGarciaBaccino()} computes \eqn{\Gamma = 8/m FF^\top}
#'   when \eqn{F} denotes metafounder allele frequencies centered at 0.5. In
#'   this interface, those frequencies are always estimated internally from
#'   \code{genotypicMatrix} and the metafounder assignment supplied in
#'   \code{metMatch},
#'   and then centered internally. The current implementation is the direct
#'   group-wise sample-frequency estimator from genotypes.
#' @export
getGammaGarciaBaccino <- function(genotypicMatrix, metMatch) {
    if (is.null(genotypicMatrix)) {
        stop("genotypicMatrix must be provided for method = 'garcia-baccino'.")
    }

    gb <- garcia_baccino_naive(genotypicMatrix, metMatch = metMatch)

    list(
        gamma = gb$gamma,
        alleleFrequencies = gb$alleleFrequencies,
        pCentered = gb$pCentered,
        frequencies = gb$alleleFrequencies,
        alleFrequency = gb$pCentered,
        metafounders = gb$metafounders,
        estimator = gb$estimator,
        method = "garcia-baccino"
    )
}


#' Estimate a Metafounder Gamma Matrix
#'
#' Compatibility wrapper around \code{getGammaLegarra()} and
#' \code{getGammaGarciaBaccino()}.
#'
#' @param ped A \code{\link{pedigree}} object. Required for
#'   \code{method = "legarra"}.
#' @param genotypicMatrix Genotype matrix with individuals in rows and markers
#'   in columns. The standard coding is \code{0/1/2}. Centered coding
#'   \code{-1/0/1} is also accepted when it can be identified unambiguously.
#'   If a matrix contains only \code{0/1} values, it is interpreted as
#'   standard \code{0/1/2} input.
#' @param groups Method-dependent metafounder definition. For
#'   \code{method = "legarra"}, \code{groups} is a single positive integer
#'   giving the number of metafounders. For
#'   \code{method = "garcia-baccino"}, \code{groups} is retained only as a
#'   legacy alias of \code{metMatch}; new code should use \code{metMatch}.
#' @param metMatch Named vector assigning each row of
#'   \code{genotypicMatrix} to a metafounder for
#'   \code{method = "garcia-baccino"}. Its names must match
#'   \code{rownames(genotypicMatrix)}.
#' @param genotypedIds Character or integer IDs for the genotyped individuals.
#'   Used only for \code{method = "legarra"}.
#' @param method One of \code{"legarra"} or \code{"garcia-baccino"}.
#' @param likelihood One of \code{"byMarker"} or \code{"byIndividual"} for the
#'   exact Legarra log-likelihood. Used only for \code{method = "legarra"}.
#' @param idSource Controls where the IDs for the rows of
#'   \code{genotypicMatrix} come from in \code{method = "legarra"}.
#' @param control Optional named list used only for \code{method = "legarra"}
#'   and \code{groups > 1}.
#' @return A list containing the estimated \code{gamma} matrix and
#'   method-specific metadata.
#' @details New code should prefer \code{getGammaLegarra()} or
#'   \code{getGammaGarciaBaccino()} directly. \code{getGamma()} is kept as a
#'   single entry point for backward compatibility and dispatches to the
#'   appropriate method-specific function.
#' @export
getGamma <- function(ped = NULL, genotypicMatrix = NULL, groups = NULL, genotypedIds = NULL,
                     method = c("legarra", "garcia-baccino"),
                     likelihood = c("byMarker", "byIndividual"),
                     idSource = c("auto", "rownames", "genotyped"),
                     control = NULL,
                     metMatch = NULL) {
    method <- match.arg(method)

    if (method == "garcia-baccino") {
        if (!is.null(groups) && !is.null(metMatch)) {
            stop("For method = 'garcia-baccino', supply metMatch and leave groups = NULL.")
        }

        return(getGammaGarciaBaccino(
            genotypicMatrix = genotypicMatrix,
            metMatch = if (is.null(metMatch)) groups else metMatch
        ))
    }

    getGammaLegarra(
        ped = ped,
        genotypicMatrix = genotypicMatrix,
        groups = groups,
        genotypedIds = genotypedIds,
        likelihood = likelihood,
        idSource = idSource,
        control = control
    )
}
