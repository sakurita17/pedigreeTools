#' Genomic Relationship Matrix Using VanRaden Methods 1 and 2
#'
#' Compute the genomic relationship matrix \eqn{G} from a genotype matrix coded
#' as 0, 1, or 2 copies of a reference allele using VanRaden (2008).
#'
#' @param genotypicMatrix Numeric matrix or \code{big.matrix} with individuals
#'   in rows and markers in columns, coded as 0, 1, or 2.
#' @param method Character string indicating whether to use VanRaden method 1
#'   (\code{"VR1"}) or method 2 (\code{"VR2"}).
#' @param alleleFreq Optional numeric vector with one allele frequency per
#'   marker. If supplied, \code{freqMode} is ignored.
#' @param freqMode Character string used only when \code{alleleFreq = NULL}.
#'   \code{"0.5"} centers each marker at 0.5, whereas \code{"observed"} uses
#'   the observed allele frequencies \code{colMeans(genotypicMatrix) / 2}.
#' @param epsilon Positive ridge value added to the diagonal of \eqn{G}.
#'
#' @details
#' The function first chooses allele frequencies \eqn{p} and centers the
#' genotype matrix as \eqn{Z = M - 2p}. With \code{method = "VR1"}, all markers
#' share the same global denominator \eqn{2 \sum p_i (1 - p_i)}. With
#' \code{method = "VR2"}, each marker is standardized by its expected variance
#' through weights \eqn{1 / [2 p_i (1 - p_i)]} before averaging across markers.
#' If row names are present, they are expected to be character and are copied to
#' the output matrix.
#'
#' @return A dense numeric matrix with the genomic relationship coefficients.
#' @references VanRaden, P. M. (2008). Efficient Methods to Compute Genomic
#'   Predictions. \emph{Journal of Dairy Science}, 91(11), 4414-4423.
#'
#' @export
getG <- function(genotypicMatrix,
                 method = c("VR1", "VR2"),
                 alleleFreq = NULL,
                 freqMode = c("0.5", "observed"),
                 epsilon = .Machine$double.eps) {
    
    method <- match.arg(method)
    freqMode <- match.arg(freqMode)

    genotypes <- prepare_genotype_matrix(genotypicMatrix, arg = "genotypicMatrix")
    genotypicMatrix <- genotypes$source

    if (!genotypes$is_bigmat && !is.numeric(genotypicMatrix)) {
        stop("genotypicMatrix must be a numeric matrix or big.matrix.")
    }
    if (!is.numeric(epsilon) || length(epsilon) != 1L || is.na(epsilon) || epsilon <= 0) {
        stop("epsilon must be a single positive numeric value.")
    }

    matrixIds <- row_ids(genotypicMatrix, arg = "genotypicMatrix", required = FALSE)
    markerCount <- genotypes$markerCount

    if (!genotypes$is_bigmat) {
        if (anyNA(genotypicMatrix)) {
            stop("genotypicMatrix must not contain NA values.")
        }
        if (any(genotypicMatrix < 0 | genotypicMatrix > 2)) {
            stop("genotypicMatrix must contain only values between 0 and 2 (inclusive).")
        }

        # If the user supplies allele frequencies, they override freqMode.
        if (is.null(alleleFreq)) {
            if (freqMode == "observed") {
                p <- colMeans(genotypicMatrix) / 2
            } else {
                p <- rep(0.5, markerCount)
            }
        } else {
            if (!is.numeric(alleleFreq)) {
                stop("alleleFreq must be a numeric vector.")
            }
            if (length(alleleFreq) != markerCount) {
                stop("alleleFreq must have length ncol(genotypicMatrix).")
            }
            p <- alleleFreq
        }

        if (anyNA(p) || any(p < 0 | p > 1)) {
            stop("alleleFreq values must be between 0 and 1.")
        }

        # Center each marker using twice its allele frequency.
        Z <- sweep(genotypicMatrix, 2L, 2 * p, FUN = "-")

        if (method == "VR1") {
            denominator <- 2 * sum(p * (1 - p))
            if (!is.finite(denominator) || denominator <= 0) {
                stop("VR1 requires a positive scaling denominator.")
            }
            G <- tcrossprod(Z) / denominator
        } else {
            if (any(p <= 0 | p >= 1)) {
                stop("VR2 requires allele frequencies strictly between 0 and 1.")
            }

            markerWeights <- 1 / (2 * p * (1 - p))
            weightedZ <- sweep(Z, 2L, sqrt(markerWeights), FUN = "*")
            G <- tcrossprod(weightedZ) / markerCount
        }

        diag(G) <- diag(G) + epsilon
        dimnames(G) <- list(matrixIds, matrixIds)

        return(G)
    }

    genotype_block_apply(genotypes, function(block, cols) {
        if (anyNA(block)) {
            stop("genotypicMatrix must not contain NA values.")
        }
        if (any(block < 0 | block > 2)) {
            stop("genotypicMatrix must contain only values between 0 and 2 (inclusive).")
        }
    })

    if (is.null(alleleFreq)) {
        if (freqMode == "observed") {
            p <- numeric(markerCount)
            genotype_block_apply(genotypes, function(block, cols) {
                p[cols] <<- colMeans(block) / 2
            })
        } else {
            p <- rep(0.5, markerCount)
        }
    } else {
        if (!is.numeric(alleleFreq)) {
            stop("alleleFreq must be a numeric vector.")
        }
        if (length(alleleFreq) != markerCount) {
            stop("alleleFreq must have length ncol(genotypicMatrix).")
        }
        p <- alleleFreq
    }

    if (anyNA(p) || any(p < 0 | p > 1)) {
        stop("alleleFreq values must be between 0 and 1.")
    }

    G <- matrix(0, nrow = nrow(genotypicMatrix), ncol = nrow(genotypicMatrix))

    if (method == "VR1") {
        denominator <- 2 * sum(p * (1 - p))
        if (!is.finite(denominator) || denominator <= 0) {
            stop("VR1 requires a positive scaling denominator.")
        }

        genotype_block_apply(genotypes, function(block, cols) {
            Z_block <- sweep(block, 2L, 2 * p[cols], FUN = "-")
            G <<- G + tcrossprod(Z_block)
        })
        G <- G / denominator
    } else {
        if (any(p <= 0 | p >= 1)) {
            stop("VR2 requires allele frequencies strictly between 0 and 1.")
        }

        markerWeights <- 1 / (2 * p * (1 - p))
        genotype_block_apply(genotypes, function(block, cols) {
            Z_block <- sweep(block, 2L, 2 * p[cols], FUN = "-")
            weightedZ <- sweep(Z_block, 2L, sqrt(markerWeights[cols]), FUN = "*")
            G <<- G + tcrossprod(weightedZ)
        })
        G <- G / markerCount
    }

    diag(G) <- diag(G) + epsilon
    dimnames(G) <- list(matrixIds, matrixIds)

    G
}
