% Manual documentation for gamma estimation entry points
\name{getGamma}
\alias{getGamma}
\alias{getGammaLegarra}
\alias{getGammaGarciaBaccino}
\title{Estimate Metafounder Gamma Matrices}
\usage{
getGammaLegarra(
  ped,
  genotypicMatrix,
  groups,
  genotypedIds = NULL,
  likelihood = c("byMarker", "byIndividual"),
  idSource = c("auto", "rownames", "genotyped"),
  control = NULL
)

getGammaGarciaBaccino(genotypicMatrix, metMatch)

getGamma(
  ped = NULL,
  genotypicMatrix = NULL,
  groups = NULL,
  genotypedIds = NULL,
  method = c("legarra", "garcia-baccino"),
  likelihood = c("byMarker", "byIndividual"),
  idSource = c("auto", "rownames", "genotyped"),
  control = NULL,
  metMatch = NULL
)
}
\arguments{
\item{ped}{A \code{\link{pedigree}} object. Required for
\code{getGammaLegarra()} and for \code{getGamma(method = "legarra")}.}

\item{genotypicMatrix}{Genotype matrix with individuals in rows and markers
in columns. The standard coding is \code{0/1/2}. Centered coding
\code{-1/0/1} is also accepted when it can be identified unambiguously. If a
matrix contains only \code{0/1} values, it is interpreted as standard
\code{0/1/2} input.}

\item{groups}{Method-dependent metafounder definition. For
\code{getGammaLegarra()} and \code{getGamma(method = "legarra")},
\code{groups} is a single positive integer giving the number of
metafounders. For \code{getGamma(method = "garcia-baccino")},
\code{groups} is retained only as a legacy alias of \code{metMatch};
new code should use \code{metMatch}.}

\item{metMatch}{Named vector assigning each row of
\code{genotypicMatrix} to a metafounder for the Garcia-Baccino workflow.
Its names must match \code{rownames(genotypicMatrix)}.}

\item{genotypedIds}{Character or integer IDs for the genotyped individuals.
Used only in the Legarra workflow. This is only required when the rows of
\code{genotypicMatrix} are not already identified by pedigree labels, or when
\code{idSource = "genotyped"} is used.}

\item{method}{One of \code{"legarra"} or \code{"garcia-baccino"} for the
\code{getGamma()} compatibility wrapper.}

\item{likelihood}{One of \code{"byMarker"} or \code{"byIndividual"} for the
exact Legarra log-likelihood. For \code{groups = 1}, it is the criterion
maximized by the scalar optimization. For \code{groups > 1}, the pseudo-EM
update does not depend on this choice. The exact log-likelihood is only
computed if requested by the user, either through \code{logLikeGamma()} or
via \code{control$compute_logLik = TRUE}.}

\item{idSource}{Controls where the IDs for the rows of
\code{genotypicMatrix} come from in the Legarra workflow.
Use \code{"rownames"} to read them from \code{rownames(genotypicMatrix)},
use \code{"genotyped"} to read them from \code{genotypedIds} in row
order, and use \code{"auto"} to prefer \code{rownames(genotypicMatrix)}
when they match \code{ped@label} and otherwise fall back to
\code{genotypedIds}. If both are available and \code{idSource = "auto"},
the row names take precedence.}

\item{control}{Optional named list used only for the multi-metafounder
Legarra pseudo-EM workflow. Supported entries are \code{maxit},
\code{tol}, \code{keep_history}, and \code{compute_logLik}. Defaults are
\code{list(maxit = 200, tol = 1e-8, keep_history = FALSE, compute_logLik = FALSE)}.
For example,
\code{control = list(maxit = 500, tol = 1e-10, compute_logLik = TRUE)}.}
}
\value{
A list containing the estimated \code{gamma} matrix and method-specific
metadata. For the Legarra workflow, the result includes the resolved
\code{genotyped} IDs and iteration metadata. The exact returned
\code{logLik} is included only when \code{control$compute_logLik = TRUE} for
\code{groups > 1}, or always for the scalar \code{groups = 1} case. For the
Garcia-Baccino workflow, the result includes the raw allele frequencies
\code{alleleFrequencies}, the centered frequencies
\code{pCentered = alleleFrequencies - 0.5}, the resolved metafounder labels
\code{metafounders}, and the estimator label \code{estimator = "naive"}.
For backward compatibility, the aliases \code{frequencies} and
\code{alleFrequency} are also returned.}
\description{
\code{getGammaLegarra()} and \code{getGammaGarciaBaccino()} are the
method-specific public entry points for estimating metafounder gamma
matrices. \code{getGamma()} is a compatibility wrapper that dispatches to one
of those two functions according to \code{method}.
}
\details{
\code{getGammaLegarra()} supports one or more metafounders. Their number is
provided through \code{groups} as a single positive integer. When possible,
the genotyped IDs are taken directly from \code{rownames(genotypicMatrix)};
otherwise they must be supplied through \code{genotypedIds}. The exact
Legarra log-likelihood is evaluated with the implicit Masuda sequence in
block form over the full marker matrix, avoiding explicit construction of
\eqn{A_{22}^{-1}} and avoiding a per-marker R loop. Genotypes coded as
\code{0/1/2} are internally recoded to \code{Z = M - 1}; centered input
coded as \code{-1/0/1} is used directly when detected unambiguously. When
\code{groups = 1}, \code{gamma} is a \eqn{1 \times 1} matrix and the
estimator uses the exact single-metafounder cubic equation from Legarra over
the scalar \eqn{\gamma_{11}} on a fixed internal interval. When
\code{groups > 1}, the estimate is obtained by the pseudo-EM
update described by Legarra (2024), while the exact log-likelihood remains
available through \code{logLikeGamma()} for verification. The
\code{likelihood} argument chooses only how that exact log-likelihood is
evaluated: by marker or by individual. The pseudo-EM iterations themselves
use only their own update and convergence criterion; exact log-likelihood
checks are optional because they can be costly for large matrices.

\code{getGammaGarciaBaccino()} computes \eqn{\Gamma = 8/m FF^\top} when
\eqn{F} denotes metafounder allele frequencies centered at 0.5. In this
interface, those frequencies are always estimated internally from
\code{genotypicMatrix} and the metafounder assignment supplied in
\code{metMatch}, and then centered internally. The current implementation is
the direct group-wise sample-frequency estimator from genotypes, returned as
\code{estimator = "naive"}.

\code{getGamma()} is kept as a single entry point for backward compatibility,
but new code should prefer the method-specific functions directly.
}
