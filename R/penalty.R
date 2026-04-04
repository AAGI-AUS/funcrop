# penalty.R -- Penalty matrix construction and mixed model reparameterisation
#
# Constructs difference penalty matrices for P-splines and provides the
# eigendecomposition/Cholesky reparameterisation that maps penalised splines
# to mixed model random effects (null space -> X, range space -> Z).
#
# Author: Maksym Gaidashenko
# Licence: GPL (>= 3)


# ---- make_penalty ----

#' Construct a Difference Penalty Matrix
#'
#' Builds the penalty matrix \eqn{P = D_d^T D_d} where \eqn{D_d} is the
#' d-th order difference matrix of dimension `(n_basis - d) x n_basis`. This
#' is the standard discrete penalty used in P-spline smoothing (Eilers & Marx,
#' 1996).
#'
#' @param n_basis Integer; number of basis functions (columns of the B-spline
#'   basis matrix). Must be >= 2.
#' @param order Integer; order of the difference penalty. Default is 2
#'   (second-order, penalising curvature). `order = 1` corresponds to a random
#'   walk penalty. Must satisfy `1 <= order < n_basis`.
#'
#' @return A sparse penalty matrix of class `dgCMatrix` (from the Matrix
#'   package), dimension `n_basis x n_basis`. The matrix is symmetric and
#'   positive semi-definite with rank `n_basis - order`.
#'
#' @details
#' For `order = 1`, the difference matrix \eqn{D_1} has entries
#' \eqn{(D_1)_{ij} = 1} if \eqn{j = i+1} and \eqn{-1} if \eqn{j = i}.
#' Higher-order difference matrices are obtained by recursive application:
#' \eqn{D_d = D_1^{(n-d+1)} D_{d-1}}, where the superscript denotes the
#' dimension of the first-order matrix.
#'
#' @references
#' Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines
#' and penalties. *Statistical Science*, 11(2), 89--121.
#'
#' @seealso [bspline_basis()] for basis construction,
#'   [make_Zspline()] for mixed model reparameterisation.
#'
#' @examples
#' # Second-order penalty for 12 basis functions
#' P <- make_penalty(12, order = 2)
#' dim(P)
#' Matrix::rankMatrix(P)
#'
#' # First-order penalty (random walk)
#' P1 <- make_penalty(12, order = 1)
#'
#' @export
make_penalty <- function(n_basis, order = 2L) {

  # ---- Input validation ----

  n_basis <- as.integer(n_basis)
  if (length(n_basis) != 1L || is.na(n_basis) || n_basis < 2L) {
    stop("`n_basis` must be a single integer >= 2.", call. = FALSE)
  }

  order <- as.integer(order)
  if (length(order) != 1L || is.na(order) || order < 1L) {
    stop("`order` must be a single positive integer.", call. = FALSE)
  }
  if (order >= n_basis) {
    stop("`order` (", order, ") must be < `n_basis` (", n_basis, ").",
         call. = FALSE)
  }

  # ---- Build difference matrix recursively ----
  # Start with identity of size n_basis, apply first-order diff `order` times
  D <- Matrix::Diagonal(n_basis)

  for (i in seq_len(order)) {
    m <- nrow(D)
    # First-order difference matrix of size (m - 1) x m
    D1 <- Matrix::bandSparse(
      n = m - 1L,
      m = m,
      k = c(0L, 1L),
      diagonals = list(rep(-1, m - 1L), rep(1, m - 1L))
    )
    D <- D1 %*% D
  }

  # Penalty: D'D -- symmetric, positive semi-definite
  P <- Matrix::crossprod(D)

  # Ensure dgCMatrix storage for consistency
  P <- methods::as(P, "CsparseMatrix")

  P
}


# ---- make_Zspline ----

#' Reparameterise B-Spline Basis for Mixed Model Representation
#'
#' Decomposes a B-spline basis into fixed-effect (X) and random-effect (Z)
#' components, enabling the equivalence between penalised splines and linear
#' mixed models (Wand, 2003; Wood, 2017). This is essential for fitting
#' P-splines via REML in mixed model software (e.g., ASReml-R, lme4).
#'
#' @param basis An object of class `"fda_basis"` (as returned by
#'   [bspline_basis()]).
#' @param constraint Character; reparameterisation method. One of:
#' \describe{
#'   \item{`"decompose"`}{(Default) Spectral decomposition of the penalty
#'     matrix. Eigenvectors corresponding to zero eigenvalues form the fixed
#'     effects (X, the null space), whilst those with positive eigenvalues form
#'     the random effects (Z, the range space), scaled by the inverse square
#'     root of their eigenvalues.}
#'   \item{`"absorb"`}{Absorb the penalty into Z via Cholesky decomposition of
#'     the Moore-Penrose pseudoinverse. The resulting Z has an identity penalty
#'     structure, which is convenient for standard mixed model software.}
#' }
#'
#' @return A list with components:
#' \describe{
#'   \item{X}{Fixed-effect design matrix (n x d), where d = `penalty_order`
#'     is the dimension of the null space of the penalty. Typically contains
#'     polynomial terms up to degree `penalty_order - 1`.}
#'   \item{Z}{Random-effect design matrix (n x r), where r = `n_basis - d`.
#'     In the `"decompose"` method, columns are orthogonal; in the `"absorb"`
#'     method, the penalty is absorbed so that the implied penalty on the
#'     random effect coefficients is the identity.}
#'   \item{penalty}{Sparse penalty matrix for the Z component. Identity matrix
#'     (dimension r x r) for the `"absorb"` method; identity for `"decompose"`
#'     when eigenvalues are appropriately absorbed into Z.}
#'   \item{rank}{Integer; rank of the original penalty matrix
#'     (`n_basis - penalty_order`).}
#'   \item{null_dim}{Integer; dimension of the null space (= `penalty_order`).}
#'   \item{constraint}{Character; the method used.}
#' }
#'
#' @details
#' The key identity exploited is: a penalised spline fit minimising
#' \deqn{\|y - B\beta\|^2 + \lambda \beta^T P \beta}
#' is equivalent to the BLUP from a mixed model
#' \deqn{y = X\alpha + Z u + \varepsilon}
#' where \eqn{u \sim N(0, \sigma_u^2 I)} and the smoothing parameter
#' \eqn{\lambda = \sigma^2 / \sigma_u^2} is estimated via REML.
#'
#' For `"decompose"`: the spectral decomposition \eqn{P = U \Lambda U^T}
#' separates the null space (zero eigenvalues, giving X) from the range space
#' (positive eigenvalues, giving Z scaled by \eqn{\Lambda^{-1/2}}).
#'
#' For `"absorb"`: uses the pseudoinverse \eqn{P^+} and its Cholesky factor
#' \eqn{L} such that \eqn{P^+ = LL^T}. Then \eqn{Z = B_{\text{range}} L}
#' absorbs the penalty, yielding identity covariance.
#'
#' @references
#' Wand, M.P. (2003). Smoothing and mixed models. *Computational Statistics*,
#' 18(2), 223--249.
#'
#' Wood, S.N. (2017). *Generalized Additive Models: An Introduction with R*
#' (2nd ed.). Chapman & Hall/CRC.
#'
#' @seealso [bspline_basis()], [make_penalty()].
#'
#' @examples
#' x <- seq(0, 1, length.out = 100)
#' basis <- bspline_basis(x, n_knots = 10, degree = 3)
#'
#' # Decompose method (default)
#' mm <- make_Zspline(basis, constraint = "decompose")
#' dim(mm$X)
#' dim(mm$Z)
#'
#' # Absorb method -- penalty becomes identity
#' mm2 <- make_Zspline(basis, constraint = "absorb")
#' dim(mm2$Z)
#'
#' @export
make_Zspline <- function(basis,
                         constraint = c("decompose", "absorb")) {

  # ---- Input validation ----

  if (!inherits(basis, "fda_basis")) {
    stop("`basis` must be an object of class 'fda_basis'.", call. = FALSE)
  }

  constraint <- match.arg(constraint)

  B <- basis$B
  P <- basis$P
  n_basis <- basis$n_basis
  d <- basis$penalty_order  # null space dimension

  if (d >= n_basis) {
    stop("Penalty order (", d, ") must be < n_basis (", n_basis,
         ") for reparameterisation.", call. = FALSE)
  }

  rank_P <- n_basis - d
  n <- nrow(B)

  # ---- Spectral decomposition of penalty ----
  # Convert to dense for eigen (penalty matrices are moderate-sized)
  P_dense <- as.matrix(P)

  eig <- eigen(P_dense, symmetric = TRUE)
  eigenvalues <- eig$values
  eigenvectors <- eig$vectors

  # Numerical tolerance for zero eigenvalues
  tol <- max(eigenvalues) * .Machine$double.eps * n_basis * 10

  # Indices: positive eigenvalues (range space) first in eigen output,

  # then zero eigenvalues (null space)
  idx_range <- which(eigenvalues > tol)
  idx_null <- which(eigenvalues <= tol)

  # Sanity check
  if (length(idx_null) != d) {
    warning(
      "Expected ", d, " zero eigenvalues but found ", length(idx_null),
      ". Proceeding with detected null space dimension.",
      call. = FALSE
    )
    d <- length(idx_null)
    rank_P <- n_basis - d
  }

  # ---- Fixed-effect component (null space) ----
  U_null <- eigenvectors[, idx_null, drop = FALSE]
  X <- B %*% U_null

  # ---- Random-effect component (range space) ----
  U_range <- eigenvectors[, idx_range, drop = FALSE]
  lambda_range <- eigenvalues[idx_range]

  if (constraint == "decompose") {
    # Scale range eigenvectors by inverse sqrt of eigenvalues
    # so that the implied penalty on u is the identity
    scaling <- 1 / sqrt(lambda_range)
    Z <- B %*% (U_range %*% diag(scaling, nrow = length(scaling)))

    penalty_out <- Matrix::Diagonal(rank_P)

  } else {
    # "absorb" method
    # B_range = B %*% U_range has penalty diag(lambda_range)
    # Absorb via: Z = B_range %*% L where L = chol(inv(diag(lambda)))
    # i.e., L = diag(1/sqrt(lambda))
    # This gives Z such that u ~ N(0, sigma_u^2 I) directly
    scaling <- 1 / sqrt(lambda_range)
    Z <- B %*% (U_range %*% diag(scaling, nrow = length(scaling)))

    penalty_out <- Matrix::Diagonal(rank_P)
  }

  # Ensure standard matrix types
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  list(
    X          = X,
    Z          = Z,
    penalty    = penalty_out,
    rank       = rank_P,
    null_dim   = d,
    constraint = constraint
  )
}
