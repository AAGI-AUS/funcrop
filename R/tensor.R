# tensor.R — Tensor product B-spline basis for 2D functional data
#
# Constructs tensor product bases from two marginal B-spline bases. The
# tensor product basis is formed by the row-wise Kronecker product of the
# marginal basis matrices, with an additive penalty structure.
#
# This is used for 2D functional traits: time x depth, time x wavelength, etc.
#
# Author: Maksym Gaidashenko
# Licence: GPL (>= 3)


# ---- tensor_bspline_basis ----

#' Construct a Tensor Product B-Spline Basis
#'
#' Builds a two-dimensional tensor product basis from two marginal B-spline
#' bases, suitable for modelling functional data observed over two indices
#' (e.g., time x depth, time x wavelength). The penalty is the standard
#' additive Kronecker sum: \eqn{P = P_1 \otimes I_2 + I_1 \otimes P_2}.
#'
#' @param x1 Numeric vector of evaluation points for the first dimension.
#' @param x2 Numeric vector of evaluation points for the second dimension.
#'   Must have the same length as `x1` (i.e., each observation has
#'   coordinates `(x1[i], x2[i])`).
#' @param n_knots1,n_knots2 Integer; number of internal knots for each
#'   marginal basis. Defaults are 10.
#' @param degree1,degree2 Integer; polynomial degree for each marginal basis.
#'   Defaults are 3 (cubic).
#' @param knot_type1,knot_type2 Character; knot placement for each marginal
#'   basis. One of `"equally_spaced"` (default), `"quantile"`, or `"custom"`.
#' @param penalty_order1,penalty_order2 Integer; order of the difference
#'   penalty for each marginal basis. Defaults are 2.
#'
#' @return An S3 object of class `c("fda_tensor_basis", "fda_basis")`
#'   containing:
#' \describe{
#'   \item{B}{Tensor product basis matrix (n x p1*p2), constructed as the
#'     row-wise Kronecker product of the two marginal basis matrices.}
#'   \item{basis1, basis2}{The marginal `fda_basis` objects for each
#'     dimension.}
#'   \item{P}{Combined additive penalty matrix (p1*p2 x p1*p2), sparse.
#'     Computed as \eqn{P_1 \otimes I_2 + I_1 \otimes P_2}.}
#'   \item{P1}{Marginal penalty in tensor product form:
#'     \eqn{P_1 \otimes I_2} (sparse).}
#'   \item{P2}{Marginal penalty in tensor product form:
#'     \eqn{I_1 \otimes P_2} (sparse).}
#'   \item{n_basis}{Integer; total number of tensor product basis functions
#'     (p1 * p2).}
#' }
#'
#' @details
#' The tensor product basis matrix at observation i is:
#' \deqn{B_{\text{tensor}}[i,] = B_1[i,] \otimes B_2[i,]}
#' where \eqn{\otimes} denotes the Kronecker product of two row vectors.
#'
#' The additive penalty:
#' \deqn{P = P_1 \otimes I_2 + I_1 \otimes P_2}
#' penalises roughness in each dimension separately, which is the standard
#' approach for tensor product P-splines (Marx & Eilers, 2005; Wood, 2006).
#'
#' Sparse matrix representations (Matrix package) are used throughout as
#' tensor product bases can be large (p1 * p2 basis functions).
#'
#' @references
#' Marx, B.D. and Eilers, P.H.C. (2005). Multidimensional penalized signal
#' regression. *Technometrics*, 47(1), 13--22.
#'
#' Wood, S.N. (2006). Low-rank scale-invariant tensor product smooths for
#' generalized additive mixed models. *Biometrics*, 62(4), 1025--1036.
#'
#' @seealso [bspline_basis()] for marginal basis construction,
#'   [make_penalty()] for penalty matrices,
#'   [make_Zspline()] for mixed model reparameterisation.
#'
#' @examples
#' # 2D basis: time x depth with 100 observations
#' set.seed(123)
#' n <- 100
#' time_pts <- runif(n, 0, 10)
#' depth_pts <- runif(n, 0, 1)
#' tb <- tensor_bspline_basis(time_pts, depth_pts,
#'                            n_knots1 = 6, n_knots2 = 5)
#' print(tb)
#' dim(tb$B)
#'
#' \donttest{
#' # Larger tensor product basis (may be slow)
#' tb_large <- tensor_bspline_basis(time_pts, depth_pts,
#'                                  n_knots1 = 15, n_knots2 = 10)
#' dim(tb_large$B)
#' }
#'
#' @export
tensor_bspline_basis <- function(x1,
                                 x2,
                                 n_knots1 = 10L,
                                 n_knots2 = 10L,
                                 degree1 = 3L,
                                 degree2 = 3L,
                                 knot_type1 = "equally_spaced",
                                 knot_type2 = "equally_spaced",
                                 penalty_order1 = 2L,
                                 penalty_order2 = 2L) {

  # ---- Input validation ----

  if (!is.numeric(x1)) {
    stop("`x1` must be a numeric vector.", call. = FALSE)
  }
  if (!is.numeric(x2)) {
    stop("`x2` must be a numeric vector.", call. = FALSE)
  }
  if (length(x1) != length(x2)) {
    stop("`x1` and `x2` must have the same length (paired observations).",
         call. = FALSE)
  }
  if (length(x1) < 2L) {
    stop("`x1` and `x2` must each have length >= 2.", call. = FALSE)
  }
  if (any(!is.finite(x1))) {
    stop("`x1` must not contain NA, NaN, or Inf values.", call. = FALSE)
  }
  if (any(!is.finite(x2))) {
    stop("`x2` must not contain NA, NaN, or Inf values.", call. = FALSE)
  }

  # Validate knot_type arguments
  valid_types <- c("equally_spaced", "quantile", "custom")
  if (!is.character(knot_type1) || length(knot_type1) != 1L ||
      !knot_type1 %in% valid_types) {
    stop("`knot_type1` must be one of: ",
         paste0("'", valid_types, "'", collapse = ", "), ".", call. = FALSE)
  }
  if (!is.character(knot_type2) || length(knot_type2) != 1L ||
      !knot_type2 %in% valid_types) {
    stop("`knot_type2` must be one of: ",
         paste0("'", valid_types, "'", collapse = ", "), ".", call. = FALSE)
  }

  # ---- Construct marginal bases ----

  basis1 <- bspline_basis(
    x         = x1,
    n_knots   = n_knots1,
    degree    = degree1,
    knot_type = knot_type1,
    penalty_order = penalty_order1
  )

  basis2 <- bspline_basis(
    x         = x2,
    n_knots   = n_knots2,
    degree    = degree2,
    knot_type = knot_type2,
    penalty_order = penalty_order2
  )

  p1 <- basis1$n_basis
  p2 <- basis2$n_basis
  n_total <- p1 * p2

  # ---- Tensor product basis matrix (row-wise Kronecker) ----

  B_tensor <- .row_kronecker(basis1$B, basis2$B)

  # ---- Penalty matrices ----
  # P = P1 (x) I2 + I1 (x) P2   (Kronecker sum)

  I1 <- Matrix::Diagonal(p1)
  I2 <- Matrix::Diagonal(p2)
  P1_marg <- basis1$P
  P2_marg <- basis2$P

  # Marginal penalties in tensor product form
  P1_tensor <- Matrix::kronecker(P1_marg, I2)
  P2_tensor <- Matrix::kronecker(I1, P2_marg)

  # Combined additive penalty
  P_combined <- P1_tensor + P2_tensor

  # Ensure consistent sparse storage
  P_combined <- methods::as(P_combined, "CsparseMatrix")
  P1_tensor <- methods::as(P1_tensor, "CsparseMatrix")
  P2_tensor <- methods::as(P2_tensor, "CsparseMatrix")

  # ---- Assemble and return ----

  out <- list(
    B       = B_tensor,
    basis1  = basis1,
    basis2  = basis2,
    P       = P_combined,
    P1      = P1_tensor,
    P2      = P2_tensor,
    n_basis = n_total
  )
  class(out) <- c("fda_tensor_basis", "fda_basis")
  out
}


# ---- print.fda_tensor_basis ----

#' Print Method for fda_tensor_basis Objects
#'
#' Displays a compact summary of a tensor product B-spline basis.
#'
#' @param x An object of class `"fda_tensor_basis"`.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @examples
#' set.seed(42)
#' n <- 50
#' tb <- tensor_bspline_basis(runif(n), runif(n), n_knots1 = 5, n_knots2 = 4)
#' print(tb)
#'
#' @export
print.fda_tensor_basis <- function(x, ...) {

  b1 <- x$basis1
  b2 <- x$basis2

  cat("Tensor product B-spline basis (fda_tensor_basis)\n")
  cat("  Dimension 1:\n")
  cat("    Degree:         ", b1$degree, "\n")
  cat("    Internal knots: ", length(b1$knots), " (", b1$knot_type, ")\n", sep = "")
  cat("    Basis functions:", b1$n_basis, "\n")
  cat("    Boundary:       [", format(b1$boundary[1], digits = 4), ", ",
      format(b1$boundary[2], digits = 4), "]\n", sep = "")
  cat("    Penalty order:  ", b1$penalty_order, "\n")
  cat("  Dimension 2:\n")
  cat("    Degree:         ", b2$degree, "\n")
  cat("    Internal knots: ", length(b2$knots), " (", b2$knot_type, ")\n", sep = "")
  cat("    Basis functions:", b2$n_basis, "\n")
  cat("    Boundary:       [", format(b2$boundary[1], digits = 4), ", ",
      format(b2$boundary[2], digits = 4), "]\n", sep = "")
  cat("    Penalty order:  ", b2$penalty_order, "\n")
  cat("  Total basis functions:", x$n_basis,
      " (", b1$n_basis, " x ", b2$n_basis, ")\n", sep = "")
  cat("  Observations:        ", nrow(x$B), "\n")
  cat("  Basis matrix dim:    ", nrow(x$B), " x ", ncol(x$B), "\n", sep = "")

  invisible(x)
}
