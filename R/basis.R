# basis.R — Core B-spline basis constructor for funcrop
#
# Constructs B-spline basis matrices with integrated difference penalty
# matrices for use in P-spline and mixed model representations.
#
# Author: Maksym Gaidashenko
# Licence: GPL (>= 3)


# ---- bspline_basis ----

#' Construct a B-Spline Basis with Difference Penalty
#'
#' Creates a B-spline basis matrix evaluated at a set of points, together with
#' a difference penalty matrix suitable for P-spline smoothing. This is the
#' core building block for all functional data representations in `funcrop`.
#'
#' @param x Numeric vector of evaluation points (e.g., time, wavelength, depth).
#'   Must contain at least 2 unique values and no `NA` or non-finite values.
#' @param n_knots Integer; number of internal knots. Default is 10. Must be
#'   >= 1. The total number of basis functions is `n_knots + degree + 1`.
#' @param degree Integer; polynomial degree of the B-spline. Default is 3
#'   (cubic). Must be >= 0.
#' @param knot_type Character; method for placing internal knots. One of
#'   `"equally_spaced"` (default), `"quantile"`, or `"custom"`. If `"custom"`,
#'   the `knots` argument must be supplied.
#' @param knots Numeric vector of internal knots. Required when
#'   `knot_type = "custom"`, ignored otherwise. Must lie strictly within the
#'   boundary knots.
#' @param boundary Numeric vector of length 2 giving the boundary knots
#'   (lower, upper). Defaults to `range(x)` extended by 1% on each side.
#' @param penalty_order Integer; order of the difference penalty matrix.
#'   Default is 2 (second-order differences, penalising curvature). Must be
#'   >= 1 and < `n_basis`.
#'
#' @return An S3 object of class `"fda_basis"` (a list) containing:
#' \describe{
#'   \item{B}{Basis matrix (n x p), where n = `length(x)` and
#'     p = `n_basis` = `n_knots + degree + 1`.}
#'   \item{knots}{Numeric vector of internal knots used.}
#'   \item{boundary}{Numeric vector of length 2; boundary knots.}
#'   \item{degree}{Integer; polynomial degree.}
#'   \item{penalty_order}{Integer; order of the difference penalty.}
#'   \item{P}{Penalty matrix (p x p); sparse `Matrix` of class `dgCMatrix`.
#'     Constructed as D'D where D is the difference matrix of the specified
#'     order.}
#'   \item{x}{Numeric vector; the original evaluation points.}
#'   \item{n_basis}{Integer; total number of basis functions.}
#'   \item{knot_type}{Character; knot placement method used.}
#' }
#'
#' @details
#' The basis matrix is computed using [splines::splineDesign()], which handles
#' boundary evaluation correctly. Knots are augmented with repeated boundary
#' knots as required by the B-spline recursion.
#'
#' The penalty matrix \eqn{P = D_d^T D_d} encodes a discrete approximation to
#' the d-th derivative penalty, where \eqn{D_d} is the d-th order difference
#' matrix operating on the B-spline coefficients (Eilers & Marx, 1996).
#'
#' @references
#' Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines
#' and penalties. *Statistical Science*, 11(2), 89--121.
#'
#' @seealso [make_penalty()] for standalone penalty construction,
#'   [make_Zspline()] for mixed model reparameterisation,
#'   [tensor_bspline_basis()] for 2D tensor product bases.
#'
#' @examples
#' # Basic cubic B-spline basis with 10 internal knots
#' x <- seq(0, 1, length.out = 100)
#' basis <- bspline_basis(x, n_knots = 10, degree = 3)
#' print(basis)
#'
#' # Quantile-based knot placement (better for skewed data)
#' set.seed(42)
#' x_skew <- sort(rbeta(200, 2, 5))
#' basis_q <- bspline_basis(x_skew, n_knots = 8, knot_type = "quantile")
#'
#' \donttest{
#' # Visualise basis functions
#' plot(basis)
#' }
#'
#' @export
bspline_basis <- function(x,
                          n_knots = 10L,
                          degree = 3L,
                          knot_type = c("equally_spaced", "quantile", "custom"),
                          knots = NULL,
                          boundary = NULL,
                          penalty_order = 2L) {


  # ---- Input validation ----

  if (!is.numeric(x)) {
    stop("`x` must be a numeric vector.", call. = FALSE)
  }
  if (length(x) < 2L) {
    stop("`x` must have length >= 2.", call. = FALSE)
  }

if (any(!is.finite(x))) {
    stop("`x` must not contain NA, NaN, or Inf values.", call. = FALSE)
  }
  if (length(unique(x)) < 2L) {
    stop("`x` must contain at least 2 unique values.", call. = FALSE)

  }

  n_knots <- as.integer(n_knots)
  if (length(n_knots) != 1L || is.na(n_knots) || n_knots < 1L) {
    stop("`n_knots` must be a single positive integer.", call. = FALSE)
  }

  degree <- as.integer(degree)
  if (length(degree) != 1L || is.na(degree) || degree < 0L) {
    stop("`degree` must be a single non-negative integer.", call. = FALSE)
  }

  knot_type <- match.arg(knot_type)

  penalty_order <- as.integer(penalty_order)
  if (length(penalty_order) != 1L || is.na(penalty_order) || penalty_order < 1L) {
    stop("`penalty_order` must be a single positive integer.", call. = FALSE)
  }

  # ---- Boundary knots ----

  if (is.null(boundary)) {
    x_range <- range(x)
    x_span <- diff(x_range)
    # Extend by 1% on each side for numerical stability at boundaries
    boundary <- x_range + c(-1, 1) * 0.01 * x_span
  } else {
    if (!is.numeric(boundary) || length(boundary) != 2L) {
      stop("`boundary` must be a numeric vector of length 2.", call. = FALSE)
    }
    if (any(!is.finite(boundary))) {
      stop("`boundary` must not contain NA, NaN, or Inf values.", call. = FALSE)
    }
    boundary <- sort(boundary)
    if (boundary[1] >= boundary[2]) {
      stop("`boundary` knots must define a non-degenerate interval.", call. = FALSE)
    }
    # Check x values lie within boundary (with small tolerance)
    tol <- 1e-10 * diff(boundary)
    if (any(x < boundary[1] - tol) || any(x > boundary[2] + tol)) {
      stop("All values in `x` must lie within (or very near) the boundary knots.",
           call. = FALSE)
    }
  }

  # ---- Internal knots ----

  if (knot_type == "custom") {
    if (is.null(knots)) {
      stop("`knots` must be supplied when knot_type = 'custom'.", call. = FALSE)
    }
    if (!is.numeric(knots)) {
      stop("`knots` must be a numeric vector.", call. = FALSE)
    }
    if (any(!is.finite(knots))) {
      stop("`knots` must not contain NA, NaN, or Inf values.", call. = FALSE)
    }
    knots <- sort(knots)
    if (any(knots <= boundary[1]) || any(knots >= boundary[2])) {
      stop("All internal `knots` must lie strictly within the boundary knots.",
           call. = FALSE)
    }
    n_knots <- length(knots)
  } else if (knot_type == "equally_spaced") {
    knots <- seq(boundary[1], boundary[2], length.out = n_knots + 2L)
    # Remove boundary knots (first and last)
    knots <- knots[-c(1L, n_knots + 2L)]
  } else {
    # quantile-based
    knots <- as.numeric(stats::quantile(
      x,
      probs = seq(0, 1, length.out = n_knots + 2L)[-c(1L, n_knots + 2L)]
    ))
  }

  # ---- Basis matrix via splineDesign ----

  n_basis <- n_knots + degree + 1L

  if (n_basis < 1L) {
    stop("Number of basis functions (n_knots + degree - 1 = ", n_knots + degree - 1L,
         ") must be >= 1.", call. = FALSE)
  }

  if (penalty_order >= n_basis) {
    stop("`penalty_order` (", penalty_order, ") must be < n_basis (", n_basis,
         ").", call. = FALSE)
  }

  # Augmented knot vector: boundary repeated (degree + 1) times + internal knots
  knot_vec <- c(
    rep(boundary[1], degree + 1L),
    knots,
    rep(boundary[2], degree + 1L)
  )

  # Evaluate B-spline basis — splineDesign handles boundary values correctly
  B <- splines::splineDesign(
    knots = knot_vec,
    x = x,
    ord = degree + 1L,
    outer.ok = TRUE
  )

  # ---- Penalty matrix ----
  P <- make_penalty(n_basis, order = penalty_order)

  # ---- Assemble and return ----

  out <- list(
    B             = B,
    knots         = knots,
    boundary      = boundary,
    degree        = degree,
    penalty_order = penalty_order,
    P             = P,
    x             = x,
    n_basis       = n_basis,
    knot_type     = knot_type
  )
  class(out) <- "fda_basis"
  out
}


# ---- print.fda_basis ----

#' Print Method for fda_basis Objects
#'
#' Displays a compact summary of a B-spline basis object.
#'
#' @param x An object of class `"fda_basis"`.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @examples
#' basis <- bspline_basis(seq(0, 1, length.out = 50), n_knots = 8)
#' print(basis)
#'
#' @export
print.fda_basis <- function(x, ...) {
  cat("B-spline basis (fda_basis)\n")
  cat("  Degree:           ", x$degree, "\n")
  cat("  Internal knots:   ", length(x$knots), " (", x$knot_type, ")\n", sep = "")
  cat("  Basis functions:  ", x$n_basis, "\n")
  cat("  Eval. points:     ", length(x$x), "\n")
  cat("  Boundary:         [", format(x$boundary[1], digits = 4), ", ",
      format(x$boundary[2], digits = 4), "]\n", sep = "")
  cat("  Penalty order:    ", x$penalty_order, "\n")
  cat("  Basis matrix dim: ", nrow(x$B), " x ", ncol(x$B), "\n", sep = "")
  invisible(x)
}


# ---- plot.fda_basis ----

#' Plot Method for fda_basis Objects
#'
#' Plots all basis functions of a B-spline basis using base R graphics. Each
#' basis function is drawn in a distinct colour.
#'
#' @param x An object of class `"fda_basis"`.
#' @param col Character vector of colours for basis functions. Defaults to a
#'   rainbow palette with length equal to `n_basis`.
#' @param lwd Numeric; line width. Default is 1.5.
#' @param main Character; plot title. Defaults to an informative string.
#' @param xlab,ylab Character; axis labels.
#' @param ... Additional arguments passed to [graphics::matplot()].
#'
#' @return Invisibly returns the input object `x`.
#'
#' @examples
#' \donttest{
#' x <- seq(0, 1, length.out = 200)
#' basis <- bspline_basis(x, n_knots = 10, degree = 3)
#' plot(basis)
#' }
#'
#' @export
plot.fda_basis <- function(x,
                           col = NULL,
                           lwd = 1.5,
                           main = NULL,
                           xlab = "x",
                           ylab = "Basis function value",
                           ...) {

  if (is.null(col)) {
    col <- grDevices::colorRampPalette(
      c("#440154", "#31688E", "#35B779", "#FDE725")
    )(x$n_basis)
  }

  if (is.null(main)) {
    main <- paste0(
      "B-spline basis (degree = ", x$degree,
      ", ", length(x$knots), " knots, ",
      x$n_basis, " functions)"
    )
  }

  # Sort by evaluation points for clean line plots
  ord <- order(x$x)
  x_sorted <- x$x[ord]
  B_sorted <- x$B[ord, , drop = FALSE]

  graphics::matplot(
    x = x_sorted,
    y = B_sorted,
    type = "l",
    lty = 1,
    lwd = lwd,
    col = col,
    main = main,
    xlab = xlab,
    ylab = ylab,
    ...
  )

  # Add knot positions as tick marks on the x-axis
  graphics::rug(x$knots, col = "grey40", lwd = 0.8)

  invisible(x)
}
