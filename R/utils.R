# funcrop: Internal utility functions
# Not exported -- for package-internal use only.

# ---- Input validation helpers ------------------------------------------------

#' Validate a numeric vector (no NA, all finite)
#' @param x Vector to validate.
#' @param name Character name for error messages.
#' @return Invisible TRUE on success; raises an error otherwise.
#' @noRd
.validate_numeric <- function(x, name) {
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric, got %s.", name, class(x)[1L]),
         call. = FALSE)
  }
  if (anyNA(x)) {
    stop(sprintf("`%s` must not contain NA values.", name), call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate and coerce to factor (no NA)
#' @param x Vector to validate (character or factor).
#' @param name Character name for error messages.
#' @return A factor.
#' @noRd
.validate_factor <- function(x, name) {
  if (!(is.character(x) || is.factor(x))) {
    stop(sprintf("`%s` must be character or factor, got %s.", name,
                 class(x)[1L]), call. = FALSE)
  }
  if (anyNA(x)) {
    stop(sprintf("`%s` must not contain NA values.", name), call. = FALSE)
  }
  if (is.character(x)) x <- as.factor(x)
  x
}

#' Validate a positive integer scalar
#' @param x Value to validate.
#' @param name Character name for error messages.
#' @return Invisible TRUE on success.
#' @noRd
.validate_positive_integer <- function(x, name) {
  if (length(x) != 1L) {
    stop(sprintf("`%s` must be a scalar (length 1), got length %d.",
                 name, length(x)), call. = FALSE)
  }
  if (!is.numeric(x) || is.na(x) || x != as.integer(x) || x < 1L) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate engine choice and check availability
#' @param engine Character string: "asreml" or "bayesreml".
#' @return The validated engine string (invisibly).
#' @noRd
.validate_engine <- function(engine) {
  if (length(engine) != 1L || !is.character(engine)) {
    stop("`engine` must be a single character string.", call. = FALSE)
  }
  engine <- tolower(engine)
  valid_engines <- c("asreml", "bayesreml", "lme4", "mgcv")
  if (!engine %in% valid_engines) {
    stop(sprintf("`engine` must be one of: %s. Got '%s'.",
                 paste(valid_engines, collapse = ", "), engine),
         call. = FALSE)
  }
  # Check availability
  if (engine == "asreml" && !.has_asreml()) {
    stop("Engine 'asreml' requested but the asreml package is not installed.\n",
         "See https://vsni.co.uk/software/asreml-r for installation.",
         call. = FALSE)
  }
  if (engine == "bayesreml" && !.has_bayesreml()) {
    stop("Engine 'bayesreml' requested but the bayesreml package is not ",
         "installed.", call. = FALSE)
  }
  if (engine == "lme4" && !.has_lme4()) {
    stop("Engine 'lme4' requested but the lme4 package is not installed.\n",
         "Install via: install.packages('lme4')", call. = FALSE)
  }
  if (engine == "mgcv" && !.has_mgcv()) {
    stop("Engine 'mgcv' requested but the mgcv package is not installed.\n",
         "Install via: install.packages('mgcv')", call. = FALSE)
  }
  invisible(engine)
}

# ---- Package availability checks --------------------------------------------

#' Check if asreml is installed
#' @return Logical scalar.
#' @noRd
.has_asreml <- function() {
  requireNamespace("asreml", quietly = TRUE)
}

#' Check if bayesreml is installed
#' @return Logical scalar.
#' @noRd
.has_bayesreml <- function() {
  requireNamespace("bayesreml", quietly = TRUE)
}

#' Check if lme4 is installed
#' @return Logical scalar.
#' @noRd
.has_lme4 <- function() {
  requireNamespace("lme4", quietly = TRUE)
}

#' Check if mgcv is installed
#' @return Logical scalar.
#' @noRd
.has_mgcv <- function() {
  requireNamespace("mgcv", quietly = TRUE)
}

# ---- Messaging ---------------------------------------------------------------

#' Verbose message wrapper
#'
#' Prints a message only when `getOption("funcrop.verbose", TRUE)` is TRUE.
#' @param ... Arguments passed to [message()].
#' @return Invisible NULL.
#' @noRd
.msg <- function(...) {
  if (isTRUE(getOption("funcrop.verbose", default = TRUE))) {
    message(...)
  }
  invisible(NULL)
}

# ---- Linear algebra utilities ------------------------------------------------

#' Row-wise Kronecker product of two matrices
#'
#' For matrices A (n x p) and B (n x q), returns an (n x pq) matrix where
#' row i equals `kronecker(A[i, ], B[i, ])`, i.e., the outer product of the
#' two row vectors reshaped into a single row. This is the key operation for
#' constructing tensor-product basis matrices from marginal bases.
#'
#' Implementation avoids explicit per-row loops by using vectorised indexing.
#'
#' @param A Numeric matrix (n x p).
#' @param B Numeric matrix (n x q).
#' @return Numeric matrix (n x pq).
#' @noRd
.row_kronecker <- function(A, B) {
  # --- Input validation -------------------------------------------------------
  if (!is.matrix(A) || !is.numeric(A)) {
    stop("`A` must be a numeric matrix.", call. = FALSE)
  }

  if (!is.matrix(B) || !is.numeric(B)) {
    stop("`B` must be a numeric matrix.", call. = FALSE)
  }
  n_a <- nrow(A)
  n_b <- nrow(B)
  if (n_a != n_b) {
    stop(sprintf("`A` and `B` must have the same number of rows (%d vs %d).",
                 n_a, n_b), call. = FALSE)
  }

  n <- n_a
  p <- ncol(A)
  q <- ncol(B)

  # Vectorised: replicate columns of A q times, multiply element-wise with

  # B tiled p times. Each row i becomes: (A[i,1]*B[i,], A[i,2]*B[i,], ...).
  # A_rep: n x pq  -- each column j of A repeated q times side-by-side
  # B_rep: n x pq  -- columns of B tiled p times
  col_idx_a <- rep(seq_len(p), each = q)
  col_idx_b <- rep(seq_len(q), times = p)

  A[, col_idx_a, drop = FALSE] * B[, col_idx_b, drop = FALSE]
}

#' Construct a difference matrix of given order
#'
#' Creates a sparse difference matrix D of dimension (n - order) x n such that
#' D %*% x computes the order-th differences of a vector x of length n. Used in
#' P-spline penalty construction where the penalty is P = D'D.
#'
#' @param n Integer, the number of columns (e.g., number of B-spline
#'   coefficients).
#' @param order Integer, the order of differencing (default 2). Must satisfy
#'   `order < n`.
#' @return A sparse Matrix (class "dgCMatrix") of dimension (n - order) x n.
#' @noRd
.difference_matrix <- function(n, order = 2L) {
  .validate_positive_integer(n, "n")
  .validate_positive_integer(order, "order")
  n <- as.integer(n)
  order <- as.integer(order)

  if (order >= n) {
    stop(sprintf("`order` (%d) must be less than `n` (%d).", order, n),
         call. = FALSE)
  }

  # Start with the identity matrix (order 0), then iterentially difference
  D <- Matrix::Diagonal(n)
  for (k in seq_len(order)) {
    # Each differencing step: take rows 2:nrow minus rows 1:(nrow-1)
    m <- nrow(D)
    D <- D[2L:m, , drop = FALSE] - D[1L:(m - 1L), , drop = FALSE]
  }
  # Ensure compressed sparse column format
  as(D, "CsparseMatrix")
}
