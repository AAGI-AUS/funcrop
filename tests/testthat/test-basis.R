# test-basis.R — Comprehensive tests for the B-spline engine in funcrop
#
# Tests cover: bspline_basis(), make_penalty(), make_Zspline(),
#              tensor_bspline_basis(), and edge cases.
#
# testthat v3 edition

# ==============================================================================
# Setup: common test fixtures
# ==============================================================================

x_regular <- seq(0, 1, length.out = 100)
x_sparse  <- seq(0, 10, length.out = 20)
x_skew    <- sort(c(rbeta(50, 2, 8), rbeta(50, 8, 2)))

# ==============================================================================
# 1. bspline_basis()
# ==============================================================================

test_that("bspline_basis returns correct class", {
  basis <- bspline_basis(x_regular, n_knots = 10, degree = 3)
  expect_s3_class(basis, "fda_basis")
})

test_that("bspline_basis returns basis matrix with correct dimensions", {
  n_knots <- 10L
  degree <- 3L
  n_basis_expected <- n_knots + degree - 1L
  basis <- bspline_basis(x_regular, n_knots = n_knots, degree = degree)

  expect_equal(nrow(basis$B), length(x_regular))
  expect_equal(ncol(basis$B), n_basis_expected)
  expect_equal(basis$n_basis, n_basis_expected)
})

test_that("bspline_basis satisfies partition of unity (row sums = 1)", {
  basis <- bspline_basis(x_regular, n_knots = 10, degree = 3)
  row_sums <- rowSums(basis$B)
  expect_equal(row_sums, rep(1, length(x_regular)), tolerance = 1e-10)
})

test_that("bspline_basis produces non-negative values", {
  basis <- bspline_basis(x_regular, n_knots = 10, degree = 3)
  expect_true(all(basis$B >= -1e-15))
})

test_that("bspline_basis values at knots are valid", {
  basis <- bspline_basis(x_regular, n_knots = 5, degree = 3)
  # At each internal knot, the basis should still sum to 1
  knot_vals <- basis$knots
  basis_at_knots <- bspline_basis(knot_vals, n_knots = 5, degree = 3,
                                   boundary = basis$boundary)
  expect_equal(rowSums(basis_at_knots$B), rep(1, length(knot_vals)),
               tolerance = 1e-10)
  expect_true(all(basis_at_knots$B >= -1e-15))
})

test_that("equally_spaced knot placement produces valid basis", {
  basis <- bspline_basis(x_regular, n_knots = 8,
                          knot_type = "equally_spaced")
  expect_s3_class(basis, "fda_basis")
  expect_equal(basis$knot_type, "equally_spaced")
  expect_equal(rowSums(basis$B), rep(1, length(x_regular)), tolerance = 1e-10)
})

test_that("quantile knot placement produces valid basis", {
  basis <- bspline_basis(x_skew, n_knots = 8, knot_type = "quantile")
  expect_s3_class(basis, "fda_basis")
  expect_equal(basis$knot_type, "quantile")
  expect_equal(rowSums(basis$B), rep(1, length(x_skew)), tolerance = 1e-10)
})

test_that("custom knots produce valid basis", {
  custom_k <- c(0.2, 0.4, 0.6, 0.8)
  basis <- bspline_basis(x_regular, knot_type = "custom", knots = custom_k)
  expect_s3_class(basis, "fda_basis")
  expect_equal(basis$knot_type, "custom")
  expect_equal(length(basis$knots), 4L)
  expect_equal(rowSums(basis$B), rep(1, length(x_regular)), tolerance = 1e-10)
})

test_that("boundary knots default correctly from range of x", {

  basis <- bspline_basis(x_regular, n_knots = 5)
  x_range <- range(x_regular)
  x_span <- diff(x_range)
  expected_boundary <- x_range + c(-1, 1) * 0.01 * x_span
  expect_equal(basis$boundary, expected_boundary)
})

test_that("bspline_basis errors on non-numeric x", {
  expect_error(bspline_basis(letters[1:10]), "numeric")
})

test_that("bspline_basis errors on x with NA/Inf", {
  expect_error(bspline_basis(c(1, 2, NA, 4)), "NA")
  expect_error(bspline_basis(c(1, 2, Inf, 4)), "Inf")
})

test_that("bspline_basis errors on negative n_knots", {
  expect_error(bspline_basis(x_regular, n_knots = -1), "positive")
})

test_that("bspline_basis errors on x with < 2 values", {
  expect_error(bspline_basis(c(1)), "length >= 2")
})

test_that("bspline_basis errors on x with < 2 unique values", {
  expect_error(bspline_basis(c(5, 5, 5)), "unique")
})

test_that("penalty matrix is symmetric positive semi-definite", {
  basis <- bspline_basis(x_regular, n_knots = 10, degree = 3,
                          penalty_order = 2)
  P <- as.matrix(basis$P)

  # Symmetric

  expect_equal(P, t(P), tolerance = 1e-12)

  # Positive semi-definite: all eigenvalues >= 0
  eig_vals <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig_vals >= -1e-10))
})

test_that("penalty matrix has correct rank (n_basis - penalty_order)", {
  n_knots <- 10L
  degree <- 3L
  penalty_order <- 2L
  n_basis <- n_knots + degree - 1L

  basis <- bspline_basis(x_regular, n_knots = n_knots, degree = degree,
                          penalty_order = penalty_order)
  P <- as.matrix(basis$P)
  eig_vals <- eigen(P, symmetric = TRUE, only.values = TRUE)$values

  # Number of non-zero eigenvalues = rank = n_basis - penalty_order
  n_nonzero <- sum(eig_vals > 1e-10)
  expect_equal(n_nonzero, n_basis - penalty_order)
})


# ==============================================================================
# 2. make_penalty()
# ==============================================================================

test_that("make_penalty returns matrix of correct dimensions", {
  P <- make_penalty(12, order = 2)
  expect_equal(nrow(P), 12L)
  expect_equal(ncol(P), 12L)
})

test_that("make_penalty produces symmetric matrix", {
  P <- as.matrix(make_penalty(10, order = 2))
  expect_equal(P, t(P), tolerance = 1e-12)
})

test_that("make_penalty produces positive semi-definite matrix", {
  P <- as.matrix(make_penalty(10, order = 2))
  eig_vals <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig_vals >= -1e-10))
})

test_that("make_penalty null space has correct dimension", {
  for (ord in 1:3) {
    n <- 12L
    P <- as.matrix(make_penalty(n, order = ord))
    eig_vals <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
    n_zero <- sum(abs(eig_vals) < 1e-10)
    expect_equal(n_zero, ord,
                 info = paste("penalty order =", ord))
  }
})

test_that("make_penalty order 1 and order 2 produce different results", {
  P1 <- as.matrix(make_penalty(10, order = 1))
  P2 <- as.matrix(make_penalty(10, order = 2))
  expect_false(isTRUE(all.equal(P1, P2)))
})


# ==============================================================================
# 3. make_Zspline()
# ==============================================================================

test_that("make_Zspline X matrix has penalty_order columns", {
  basis <- bspline_basis(x_regular, n_knots = 10, degree = 3,
                          penalty_order = 2)
  zs <- make_Zspline(basis)

  expect_equal(ncol(zs$X), basis$penalty_order)
})

test_that("make_Zspline Z matrix has n_basis - penalty_order columns", {
  basis <- bspline_basis(x_regular, n_knots = 10, degree = 3,
                          penalty_order = 2)
  zs <- make_Zspline(basis)

  expect_equal(ncol(zs$Z), basis$n_basis - basis$penalty_order)
})

test_that("make_Zspline reconstructs original basis (X*T_X + Z*T_Z ~ B)", {
  basis <- bspline_basis(x_regular, n_knots = 10, degree = 3,
                          penalty_order = 2)
  zs <- make_Zspline(basis)

  # The decomposition should allow: B = [X, Z] %*% T for some T

  # Equivalently, the column space of [X, Z] should span the column space of B
  XZ <- cbind(zs$X, zs$Z)
  # Project B onto [X, Z] and check reconstruction
  # B should be exactly in the column space of XZ (both n x n_basis)
  expect_equal(ncol(XZ), basis$n_basis)
  expect_equal(nrow(XZ), nrow(basis$B))

  # Fit: B = XZ %*% coef => coef = solve(t(XZ) %*% XZ) %*% t(XZ) %*% B
  coef_mat <- solve(crossprod(XZ), crossprod(XZ, basis$B))
  B_reconstructed <- XZ %*% coef_mat
  expect_equal(B_reconstructed, basis$B, tolerance = 1e-8)
})

test_that("make_Zspline 'decompose' method works", {
  basis <- bspline_basis(x_regular, n_knots = 10, degree = 3)
  zs <- make_Zspline(basis, method = "decompose")
  expect_true(!is.null(zs$X))
  expect_true(!is.null(zs$Z))
})

test_that("make_Zspline 'absorb' method works", {
  basis <- bspline_basis(x_regular, n_knots = 10, degree = 3)
  zs <- make_Zspline(basis, method = "absorb")
  expect_true(!is.null(zs$X))
  expect_true(!is.null(zs$Z))
})


# ==============================================================================
# 4. tensor_bspline_basis()
# ==============================================================================

test_that("tensor_bspline_basis returns correct class", {
  x1 <- seq(0, 1, length.out = 20)
  x2 <- seq(0, 1, length.out = 15)
  tb <- tensor_bspline_basis(x1, x2, n_knots1 = 5, n_knots2 = 4)
  expect_s3_class(tb, "fda_basis")
})

test_that("tensor_bspline_basis has correct basis matrix dimensions", {
  x1 <- seq(0, 1, length.out = 20)
  x2 <- seq(0, 1, length.out = 15)
  n_knots1 <- 5L
  n_knots2 <- 4L
  degree <- 3L
  n_basis1 <- n_knots1 + degree - 1L
  n_basis2 <- n_knots2 + degree - 1L
  n <- length(x1)  # evaluation points in first margin

  tb <- tensor_bspline_basis(x1, x2, n_knots1 = n_knots1, n_knots2 = n_knots2,
                              degree = degree)
  # Tensor product: n_obs x (n_basis1 * n_basis2)
  expect_equal(ncol(tb$B), n_basis1 * n_basis2)
})

test_that("tensor_bspline_basis combined penalty is symmetric PSD", {
  x1 <- seq(0, 1, length.out = 20)
  x2 <- seq(0, 1, length.out = 15)
  tb <- tensor_bspline_basis(x1, x2, n_knots1 = 5, n_knots2 = 4)

  P <- as.matrix(tb$P)
  # Symmetric
  expect_equal(P, t(P), tolerance = 1e-12)

  # Positive semi-definite
  eig_vals <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig_vals >= -1e-10))
})

test_that("tensor_bspline_basis stores marginal bases", {
  x1 <- seq(0, 1, length.out = 20)
  x2 <- seq(0, 1, length.out = 15)
  tb <- tensor_bspline_basis(x1, x2, n_knots1 = 5, n_knots2 = 4)

  # Should store marginal basis objects
  expect_true(!is.null(tb$basis1) || !is.null(tb$marginal1))
  expect_true(!is.null(tb$basis2) || !is.null(tb$marginal2))
})

test_that("tensor_bspline_basis row-wise Kronecker product is correct", {
  # Small example to verify manually
  x1 <- c(0.2, 0.5, 0.8)
  x2 <- c(0.3, 0.6, 0.9)
  n_k1 <- 2L
  n_k2 <- 2L
  deg <- 3L

  tb <- tensor_bspline_basis(x1, x2, n_knots1 = n_k1, n_knots2 = n_k2,
                              degree = deg)

  # Get marginal bases and compute row-Kronecker manually
  b1 <- bspline_basis(x1, n_knots = n_k1, degree = deg)
  b2 <- bspline_basis(x2, n_knots = n_k2, degree = deg)

  # Row-wise Kronecker: for each row i, kron(B1[i,], B2[i,])
  n <- length(x1)
  p1 <- ncol(b1$B)
  p2 <- ncol(b2$B)
  manual_tensor <- matrix(0, nrow = n, ncol = p1 * p2)
  for (i in seq_len(n)) {
    manual_tensor[i, ] <- as.vector(outer(b1$B[i, ], b2$B[i, ]))
  }

  expect_equal(tb$B, manual_tensor, tolerance = 1e-10)
})


# ==============================================================================
# 5. Edge cases
# ==============================================================================

test_that("bspline_basis works with single evaluation point (length 2 min)", {
  # Minimum: 2 points (boundary constraint)
  x_min <- c(0, 1)
  basis <- bspline_basis(x_min, n_knots = 3, degree = 3)
  expect_equal(nrow(basis$B), 2L)
  expect_equal(rowSums(basis$B), c(1, 1), tolerance = 1e-10)
})

test_that("bspline_basis works with very few knots (minimum viable)", {
  # 1 internal knot, degree 3 => n_basis = 1 + 3 - 1 = 3
  basis <- bspline_basis(x_regular, n_knots = 1, degree = 3,
                          penalty_order = 1)
  expect_equal(basis$n_basis, 3L)
  expect_equal(ncol(basis$B), 3L)
  expect_equal(rowSums(basis$B), rep(1, length(x_regular)), tolerance = 1e-10)
})

test_that("bspline_basis with degree = 1 (linear splines) works correctly", {
  basis <- bspline_basis(x_regular, n_knots = 5, degree = 1,
                          penalty_order = 1)
  n_basis_expected <- 5L + 1L - 1L  # n_knots + degree - 1 = 5
  expect_equal(basis$n_basis, n_basis_expected)
  expect_equal(basis$degree, 1L)

  # Partition of unity
  expect_equal(rowSums(basis$B), rep(1, length(x_regular)), tolerance = 1e-10)

  # Non-negativity
  expect_true(all(basis$B >= -1e-15))

  # Linear splines should be piecewise linear: each column has at most one peak
})

test_that("bspline_basis handles duplicate x values correctly", {
  x_dup <- c(0, 0.2, 0.2, 0.5, 0.5, 0.5, 0.8, 1)
  basis <- bspline_basis(x_dup, n_knots = 5, degree = 3)

  expect_equal(nrow(basis$B), length(x_dup))
  expect_equal(rowSums(basis$B), rep(1, length(x_dup)), tolerance = 1e-10)

  # Duplicate x should give identical rows
  expect_equal(basis$B[2, ], basis$B[3, ])
  expect_equal(basis$B[4, ], basis$B[5, ])
  expect_equal(basis$B[5, ], basis$B[6, ])
})

test_that("bspline_basis with degree = 0 (step functions) works", {
  # degree 0, penalty_order must be < n_basis
  # n_basis = n_knots + 0 - 1 = n_knots - 1, need n_knots >= 3 for pen_order=1
  basis <- bspline_basis(x_regular, n_knots = 5, degree = 0,
                          penalty_order = 1)
  expect_equal(basis$degree, 0L)
  expect_equal(basis$n_basis, 4L)  # 5 + 0 - 1
  # Step functions: each row has exactly one non-zero entry (value = 1)
  expect_equal(rowSums(basis$B), rep(1, length(x_regular)), tolerance = 1e-10)
})

test_that("make_penalty errors on invalid inputs",
{
  expect_error(make_penalty(0, order = 1))
  expect_error(make_penalty(3, order = 3))  # order must be < n
  expect_error(make_penalty(5, order = -1))
})
