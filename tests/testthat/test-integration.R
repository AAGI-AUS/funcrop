# test-integration.R — Integration tests for funcrop pipeline
#
# Tests the end-to-end workflow from basis construction through data
# preparation, and validates mathematical properties. Does NOT require
# ASReml or bayesreml — uses mocked backends or tests only backend-free
# components.
#
# testthat v3 edition

# ==============================================================================
# Setup: simulated data fixtures
# ==============================================================================

set.seed(42L)

# Simulated functional data: 10 varieties, 8 time points, known logistic curves
n_var    <- 10L
n_time   <- 8L
n_plots  <- n_var * 3L  # 3 reps
time_pts <- seq(10, 45, length.out = n_time)

# Known parameters per variety
wmax   <- runif(n_var, 20, 40)
rate   <- runif(n_var, 0.15, 0.35)
t_mid  <- runif(n_var, 22, 33)

# Generate data
sim_rows <- vector("list", n_plots)
for (p in seq_len(n_plots)) {
  v <- ((p - 1L) %% n_var) + 1L
  curve <- wmax[v] / (1 + exp(-rate[v] * (time_pts - t_mid[v])))
  sim_rows[[p]] <- data.table::data.table(
    plot_id = paste0("P", sprintf("%03d", p)),
    variety = paste0("V", sprintf("%02d", v)),
    time    = time_pts,
    value   = curve + rnorm(n_time, sd = 1.5)
  )
}
sim_dt <- data.table::rbindlist(sim_rows)


# ==============================================================================
# 1. Full pipeline test: basis -> data -> profiles structure
# ==============================================================================

test_that("full pipeline: basis construction through fda_data preparation", {
  # Step 1: construct basis
  basis <- bspline_basis(time_pts, n_knots = 4L, degree = 3L)
  expect_s3_class(basis, "fda_basis")
  expect_equal(nrow(basis$B), n_time)

  # Step 2: construct fda_data
  fd <- fda_data(
    time  = sim_dt$time,
    value = sim_dt$value,
    id    = sim_dt$plot_id,
    group = sim_dt$variety
  )
  expect_s3_class(fd, "fda_data")

  meta <- attr(fd, "fda_meta")
  expect_equal(meta$n_ids, n_plots)
  expect_equal(meta$n_timepoints, n_time)
  expect_true(meta$has_group)

  # Step 3: verify basis evaluated at data time points matches
  basis_at_data <- bspline_basis(
    sort(unique(sim_dt$time)),
    n_knots  = 4L,
    degree   = 3L,
    boundary = basis$boundary
  )
  expect_equal(nrow(basis_at_data$B), n_time)
  expect_equal(ncol(basis_at_data$B), basis$n_basis)

  # Step 4: mixed model reparameterisation
  mm <- make_Zspline(basis)
  expect_equal(ncol(mm$X) + ncol(mm$Z), basis$n_basis)
  expect_equal(nrow(mm$X), n_time)
  expect_equal(nrow(mm$Z), n_time)
})


# ==============================================================================
# 2. B-spline inner product matrix properties
# ==============================================================================

test_that("B-spline inner product matrix J is symmetric", {
  x_fine <- seq(0, 1, length.out = 500)
  basis  <- bspline_basis(x_fine, n_knots = 8L, degree = 3L)
  B <- basis$B

  # Approximate inner product J = (1/n) * t(B) %*% B * domain_width
  domain_width <- diff(basis$boundary)
  dx <- domain_width / (length(x_fine) - 1L)
  J <- crossprod(B) * dx

  expect_equal(J, t(J), tolerance = 1e-12)
})

test_that("B-spline inner product matrix J is positive definite", {
  x_fine <- seq(0, 1, length.out = 500)
  basis  <- bspline_basis(x_fine, n_knots = 8L, degree = 3L)
  B <- basis$B

  domain_width <- diff(basis$boundary)
  dx <- domain_width / (length(x_fine) - 1L)
  J <- crossprod(B) * dx

  eig_vals <- eigen(J, symmetric = TRUE, only.values = TRUE)$values
  # All eigenvalues should be positive (not just non-negative)
  expect_true(all(eig_vals > 1e-10))
})

test_that("B-spline partition of unity: sum of integrals equals domain width", {
  x_fine <- seq(0, 1, length.out = 1000)
  basis  <- bspline_basis(x_fine, n_knots = 8L, degree = 3L)
  B <- basis$B

  domain_width <- diff(basis$boundary)
  dx <- domain_width / (length(x_fine) - 1L)

  # Sum of all basis function integrals should equal the domain width
  # integral of sum_j B_j(x) dx = integral of 1 dx = domain_width
  # (by partition of unity)
  total_integral <- sum(colSums(B) * dx)
  expect_equal(total_integral, domain_width, tolerance = 0.01)
})


# ==============================================================================
# 3. Functional covariate computation
# ==============================================================================

test_that("functional covariate C = Alpha %*% J has correct dimensions", {
  x_fine <- seq(0, 1, length.out = 200)
  basis  <- bspline_basis(x_fine, n_knots = 6L, degree = 3L)
  n_basis <- basis$n_basis

  # Simulate coefficient matrix Alpha (n_ids x n_basis)
  n_ids <- 15L
  set.seed(123)
  Alpha <- matrix(rnorm(n_ids * n_basis), nrow = n_ids, ncol = n_basis)

  # Inner product matrix
  dx <- diff(basis$boundary) / (length(x_fine) - 1L)
  J <- crossprod(basis$B) * dx

  # Functional covariate: C = Alpha %*% J

  C <- Alpha %*% J
  expect_equal(dim(C), c(n_ids, n_basis))

  # Each row of C is the functional covariate for one individual
  # C should not be all zeros if Alpha is non-zero
  expect_true(any(abs(C) > 1e-10))
})

test_that("functional covariate with known Alpha is mathematically correct", {
  x_fine <- seq(0, 1, length.out = 300)
  basis  <- bspline_basis(x_fine, n_knots = 5L, degree = 3L)
  B <- basis$B
  n_basis <- basis$n_basis

  # Known Alpha: identity (each individual uses one basis function)
  Alpha <- diag(n_basis)

  dx <- diff(basis$boundary) / (length(x_fine) - 1L)
  J <- crossprod(B) * dx
  C <- Alpha %*% J

  # C should equal J itself when Alpha = I
  expect_equal(C, J, tolerance = 1e-10)
})


# ==============================================================================
# 4. Genomic relationship matrix — SKIPPED if functions not yet implemented
# ==============================================================================

test_that("make_genomic_matrix produces symmetric matrix on test data", {
  skip_if_not(
    exists("make_genomic_matrix", envir = asNamespace("funcrop"),
           mode = "function"),
    message = "make_genomic_matrix not yet implemented"
  )

  # Small test marker data: 10 individuals, 50 markers
  set.seed(99)
  M <- matrix(sample(0:2, 500, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
              nrow = 10, ncol = 50)
  rownames(M) <- paste0("G", 1:10)

  G <- make_genomic_matrix(M)
  expect_true(isSymmetric(G, tol = 1e-10))
  expect_equal(nrow(G), 10L)
  expect_equal(ncol(G), 10L)
})

test_that("make_genomic_matrix mean diagonal approx 1 when scaled", {
  skip_if_not(
    exists("make_genomic_matrix", envir = asNamespace("funcrop"),
           mode = "function"),
    message = "make_genomic_matrix not yet implemented"
  )

  set.seed(99)
  M <- matrix(sample(0:2, 500, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
              nrow = 10, ncol = 50)
  rownames(M) <- paste0("G", 1:10)

  G <- make_genomic_matrix(M, scale = TRUE)
  expect_equal(mean(diag(G)), 1.0, tolerance = 0.3)
})

test_that("make_genomic_matrix methods produce different results", {
  skip_if_not(
    exists("make_genomic_matrix", envir = asNamespace("funcrop"),
           mode = "function"),
    message = "make_genomic_matrix not yet implemented"
  )

  set.seed(99)
  M <- matrix(sample(0:2, 500, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
              nrow = 10, ncol = 50)
  rownames(M) <- paste0("G", 1:10)

  G1 <- make_genomic_matrix(M, method = "vanraden1")
  G2 <- make_genomic_matrix(M, method = "vanraden2")
  expect_false(isTRUE(all.equal(G1, G2)))
})

test_that("make_genomic_matrix MAF filter works", {
  skip_if_not(
    exists("make_genomic_matrix", envir = asNamespace("funcrop"),
           mode = "function"),
    message = "make_genomic_matrix not yet implemented"
  )

  set.seed(99)
  # Include some nearly monomorphic markers
  M <- matrix(sample(0:2, 500, replace = TRUE, prob = c(0.01, 0.04, 0.95)),
              nrow = 10, ncol = 50)
  rownames(M) <- paste0("G", 1:10)

  # With stringent MAF filter should use fewer markers
  G_strict <- make_genomic_matrix(M, min_maf = 0.1)
  G_loose  <- make_genomic_matrix(M, min_maf = 0.01)
  # Both should still be valid symmetric matrices
  expect_true(isSymmetric(G_strict, tol = 1e-10))
  expect_true(isSymmetric(G_loose,  tol = 1e-10))
})


# ==============================================================================
# 5. Pedigree relationship matrix — SKIPPED if functions not yet implemented
# ==============================================================================

test_that("make_pedigree_matrix produces valid A matrix for simple pedigree", {
  skip_if_not(
    exists("make_pedigree_matrix", envir = asNamespace("funcrop"),
           mode = "function"),
    message = "make_pedigree_matrix not yet implemented"
  )

  # Simple 3-generation pedigree
  ped <- data.frame(
    id     = c("F1", "M1", "O1", "O2"),
    sire   = c(NA,   NA,   "F1", "F1"),
    dam    = c(NA,   NA,   "M1", "M1"),
    stringsAsFactors = FALSE
  )

  A <- make_pedigree_matrix(ped)

  # Diagonal elements >= 1 (1 + inbreeding coefficient)
  expect_true(all(diag(A) >= 1))

  # Symmetric
  expect_true(isSymmetric(A, tol = 1e-10))

  # Dimensions
  expect_equal(nrow(A), 4L)
  expect_equal(ncol(A), 4L)

  # Known values:
  # Parent-offspring: A[F1, O1] = 0.5
  expect_equal(A["F1", "O1"], 0.5, tolerance = 1e-10)
  # Full sibs: A[O1, O2] = 0.5
  expect_equal(A["O1", "O2"], 0.5, tolerance = 1e-10)
  # Founders unrelated: A[F1, M1] = 0
  expect_equal(A["F1", "M1"], 0, tolerance = 1e-10)
})

test_that("make_pedigree_matrix diagonal for non-inbred equals 1", {
  skip_if_not(
    exists("make_pedigree_matrix", envir = asNamespace("funcrop"),
           mode = "function"),
    message = "make_pedigree_matrix not yet implemented"
  )

  ped <- data.frame(
    id   = c("A", "B", "C"),
    sire = c(NA,  NA,  "A"),
    dam  = c(NA,  NA,  "B"),
    stringsAsFactors = FALSE
  )
  A <- make_pedigree_matrix(ped)
  # Non-inbred: all diagonals exactly 1
  expect_equal(unname(diag(A)), rep(1, 3), tolerance = 1e-10)
})


# ==============================================================================
# 6. check_relationship_matrix — SKIPPED if not yet implemented
# ==============================================================================

test_that("check_relationship_matrix passes PD matrix unchanged", {
  skip_if_not(
    exists("check_relationship_matrix", envir = asNamespace("funcrop"),
           mode = "function"),
    message = "check_relationship_matrix not yet implemented"
  )

  # Construct a known PD matrix
  set.seed(77)
  n <- 5L
  L <- matrix(rnorm(n * n), n, n)
  A <- crossprod(L)  # guaranteed PD
  A <- (A + t(A)) / 2  # ensure exact symmetry

  result <- check_relationship_matrix(A)
  expect_equal(result, A, tolerance = 1e-10)
})

test_that("check_relationship_matrix bends non-PD matrix", {
  skip_if_not(
    exists("check_relationship_matrix", envir = asNamespace("funcrop"),
           mode = "function"),
    message = "check_relationship_matrix not yet implemented"
  )

  # Construct a non-PD symmetric matrix
  A <- matrix(c(1, 0.9, 0.9,
                0.9, 1, 0.9,
                0.9, 0.9, -0.5), 3, 3)
  A <- (A + t(A)) / 2

  result <- check_relationship_matrix(A)

  # Result should be PD
  eig_vals <- eigen(result, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig_vals > 0))
  # Result should be symmetric
  expect_true(isSymmetric(result, tol = 1e-10))
})

test_that("check_relationship_matrix enforces symmetry", {
  skip_if_not(
    exists("check_relationship_matrix", envir = asNamespace("funcrop"),
           mode = "function"),
    message = "check_relationship_matrix not yet implemented"
  )

  # Slightly asymmetric PD matrix
  set.seed(77)
  n <- 4L
  L <- matrix(rnorm(n * n), n, n)
  A <- crossprod(L)
  A[1, 2] <- A[1, 2] + 0.001  # break symmetry slightly

  result <- check_relationship_matrix(A)
  expect_true(isSymmetric(result, tol = 1e-10))
})


# ==============================================================================
# 7. Backend dispatcher
# ==============================================================================

test_that("funcrop_engines returns character vector", {
  engines <- funcrop_engines()
  expect_type(engines, "character")
  # Could be empty if no engine installed, but must be character
  if (length(engines) > 0L) {
    expect_true(all(engines %in% c("mgcv", "lme4", "asreml", "bayesreml")))
  }
})

test_that(".resolve_engine('auto') returns valid engine or errors gracefully", {
  # Mock all engines as unavailable, expect informative error
  local_mocked_bindings(
    .has_asreml    = function() FALSE,
    .has_bayesreml = function() FALSE,
    .has_lme4      = function() FALSE,
    .has_mgcv      = function() FALSE,
    .package = "funcrop"
  )

  withr::local_options(list(funcrop.engine = NULL))

  expect_error(
    .resolve_engine("auto"),
    "No estimation engine"
  )
})

test_that(".resolve_engine('auto') prefers mgcv (open-source first)", {
  local_mocked_bindings(
    .has_asreml    = function() TRUE,
    .has_bayesreml = function() FALSE,
    .has_lme4      = function() TRUE,
    .has_mgcv      = function() TRUE,
    .package = "funcrop"
  )

  withr::local_options(list(funcrop.engine = NULL))

  engine <- .resolve_engine("auto")
  expect_equal(engine, "mgcv")
})

test_that(".resolve_engine('auto') falls back to lme4 then asreml", {
  local_mocked_bindings(
    .has_asreml    = function() TRUE,
    .has_bayesreml = function() FALSE,
    .has_lme4      = function() FALSE,
    .has_mgcv      = function() FALSE,
    .package = "funcrop"
  )

  withr::local_options(list(funcrop.engine = NULL))

  engine <- .resolve_engine("auto")
  expect_equal(engine, "asreml")
})

test_that(".resolve_engine('auto') falls back to bayesreml", {
  local_mocked_bindings(
    .has_asreml    = function() FALSE,
    .has_bayesreml = function() TRUE,
    .has_lme4      = function() FALSE,
    .has_mgcv      = function() FALSE,
    .package = "funcrop"
  )

  withr::local_options(list(funcrop.engine = NULL))

  engine <- .resolve_engine("auto")
  expect_equal(engine, "bayesreml")
})

test_that("funcrop_default_engine set/get works with mocked backend", {
  local_mocked_bindings(
    .has_asreml    = function() FALSE,
    .has_bayesreml = function() TRUE,
    .has_lme4      = function() FALSE,
    .has_mgcv      = function() FALSE,
    .package = "funcrop"
  )

  withr::local_options(list(funcrop.engine = NULL, funcrop.verbose = FALSE))

  # Set bayesreml as default
  result <- funcrop_default_engine("bayesreml")

  # Query current default
  current <- funcrop_default_engine()
  expect_equal(current, "bayesreml")
})

test_that("funcrop_default_engine errors on unavailable engine", {
  local_mocked_bindings(
    .has_asreml    = function() FALSE,
    .has_bayesreml = function() FALSE,
    .has_lme4      = function() FALSE,
    .has_mgcv      = function() FALSE,
    .package = "funcrop"
  )

  withr::local_options(list(funcrop.engine = NULL))

  expect_error(
    funcrop_default_engine("asreml"),
    "not installed"
  )
})

test_that("funcrop_engines lists all available backends", {
  local_mocked_bindings(
    .has_asreml    = function() TRUE,
    .has_bayesreml = function() TRUE,
    .has_lme4      = function() TRUE,
    .has_mgcv      = function() TRUE,
    .package = "funcrop"
  )

  engines <- funcrop_engines()
  expect_true("mgcv" %in% engines)
  expect_true("lme4" %in% engines)
  expect_true("asreml" %in% engines)
  expect_true("bayesreml" %in% engines)
})


# ==============================================================================
# 8. Edge cases
# ==============================================================================

test_that("fda_data with single variety in pipeline context", {
  one_plot <- sim_dt[variety == "V01"][1:n_time]
  fd <- fda_data(
    time  = one_plot$time,
    value = one_plot$value,
    id    = one_plot$plot_id,
    group = one_plot$variety
  )
  expect_s3_class(fd, "fda_data")
  meta <- attr(fd, "fda_meta")
  expect_equal(meta$n_ids, 1L)
})

test_that("basis with minimum knots (n_knots = 1) works in pipeline", {
  basis <- bspline_basis(time_pts, n_knots = 1L, degree = 3L,
                          penalty_order = 1L)
  expect_s3_class(basis, "fda_basis")
  # n_basis = 1 + 3 + 1 = 5
  expect_equal(basis$n_basis, 5L)

  # Still satisfies partition of unity
  expect_equal(rowSums(basis$B), rep(1, n_time), tolerance = 1e-10)

  # Can still reparameterise
  mm <- make_Zspline(basis)
  expect_equal(ncol(mm$X), 1L)  # null_dim = penalty_order = 1
  expect_equal(ncol(mm$Z), 4L)  # rank = n_basis - 1 = 4
})

test_that("tensor basis with minimum knots works", {
  x1 <- seq(0, 1, length.out = 10)
  x2 <- seq(0, 1, length.out = 10)
  tb <- tensor_bspline_basis(x1, x2, n_knots1 = 1L, n_knots2 = 1L,
                              degree1 = 3L, degree2 = 3L)
  expect_s3_class(tb, "fda_tensor_basis")
  # Each marginal: 1 + 3 + 1 = 5 basis functions => tensor: 25
  expect_equal(tb$n_basis, 25L)
  expect_equal(ncol(tb$B), 25L)
})

test_that("penalty matrix for small n_basis (n = 2, order = 1)", {
  P <- make_penalty(2L, order = 1L)
  expect_equal(nrow(P), 2L)
  expect_equal(ncol(P), 2L)
  P_dense <- as.matrix(P)
  # D = [-1, 1], P = D'D = [[1,-1],[-1,1]]
  expect_equal(P_dense, matrix(c(1, -1, -1, 1), 2, 2), tolerance = 1e-12)
})

test_that("make_Zspline decompose and absorb give same column spaces", {
  x <- seq(0, 1, length.out = 50)
  basis <- bspline_basis(x, n_knots = 6L, degree = 3L)

  mm_dec <- make_Zspline(basis, constraint = "decompose")
  mm_abs <- make_Zspline(basis, constraint = "absorb")

  # Both should have same X and Z dimensions

  expect_equal(ncol(mm_dec$X), ncol(mm_abs$X))
  expect_equal(ncol(mm_dec$Z), ncol(mm_abs$Z))

  # Column space of [X, Z] should span same space for both methods
  # Test: B can be reconstructed from either
  XZ_dec <- cbind(mm_dec$X, mm_dec$Z)
  XZ_abs <- cbind(mm_abs$X, mm_abs$Z)

  coef_dec <- solve(crossprod(XZ_dec), crossprod(XZ_dec, basis$B))
  coef_abs <- solve(crossprod(XZ_abs), crossprod(XZ_abs, basis$B))

  B_rec_dec <- XZ_dec %*% coef_dec
  B_rec_abs <- XZ_abs %*% coef_abs

  expect_equal(B_rec_dec, basis$B, tolerance = 1e-8)
  expect_equal(B_rec_abs, basis$B, tolerance = 1e-8)
})

test_that(".row_kronecker produces correct output for small matrices", {
  A <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  B <- matrix(c(5, 6, 7, 8), nrow = 2, ncol = 2)

  result <- funcrop:::.row_kronecker(A, B)

  # Row 1: kron(c(1,3), c(5,7)) = c(1*5, 1*7, 3*5, 3*7) = c(5,7,15,21)
  expect_equal(result[1, ], c(5, 7, 15, 21))
  # Row 2: kron(c(2,4), c(6,8)) = c(2*6, 2*8, 4*6, 4*8) = c(12,16,24,32)
  expect_equal(result[2, ], c(12, 16, 24, 32))
})

test_that(".validate_numeric rejects non-numeric input", {
  expect_error(
    funcrop:::.validate_numeric("hello", "test_var"),
    "numeric"
  )
})

test_that(".validate_numeric rejects NA values", {
  expect_error(
    funcrop:::.validate_numeric(c(1, NA, 3), "test_var"),
    "NA"
  )
})

test_that(".validate_numeric rejects Inf values", {
  expect_error(
    funcrop:::.validate_numeric(c(1, Inf), "test_var"),
    "finite"
  )
})

test_that(".validate_numeric accepts valid input", {
  expect_invisible(funcrop:::.validate_numeric(c(1, 2, 3), "test_var"))
})
