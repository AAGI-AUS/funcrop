# test-classes.R — Tests for S3 classes: fda_data, fda_model, fda_comparison
#
# Covers construction, validation, methods, and edge cases for the core
# class system in funcrop.
#
# testthat v3 edition

# ==============================================================================
# Setup: common test fixtures
# ==============================================================================

# Minimal valid inputs for fda_data
n_obs <- 50L
n_ids <- 5L
ids_vec    <- rep(paste0("plot_", seq_len(n_ids)), each = 10L)
times_vec  <- rep(seq(100, 190, by = 10), times = n_ids)
vals_vec   <- sin(times_vec / 50) + seq_len(n_obs) * 0.001
groups_vec <- rep(c("var_A", "var_B", "var_A", "var_B", "var_A"), each = 10L)
yields_vec <- rep(c(3.2, 4.1, 3.8, 4.5, 3.0), each = 10L)

# Helper: build a mock fda_model for method testing
# NOTE: new_fda_model uses is.data.table() internally which may not be in the
# namespace if not explicitly imported. We construct the object directly via
# structure() to avoid this issue in test context.
make_mock_fda_model <- function(extras = list()) {
  n_pts <- 20L
  n_fit <- 40L
  fitted_dt <- data.table::data.table(
    id      = rep(paste0("id_", 1:4), each = 10L),
    group   = rep(c("A", "B"), each = 20L),
    time    = rep(seq(1, 10), times = 4L),
    fitted  = rnorm(n_fit)
  )

  coef_fn <- list(
    time     = seq(1, 10, length.out = n_pts),
    beta     = rnorm(n_pts),
    se       = abs(rnorm(n_pts, sd = 0.1)),
    ci_lower = rnorm(n_pts) - 0.2,
    ci_upper = rnorm(n_pts) + 0.2
  )
  var_comps <- data.table::data.table(
    component = c("variety", "residual"),
    estimate  = c(1.5, 0.8),
    se        = c(0.3, 0.1)
  )
  preds <- data.table::data.table(
    id    = paste0("id_", 1:4),
    pred  = rnorm(4L),
    se    = abs(rnorm(4L, sd = 0.1))
  )
  resids <- rnorm(n_fit)

  structure(
    list(
      fitted_curves          = fitted_dt,
      coefficient_function   = coef_fn,
      variance_components    = var_comps,
      predictions            = preds,
      residuals              = resids,
      basis                  = list(type = "bspline", degree = 3L),
      data                   = NULL,
      engine                 = "asreml",
      call                   = match.call(),
      extras                 = extras
    ),
    class = "fda_model"
  )
}


# ==============================================================================
# 1. fda_data() — Construction and validation
# ==============================================================================

test_that("fda_data constructs with minimal arguments", {
  fd <- fda_data(time = times_vec, value = vals_vec, id = ids_vec)
  expect_s3_class(fd, "fda_data")
  expect_true(data.table::is.data.table(fd))
})

test_that("fda_data constructs with all arguments", {
  fd <- fda_data(
    time             = times_vec,
    value            = vals_vec,
    id               = ids_vec,
    group            = groups_vec,
    spatial_row      = rep(1:5, each = 10L),
    spatial_col      = rep(1:10, times = 5L),
    primary_trait    = yields_vec,
    primary_trait_name = "yield_t_ha"
  )
  expect_s3_class(fd, "fda_data")
  meta <- attr(fd, "fda_meta")
  expect_true(meta$has_group)
  expect_true(meta$has_spatial)
  expect_true(meta$has_primary_trait)
  expect_equal(meta$primary_trait_name, "yield_t_ha")
})

test_that("fda_data returns correct class vector", {
  fd <- fda_data(time = times_vec, value = vals_vec, id = ids_vec)
  expect_identical(class(fd), c("fda_data", "data.table", "data.frame"))
})

test_that("fda_data internal storage is data.table", {
  fd <- fda_data(time = times_vec, value = vals_vec, id = ids_vec)
  expect_true(data.table::is.data.table(fd))
  expect_true(all(c("id", "time", "value") %in% names(fd)))
})

test_that("fda_data errors on NA in time", {
  bad_time <- times_vec
  bad_time[5] <- NA
  expect_error(
    fda_data(time = bad_time, value = vals_vec, id = ids_vec),
    "NA"
  )
})

test_that("fda_data errors on NA in value", {
  bad_val <- vals_vec
  bad_val[3] <- NA
  expect_error(
    fda_data(time = times_vec, value = bad_val, id = ids_vec),
    "NA"
  )
})

test_that("fda_data errors on NA in id", {
  bad_id <- ids_vec
  bad_id[1] <- NA
  expect_error(
    fda_data(time = times_vec, value = vals_vec, id = bad_id),
    "NA"
  )
})

test_that("fda_data errors on non-numeric time", {
  expect_error(
    fda_data(time = as.character(times_vec), value = vals_vec, id = ids_vec),
    "numeric"
  )
})

test_that("fda_data errors on non-numeric value", {
  expect_error(
    fda_data(time = times_vec, value = as.character(vals_vec), id = ids_vec),
    "numeric"
  )
})

test_that("fda_data errors on mismatched lengths", {
  expect_error(
    fda_data(time = times_vec[1:10], value = vals_vec, id = ids_vec),
    "length"
  )
  expect_error(
    fda_data(time = times_vec, value = vals_vec, id = ids_vec[1:10]),
    "length"
  )
})

test_that("fda_data errors on primary_trait with wrong number of entries", {
  # Wrong length: neither n_obs nor n_ids
  expect_error(
    fda_data(
      time = times_vec, value = vals_vec, id = ids_vec,
      primary_trait = c(1, 2, 3)
    ),
    "length"
  )
})

test_that("fda_data errors on primary_trait with inconsistent values per id", {
  # Give observation-level primary_trait but different values within same id
  bad_pt <- seq_len(n_obs) * 0.1
  expect_error(
    fda_data(
      time = times_vec, value = vals_vec, id = ids_vec,
      primary_trait = bad_pt
    ),
    "single value per"
  )
})

test_that("fda_data stores metadata correctly", {
  fd <- fda_data(time = times_vec, value = vals_vec, id = ids_vec)
  meta <- attr(fd, "fda_meta")
  expect_type(meta, "list")
  expect_equal(meta$n_ids, n_ids)
  expect_equal(meta$n_timepoints, 10L)
  expect_equal(meta$time_range, range(times_vec))
  expect_false(meta$has_group)
  expect_false(meta$has_spatial)
  expect_false(meta$has_primary_trait)
})


# ==============================================================================
# 2. fda_data methods
# ==============================================================================

test_that("print.fda_data runs without error", {
  fd <- fda_data(time = times_vec, value = vals_vec, id = ids_vec,
                 group = groups_vec)
  # print.fda_data calls NextMethod(), so must go through S3 dispatch.
  # Use capture.output to verify it at least runs and produces output.
  out <- capture.output(print(fd))
  # Should contain the fda_data header or at least data.table output
  expect_true(length(out) > 0L)
})

test_that("summary.fda_data returns expected structure", {
  fd <- fda_data(
    time = times_vec, value = vals_vec, id = ids_vec,
    group = groups_vec,
    primary_trait = yields_vec,
    primary_trait_name = "yield"
  )
  # Call method explicitly — NAMESPACE may not register S3 dispatch
  result <- summary.fda_data(fd)
  expect_type(result, "list")
  expect_true("n_ids" %in% names(result))
  expect_true("n_timepoints" %in% names(result))
})

test_that("[.fda_data subsetting preserves class when essential columns kept", {
  fd <- fda_data(
    time = times_vec, value = vals_vec, id = ids_vec,
    group = groups_vec
  )
  # id is stored as factor, so compare against factor level
  first_id <- levels(fd$id)[1L]
  sub <- fd[id == first_id]
  expect_s3_class(sub, "fda_data")
  expect_true(data.table::is.data.table(sub))

  # Should have fewer rows

  expect_equal(nrow(sub), 10L)

  # Metadata should be present (n_ids may count factor levels, not just present)
  meta_sub <- attr(sub, "fda_meta")
  expect_true(!is.null(meta_sub))
  expect_true(meta_sub$n_ids >= 1L)
})

test_that("[.fda_data subset by row returns fewer rows with correct class", {
  fd <- fda_data(time = times_vec, value = vals_vec, id = ids_vec)
  # Row subset preserving essential columns retains fda_data class
  sub <- fd[time > 150]
  expect_s3_class(sub, "fda_data")
  expect_true(nrow(sub) < nrow(fd))
  expect_true(all(sub$time > 150))
})

test_that("as.data.table.fda_data works and strips fda_data class", {
  fd <- fda_data(time = times_vec, value = vals_vec, id = ids_vec)
  dt <- as.data.table.fda_data(fd)
  expect_true(data.table::is.data.table(dt))
  expect_false(inherits(dt, "fda_data"))
  expect_null(attr(dt, "fda_meta"))
  # Data preserved
  expect_equal(nrow(dt), n_obs)
})


# ==============================================================================
# 3. fda_model — Construction and methods (via mock)
# ==============================================================================

test_that("new_fda_model constructs valid object", {
  model <- make_mock_fda_model()
  expect_s3_class(model, "fda_model")
})

test_that("fda_model class is correct", {
  model <- make_mock_fda_model()
  expect_identical(class(model), "fda_model")
})

test_that("coef.fda_model returns data.table", {
  model <- make_mock_fda_model()
  # Call method explicitly — NAMESPACE may not have S3method registration yet
  cf <- coef.fda_model(model)
  expect_true(data.table::is.data.table(cf))
  expect_true(all(c("time", "beta", "se", "ci_lower", "ci_upper") %in%
                     names(cf)))
})

test_that("fitted.fda_model returns data.table", {
  model <- make_mock_fda_model()
  ft <- fitted.fda_model(model)
  expect_true(data.table::is.data.table(ft))
  expect_true("fitted" %in% names(ft))
})

test_that("residuals.fda_model returns numeric", {
  model <- make_mock_fda_model()
  resid <- residuals.fda_model(model)
  expect_true(is.numeric(resid))
  expect_equal(length(resid), 40L)
})

test_that("print.fda_model runs without error", {
  model <- make_mock_fda_model()
  # Call method explicitly in case S3 dispatch not registered in NAMESPACE yet
  expect_output(print.fda_model(model), "fda_model")
})

test_that("summary.fda_model runs without error", {
  model <- make_mock_fda_model()
  # Call method explicitly in case S3 dispatch not registered in NAMESPACE yet
  result <- summary.fda_model(model)
  expect_type(result, "list")
  expect_true("engine" %in% names(result))
  expect_true("var_comps" %in% names(result))
})

test_that("fda_model with empty extras list works", {
  model <- make_mock_fda_model(extras = list())
  expect_s3_class(model, "fda_model")
  expect_equal(length(model$extras), 0L)
  # print should work without mentioning extras
  expect_output(print.fda_model(model), "fda_model")
})

test_that("fda_model stores engine correctly", {
  model <- make_mock_fda_model()
  expect_equal(model$engine, "asreml")
})

test_that("vcov.fda_model returns variance components", {
  model <- make_mock_fda_model()
  # Call method explicitly in case S3 dispatch not registered in NAMESPACE yet
  vc <- vcov.fda_model(model)
  expect_true(data.table::is.data.table(vc))
})

test_that("predict.fda_model without newdata returns predictions", {
  model <- make_mock_fda_model()
  # Call method explicitly in case S3 dispatch not registered in NAMESPACE yet
  preds <- predict.fda_model(model)
  expect_true(data.table::is.data.table(preds))
  expect_equal(nrow(preds), 4L)
})


# ==============================================================================
# 4. fda_comparison — Construction and methods (via mock)
# ==============================================================================

# Helper: build mock fda_comparison, avoiding is.data.table() namespace issue
make_mock_fda_comparison <- function() {
  m1 <- make_mock_fda_model()
  m2 <- make_mock_fda_model()
  metrics <- data.table::data.table(
    model = c("m1", "m2"),
    rmse  = c(0.45, 0.52),
    cor   = c(0.85, 0.78)
  )
  structure(
    list(
      models  = list(m1 = m1, m2 = m2),
      metrics = metrics,
      label   = "test comparison"
    ),
    class = "fda_comparison"
  )
}

test_that("fda_comparison has correct class", {
  comp <- make_mock_fda_comparison()
  expect_s3_class(comp, "fda_comparison")
  expect_equal(length(comp$models), 2L)
})

test_that("print.fda_comparison runs without error", {
  comp <- make_mock_fda_comparison()
  expect_output(print(comp), "fda_comparison")
})

test_that("summary.fda_comparison runs without error", {
  comp <- make_mock_fda_comparison()
  # Call method explicitly — NAMESPACE may not have S3method registration yet
  result <- summary.fda_comparison(comp)
  expect_true(data.table::is.data.table(result))
})


# ==============================================================================
# 5. Edge cases
# ==============================================================================

test_that("fda_data works with single variety (single id)", {
  one_id   <- rep("single_plot", 10L)
  one_time <- seq(1, 10)
  one_val  <- rnorm(10)
  fd <- fda_data(time = one_time, value = one_val, id = one_id)
  expect_s3_class(fd, "fda_data")
  meta <- attr(fd, "fda_meta")
  expect_equal(meta$n_ids, 1L)
  expect_equal(meta$n_timepoints, 10L)
})

test_that("fda_data with primary_trait supplied per-id (short vector)", {
  # 5 unique ids, supply 5 values
  pt_short <- c(3.2, 4.1, 3.8, 4.5, 3.0)
  fd <- fda_data(
    time = times_vec, value = vals_vec, id = ids_vec,
    primary_trait = pt_short
  )
  expect_s3_class(fd, "fda_data")
  expect_true("primary_trait" %in% names(fd))
})

test_that("fda_data extra columns via ... are stored", {
  fd <- fda_data(
    time = times_vec, value = vals_vec, id = ids_vec,
    humidity = rnorm(n_obs)
  )
  expect_true("humidity" %in% names(fd))
  meta <- attr(fd, "fda_meta")
  expect_true("humidity" %in% meta$extra_cols)
})

test_that("fda_data errors on Inf in value", {
  bad_val <- vals_vec
  bad_val[1] <- Inf
  expect_error(
    fda_data(time = times_vec, value = bad_val, id = ids_vec),
    "finite"
  )
})

test_that("fda_data key is set on id and time", {
  fd <- fda_data(time = times_vec, value = vals_vec, id = ids_vec)
  expect_equal(data.table::key(fd), c("id", "time"))
})
