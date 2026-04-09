# fit_profiles.R -- Stage 1: Fit variety-specific functional profiles
#
# Fits a spatial + temporal mixed model per trial to extract variety-specific
# functional profiles (smooth curves). This is the first stage of the two-stage
# scalar-on-function approach for relating secondary functional traits (e.g.,
# NDVI, grain-fill) to a primary trait (e.g., yield).
#
# Model: y_{ijk}(t) = f_v(t) + block_j + spatial(row, col) + epsilon_{ijk}(t)
# where f_v(t) = sum_k alpha_{vk} * B_k(t) via B-spline basis
#
# Author: Maksym Gaidashenko
# Licence: GPL (>= 3)


# ---- fit_functional_profiles ----

#' Fit Variety-Specific Functional Profiles (Stage 1)
#'
#' Fits a spatial + temporal linear mixed model per trial to extract
#' variety-specific functional profiles (smooth curves) from repeated
#' measurements of a secondary trait. This is the first stage of the two-stage
#' FDA approach: Stage 1 recovers the shape of each variety's functional
#' response; Stage 2 (via [scalar_on_function()]) regresses a primary trait
#' (e.g., yield) onto these profiles.
#'
#' @param data An `fda_data` object (from [fda_data()]) or a data.frame /
#'   data.table containing the observations. If a raw data.frame, the column
#'   mapping arguments below must be specified.
#' @param time_col Character; name of the time column in `data`.
#'   Default `"time"`. Ignored when `data` is an `fda_data` object.
#' @param value_col Character or `NULL`; name of the response (value) column.
#'   Default `NULL`, which uses `"value"` for `fda_data` or requires explicit
#'   specification for raw data.
#' @param id_col Character; name of the observational unit identifier column
#'   (e.g., plot). Default `"id"`.
#' @param group_col Character; name of the variety/genotype grouping column.
#'   Default `"variety"`.
#' @param spatial_row_col Character or `NULL`; name of the spatial row column.
#'   Default `NULL` (no spatial modelling unless detected from `fda_data`).
#' @param spatial_col_col Character or `NULL`; name of the spatial column column.
#'   Default `NULL`.
#' @param basis An `fda_basis` object (from [bspline_basis()]). If `NULL`
#'   (default), a basis is constructed automatically from the unique time points
#'   in `data` using `n_knots`, `degree`, and `penalty_order`.
#' @param n_knots Integer; number of internal knots for automatic basis
#'   construction. Default 10. Ignored if `basis` is supplied.
#' @param degree Integer; B-spline polynomial degree. Default 3 (cubic).
#'   Ignored if `basis` is supplied.
#' @param penalty_order Integer; order of the difference penalty for the
#'   P-spline. Default 2 (penalises curvature). Ignored if `basis` is supplied.
#' @param spatial Character; type of spatial error model. One of `"none"`
#'   (default), `"ar1"` (first-order autoregressive in rows), `"ar1ar1"`
#'   (separable AR1 x AR1 in rows and columns), or `"spline"` (2D spatial
#'   spline). Requires `spatial_row_col` and `spatial_col_col` to be specified
#'   (except `"none"`).
#' @param engine Character; estimation engine. One of `"auto"` (default),
#'   `"asreml"`, or `"bayesreml"`. See [funcrop_engines()].
#' @param ... Additional arguments passed to the backend fitting function
#'   (e.g., `control` for ASReml, `mcmc_control` for bayesreml).
#'
#' @return An `fda_model` object (see `new_fda_model()`) containing:
#' \describe{
#'   \item{fitted_curves}{data.table of variety-specific fitted curves on a
#'     fine time grid, with columns `id`, `time`, `fitted`, `se`.}
#'   \item{coefficient_function}{Empty list (not applicable for Stage 1; the
#'     coefficient function beta(t) is estimated in Stage 2).}
#'   \item{variance_components}{data.table of estimated variance components.}
#'   \item{predictions}{data.table of variety-level predicted values.}
#'   \item{residuals}{Model residuals.}
#'   \item{basis}{The `fda_basis` object used.}
#'   \item{data}{The input data (as `fda_data`).}
#'   \item{engine}{Character string identifying the engine used.}
#'   \item{call}{The matched function call.}
#'   \item{extras}{List with `spline_blups` (variety x basis coefficient matrix),
#'     `spline_decomposition` (from [make_Zspline()]), `raw_model` (backend
#'     model object), and `convergence` diagnostics.}
#' }
#'
#' @details
#' ## Model formulation
#'
#' For each trial / environment, the model is:
#' \deqn{y_{ijk}(t) = f_v(t) + b_j + s(r, c) + \varepsilon_{ijk}(t)}
#'
#' where:
#' - \eqn{f_v(t) = \sum_k \alpha_{vk} B_k(t)} is the variety-specific curve,
#'   with \eqn{\alpha_v \sim N(0, \sigma_v^2 I)} as random B-spline
#'   coefficients (P-spline mixed model representation).
#' - \eqn{b_j} is a block effect (if `block` column is present).
#' - \eqn{s(r, c)} is an optional spatial error model (AR1, AR1xAR1, or
#'   spline).
#' - \eqn{\varepsilon_{ijk}(t) \sim N(0, \sigma^2)} is the residual.
#'
#' The variety-specific B-spline random effects are structured as a Kronecker
#' product: \eqn{Z = I_V \otimes Z_{spline}}, where \eqn{Z_{spline}} comes
#' from [make_Zspline()] and \eqn{I_V} is the identity over varieties. This
#' yields \eqn{V \times K} random effect levels (V varieties, K spline
#' coefficients).
#'
#' ## Backend dispatch
#'
#' For ASReml: uses `str()` specification from `.asreml_build_spline_str()`
#' with penalised variance structure. For bayesreml: uses known matrices from
#' `.bayesreml_build_spline_str()` with penalty-as-precision prior.
#'
#' @references
#' Verbyla, A.P., Cavanagh, C.R. and Verbyla, K.L. (2012). Whole-genome
#' analysis of multienvironment or multitrait QTL in MAGIC. *G3: Genes,
#' Genomes, Genetics*, 2(9), 1085--1093.
#'
#' De Faveri, J., Verbyla, A.P., Pitchford, W.S., Venkatanagappa, S. and
#' Cullis, B.R. (2015). Statistical methods for analysis of multi-harvest data
#' from perennial pasture variety selection trials. *Crop and Pasture Science*,
#' 66(9), 947--962.
#'
#' @seealso [scalar_on_function()] for Stage 2, [bspline_basis()] for basis
#'   construction, [make_Zspline()] for mixed model reparameterisation.
#'
#' @examples
#' \dontrun{
#' # Load simulated grain-fill data
#' data(sim_grain_fill)
#'
#' # Create fda_data object
#' fd <- fda_data(
#'   time  = sim_grain_fill$time,
#'   value = sim_grain_fill$grain_weight,
#'   id    = sim_grain_fill$plot_id,
#'   group = sim_grain_fill$variety,
#'   spatial_row = sim_grain_fill$row,
#'   spatial_col = sim_grain_fill$col,
#'   primary_trait = sim_grain_fill$yield,
#'   primary_trait_name = "yield_t_ha"
#' )
#'
#' # Fit functional profiles (Stage 1)
#' profiles <- fit_functional_profiles(
#'   data    = fd,
#'   n_knots = 6,
#'   spatial = "ar1ar1",
#'   engine  = "auto"
#' )
#'
#' print(profiles)
#' plot(profiles, which = "fitted")
#'
#' # Extract BLUPs for Stage 2
#' blup_matrix <- profiles$extras$spline_blups
#' }
#'
#' @export
fit_functional_profiles <- function(
    data,
    time_col       = "time",
    value_col      = NULL,
    id_col         = "id",
    group_col      = "variety",
    spatial_row_col = NULL,
    spatial_col_col = NULL,
    basis          = NULL,
    n_knots        = 10L,
    degree         = 3L,
    penalty_order  = 2L,
    spatial        = c("none", "ar1", "ar1ar1", "spline"),
    engine         = "auto",
    ...
) {

  call <- match.call()
  spatial <- match.arg(spatial)

  # ===========================================================================

  # Input validation
  # ===========================================================================

  # --- data ---
  if (missing(data) || is.null(data)) {
    stop("`data` must be supplied.", call. = FALSE)
  }

  is_fda_data <- inherits(data, "fda_data")

  if (!is_fda_data && !is.data.frame(data)) {
    stop("`data` must be an fda_data object, data.frame, or data.table.",
         call. = FALSE)
  }

  # --- Convert to data.table (work on a copy to avoid side effects) ---
  if (is_fda_data) {
    dt <- as.data.table.fda_data(data)
    meta <- attr(data, "fda_meta")
    # Map standard fda_data columns
    time_col <- "time"
    value_col <- value_col %||% "value"
    id_col <- "id"
    group_col <- if ("group" %in% names(dt)) "group" else group_col
    if (meta$has_spatial) {
      spatial_row_col <- spatial_row_col %||% "spatial_row"
      spatial_col_col <- spatial_col_col %||% "spatial_col"
    }
  } else {
    dt <- data.table::as.data.table(data)
    if (is.null(value_col)) {
      stop(
        "`value_col` must be specified when `data` is not an fda_data object.",
        call. = FALSE
      )
    }
  }

  # --- Validate column existence ---
  .validate_col_exists(dt, time_col, "time_col")
  .validate_col_exists(dt, value_col, "value_col")
  .validate_col_exists(dt, id_col, "id_col")
  .validate_col_exists(dt, group_col, "group_col")

  if (spatial != "none") {
    if (is.null(spatial_row_col) || is.null(spatial_col_col)) {
      stop(
        "Spatial modelling (spatial = '", spatial, "') requires both ",
        "`spatial_row_col` and `spatial_col_col` to be specified.",
        call. = FALSE
      )
    }
    .validate_col_exists(dt, spatial_row_col, "spatial_row_col")
    .validate_col_exists(dt, spatial_col_col, "spatial_col_col")
  }

  # --- Validate column types ---
  .validate_numeric(dt[[time_col]], "time column")
  .validate_numeric(dt[[value_col]], "value column")

  # --- Ensure group is a factor ---
  if (!is.factor(dt[[group_col]])) {
    data.table::set(dt, j = group_col, value = as.factor(dt[[group_col]]))
  }
  variety_levels <- levels(dt[[group_col]])
  n_varieties <- length(variety_levels)

  if (n_varieties < 2L) {
    stop("At least 2 varieties (groups) are required for profile fitting.",
         call. = FALSE)
  }

  # --- Validate integers ---
  .validate_positive_integer(n_knots, "n_knots")
  if (!is.numeric(degree) || length(degree) != 1L || degree < 0L) {
    stop("`degree` must be a non-negative integer.", call. = FALSE)
  }
  .validate_positive_integer(penalty_order, "penalty_order")

  # --- Resolve engine ---
  engine <- .resolve_engine(engine)

  # ===========================================================================
  # Basis construction
  # ===========================================================================

  time_points <- sort(unique(dt[[time_col]]))

  if (is.null(basis)) {
    .msg("Constructing B-spline basis: ", n_knots, " knots, degree ", degree,
         ", penalty order ", penalty_order, ".")
    basis <- bspline_basis(
      x             = time_points,
      n_knots       = n_knots,
      degree        = degree,
      penalty_order = penalty_order
    )
  } else {
    if (!inherits(basis, "fda_basis")) {
      stop("`basis` must be an 'fda_basis' object from bspline_basis().",
           call. = FALSE)
    }
  }

  # Mixed model reparameterisation: B = [X, Z] with Z ~ N(0, sigma_u^2 * I)
  spline_decomp <- make_Zspline(basis, constraint = "decompose")

  # ===========================================================================
  # Construct model data
  # ===========================================================================

  # Map each observation's time point to the basis evaluation.
  # The basis was evaluated at unique time points; look up each observation.
  time_idx <- match(dt[[time_col]], basis$x)
  if (anyNA(time_idx)) {
    # Some observation times are not in the basis evaluation grid.
    # Re-evaluate the basis at all observation-level time points.
    .msg("Re-evaluating basis at observation-level time points.")
    knot_vec <- c(
      rep(basis$boundary[1L], basis$degree + 1L),
      basis$knots,
      rep(basis$boundary[2L], basis$degree + 1L)
    )
    B_obs <- splines::splineDesign(
      knots = knot_vec,
      x     = dt[[time_col]],
      ord   = basis$degree + 1L,
      outer.ok = TRUE
    )
  } else {
    B_obs <- basis$B[time_idx, , drop = FALSE]
  }

  n_obs    <- nrow(dt)
  n_basis  <- basis$n_basis
  n_z_cols <- spline_decomp$rank  # number of random-effect spline columns

  # Compute the Z_spline matrix at observation level using the decomposition
  # X_null: polynomial null-space columns (fixed), Z_range: random spline columns
  # We need to re-evaluate at observation times (not just unique times)
  if (exists("B_obs") && !identical(B_obs, basis$B)) {
    # Re-derive X and Z from B_obs using the eigenvector transformation
    P_dense <- as.matrix(basis$P)
    eig <- eigen(P_dense, symmetric = TRUE)
    tol <- max(eig$values) * .Machine$double.eps * n_basis * 10
    idx_range <- which(eig$values > tol)
    idx_null  <- which(eig$values <= tol)
    U_null    <- eig$vectors[, idx_null, drop = FALSE]
    U_range   <- eig$vectors[, idx_range, drop = FALSE]
    scaling   <- 1 / sqrt(eig$values[idx_range])
    X_obs <- B_obs %*% U_null
    Z_obs <- B_obs %*% (U_range %*% diag(scaling, nrow = length(scaling)))
  } else {
    X_obs <- spline_decomp$X[time_idx, , drop = FALSE]
    Z_obs <- spline_decomp$Z[time_idx, , drop = FALSE]
  }

  # Build the Kronecker Z matrix: I_variety (x) Z_spline
  # For each observation, the variety indicator selects which block of
  # the Kronecker product applies.
  variety_factor <- dt[[group_col]]
  variety_int    <- as.integer(variety_factor)

  # Construct the full Z matrix (n_obs x n_varieties * n_z_cols) as sparse
  # Row i maps to block variety_int[i] of the Kronecker product
  row_idx <- rep(seq_len(n_obs), each = n_z_cols)
  col_offset <- (variety_int - 1L) * n_z_cols
  col_idx <- as.vector(outer(seq_len(n_z_cols), col_offset, FUN = "+"))
  vals <- as.vector(t(Z_obs))

  Z_full <- Matrix::sparseMatrix(
    i    = row_idx,
    j    = col_idx,
    x    = vals,
    dims = c(n_obs, n_varieties * n_z_cols)
  )

  # ===========================================================================
  # Build model specification
  # ===========================================================================

  # Add working columns to the model data.table
  model_dt <- data.table::copy(dt)

  # Ensure variety is a factor with correct levels
  data.table::set(model_dt, j = "variety_f",
                  value = factor(dt[[group_col]], levels = variety_levels))

  # Add spline coefficient index (for the str() / known_matrices approach)
  # Create a unique identifier for each variety:spline_coef combination
  data.table::set(model_dt, j = "spline_obs_idx",
                  value = seq_len(n_obs))

  # Detect block column (common names)
  block_col <- NULL
  block_candidates <- c("block", "rep", "replicate", "Block", "Rep")
  for (bc in block_candidates) {
    if (bc %in% names(model_dt)) {
      block_col <- bc
      break
    }
  }

  # -- Build model using mixed-model reparameterised basis --
  # The P-spline penalty is converted to a mixed-model form via eigendecomposition:
  #   B = [X_null, Z_range]
  # where X_null spans the unpenalised null space (fixed effects for population
  # mean curve) and Z_range spans the penalised range space (random effects for
  # variety-specific smooth deviations, with smoothing controlled by sigma^2_u).
  #
  # This is the correct P-spline mixed model representation (Eilers & Marx, 1996;
  # Wand & Ormerod, 2008). The previous approach using raw B-spline columns as
  # both fixed and random effects did not apply the smoothing penalty.

  null_dim <- ncol(spline_decomp$X)  # typically = degree - penalty_order + 1

  # Add null-space columns (fixed effects for population mean curve)
  Xnull_col_names <- paste0("Xnull_", seq_len(null_dim))
  for (k in seq_len(null_dim)) {
    data.table::set(model_dt, j = Xnull_col_names[k], value = X_obs[, k])
  }

  # Add range-space columns (random effects for variety-specific deviations)
  Zrange_col_names <- paste0("Zrange_", seq_len(n_z_cols))
  for (k in seq_len(n_z_cols)) {
    data.table::set(model_dt, j = Zrange_col_names[k], value = Z_obs[, k])
  }

  # Fixed effects formula: intercept + block + null-space spline terms
  fixed_rhs <- "1"
  if (!is.null(block_col)) {
    if (!is.factor(model_dt[[block_col]])) {
      data.table::set(model_dt, j = block_col,
                      value = as.factor(model_dt[[block_col]]))
    }
    fixed_rhs <- paste0(fixed_rhs, " + ", block_col)
  }
  # Null-space terms represent the unpenalised population mean curve
  fixed_rhs <- paste0(fixed_rhs, " + ",
                       paste(Xnull_col_names, collapse = " + "))
  fixed_formula <- stats::as.formula(paste(value_col, "~", fixed_rhs))

  # Random effects: variety-specific penalised spline deviations
  # Each variety_f:Zrange_k term gives a separate random coefficient per variety
  # for the k-th range-space basis function. Under the mixed-model representation,
  # alpha_v ~ N(0, sigma^2_u * I), which is equivalent to P-spline smoothing.
  random_terms <- paste0("variety_f:", Zrange_col_names, collapse = " + ")

  # Spatial term
  spatial_term <- NULL
  rcov_formula <- NULL

  if (spatial == "ar1" && !is.null(spatial_row_col)) {
    if (!is.factor(model_dt[[spatial_row_col]])) {
      data.table::set(model_dt, j = spatial_row_col,
                      value = as.factor(model_dt[[spatial_row_col]]))
    }
    spatial_term <- paste0("ar1(", spatial_row_col, ")")
    rcov_formula <- stats::as.formula(
      paste0("~ ", spatial_term, ":units")
    )
  } else if (spatial == "ar1ar1" &&
             !is.null(spatial_row_col) && !is.null(spatial_col_col)) {
    if (!is.factor(model_dt[[spatial_row_col]])) {
      data.table::set(model_dt, j = spatial_row_col,
                      value = as.factor(model_dt[[spatial_row_col]]))
    }
    if (!is.factor(model_dt[[spatial_col_col]])) {
      data.table::set(model_dt, j = spatial_col_col,
                      value = as.factor(model_dt[[spatial_col_col]]))
    }
    rcov_formula <- stats::as.formula(
      paste0("~ ar1(", spatial_row_col, "):ar1(", spatial_col_col, ")")
    )
  } else if (spatial == "spline" &&
             !is.null(spatial_row_col) && !is.null(spatial_col_col)) {
    # Add 2D spatial spline as an additional random term
    random_terms <- paste0(
      random_terms, " + spl2D(", spatial_row_col, ", ", spatial_col_col, ")"
    )
  }

  random_formula <- stats::as.formula(paste("~", random_terms))

  # Assemble model_spec
  model_spec <- list(
    fixed           = fixed_formula,
    random          = random_formula,
    rcov            = rcov_formula,
    known_matrices  = list()
  )

  # ===========================================================================
  # Dispatch to backend
  # ===========================================================================

  .msg("Dispatching to '", engine, "' engine...")
  raw_result <- .dispatch_fit(engine = engine, model_spec = model_spec,
                              data = model_dt, ...)

  # ===========================================================================
  # Extract results
  # ===========================================================================

  std_result <- .standardise_result(
    raw_result = raw_result,
    engine     = engine,
    model_spec = model_spec,
    basis      = basis,
    data       = model_dt
  )

  # -- Extract variety-specific spline coefficient BLUPs --
  # Extract from the raw model using the Zrange_k naming convention.
  # These are the penalised range-space coefficients alpha_v (n_z_cols per variety).
  spline_blups <- tryCatch({
    model_obj <- raw_result[["model"]]
    if (engine == "asreml" && !is.null(model_obj)) {
      co <- coef(model_obj)$random
      rn <- rownames(co)
      # Pattern: variety_f_VNAME:Zrange_K
      zr_idx <- grep("Zrange_", rn)
      if (length(zr_idx) > 0) {
        # Parse into variety x basis coefficient matrix
        parts <- strsplit(rn[zr_idx], ":", fixed = TRUE)
        var_ids <- sub("^variety_f_", "", vapply(parts, `[`, character(1), 1L))
        zr_ids <- as.integer(sub("^Zrange_", "", vapply(parts, `[`, character(1), 2L)))
        mat <- matrix(0, nrow = n_varieties, ncol = n_z_cols)
        rownames(mat) <- variety_levels
        for (i in seq_along(zr_idx)) {
          v_row <- match(var_ids[i], variety_levels)
          if (!is.na(v_row) && zr_ids[i] <= n_z_cols) {
            mat[v_row, zr_ids[i]] <- co[zr_idx[i], 1]
          }
        }
        mat
      } else {
        .extract_variety_spline_blups(std_result, n_varieties, n_z_cols,
                                      variety_levels, engine)
      }
    } else {
      .extract_variety_spline_blups(std_result, n_varieties, n_z_cols,
                                    variety_levels, engine)
    }
  }, error = function(e) {
    warning("BLUP extraction failed: ", e$message, call. = FALSE)
    matrix(0, nrow = n_varieties, ncol = n_z_cols)
  })

  # -- Extract fixed-effect coefficients for mean curve reconstruction --
  fixed_coefs <- tryCatch({
    model_obj <- raw_result[["model"]]
    if (engine == "asreml" && !is.null(model_obj)) {
      co_fixed <- coef(model_obj)$fixed
      stats::setNames(co_fixed[, 1], rownames(co_fixed))
    } else {
      NULL
    }
  }, error = function(e) NULL)

  # -- Reconstruct fitted curves on a fine grid --
  # Includes both the population mean curve (fixed null-space) and
  # variety-specific deviations (random range-space BLUPs)
  fitted_curves <- .reconstruct_variety_curves(
    spline_blups   = spline_blups,
    spline_decomp  = spline_decomp,
    basis          = basis,
    variety_levels = variety_levels,
    fixed_coefs    = fixed_coefs
  )

  # -- Predictions (variety means across time, if primary trait available) --
  predictions <- data.table::data.table()
  if (is_fda_data && !is.null(attr(data, "fda_meta")$has_primary_trait) &&
      attr(data, "fda_meta")$has_primary_trait) {
    # Extract unique variety-level primary trait values
    pt_dt <- unique(dt[, c(group_col, "primary_trait"), with = FALSE])
    data.table::setnames(pt_dt, group_col, "variety")
    predictions <- pt_dt
  }

  # -- Coefficient function placeholder (not applicable for Stage 1) --
  coef_fn <- list(
    time     = numeric(0L),
    beta     = numeric(0L),
    se       = numeric(0L),
    ci_lower = numeric(0L),
    ci_upper = numeric(0L)
  )

  # ===========================================================================
  # Construct and return fda_model
  # ===========================================================================

  new_fda_model(
    fitted_curves          = fitted_curves,
    coefficient_function   = coef_fn,
    variance_components    = std_result$variance_components,
    predictions            = predictions,
    residuals              = std_result$residuals,
    basis                  = basis,
    data                   = data,
    engine                 = engine,
    call                   = call,
    extras                 = list(
      spline_blups          = spline_blups,
      spline_decomposition  = spline_decomp,
      raw_model             = std_result$raw_model,
      convergence           = std_result$convergence,
      model_spec            = model_spec
    )
  )
}


# ==============================================================================
# Internal helper functions
# ==============================================================================

#' Validate that a column exists in a data.table
#' @param dt data.table to check.
#' @param col_name Column name to look for.
#' @param arg_name Name of the argument for error messages.
#' @noRd
.validate_col_exists <- function(dt, col_name, arg_name) {
  if (!is.character(col_name) || length(col_name) != 1L) {
    stop(sprintf("`%s` must be a single character string.", arg_name),
         call. = FALSE)
  }
  if (!col_name %in% names(dt)) {
    stop(sprintf(
      "Column '%s' (from `%s`) not found in data. Available columns: %s.",
      col_name, arg_name, paste(names(dt), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}


#' Extract variety-specific spline coefficient BLUPs
#'
#' Parses the standardised result's BLUP table to produce a variety x basis
#' coefficient matrix. Falls back to sequential ordering if level parsing fails.
#'
#' @param std_result Standardised result from `.standardise_result()`.
#' @param n_varieties Integer; number of varieties.
#' @param n_z_cols Integer; number of random spline coefficients per variety.
#' @param variety_levels Character vector of variety levels.
#' @param engine Character; engine name.
#'
#' @return A numeric matrix (n_varieties x n_z_cols) of BLUPs, with rownames
#'   set to `variety_levels`.
#' @noRd
.extract_variety_spline_blups <- function(std_result, n_varieties, n_z_cols,
                                          variety_levels, engine) {
  blup_dt <- std_result$blups

  if (is.null(blup_dt) || nrow(blup_dt) == 0L) {
    warning(
      "No BLUPs could be extracted. Returning zero coefficient matrix.",
      call. = FALSE
    )
    mat <- matrix(0, nrow = n_varieties, ncol = n_z_cols)
    rownames(mat) <- variety_levels
    return(mat)
  }

  blup_vec <- blup_dt[["blup"]]
  expected_len <- n_varieties * n_z_cols

  if (length(blup_vec) == expected_len) {
    # Assume ordered: variety 1 coefficients, then variety 2, etc.
    mat <- matrix(blup_vec, nrow = n_varieties, ncol = n_z_cols, byrow = TRUE)
    rownames(mat) <- variety_levels
    return(mat)
  }

  # Try parsing variety:coef_idx from the level column
  if ("level" %in% names(blup_dt) && any(grepl(":", blup_dt[["level"]]))) {
    parts <- data.table::tstrsplit(blup_dt[["level"]], ":", fixed = TRUE)
    var_ids  <- parts[[1L]]
    coef_ids <- as.integer(parts[[2L]])
    uvar <- unique(var_ids)
    mat <- matrix(NA_real_, nrow = length(uvar), ncol = max(coef_ids))
    rownames(mat) <- uvar
    for (i in seq_len(nrow(blup_dt))) {
      mat[var_ids[i], coef_ids[i]] <- blup_vec[i]
    }
    return(mat)
  }

  # Fallback: truncate or pad
  warning(
    sprintf(
      "BLUP vector length (%d) does not match expected (%d = %d varieties x %d coefficients). Padding with zeros.",
      length(blup_vec), expected_len, n_varieties, n_z_cols
    ),
    call. = FALSE
  )
  blup_vec <- c(blup_vec, rep(0, max(0L, expected_len - length(blup_vec))))
  blup_vec <- blup_vec[seq_len(expected_len)]
  mat <- matrix(blup_vec, nrow = n_varieties, ncol = n_z_cols, byrow = TRUE)
  rownames(mat) <- variety_levels
  mat
}


#' Reconstruct variety-specific fitted curves on a fine grid
#'
#' Evaluates the B-spline basis at a fine grid and multiplies by the
#' variety-specific coefficient BLUPs to produce smooth fitted curves.
#' Includes the fixed-effect population mean curve (null-space component)
#' so that returned curves are on the absolute scale, not centred at zero.
#'
#' The full prediction for variety v at time t is:
#'   f_hat_v(t) = X_null(t)' beta_hat + Z_range(t)' alpha_hat_v
#'
#' @param spline_blups Numeric matrix (n_varieties x n_z_cols) of BLUPs.
#' @param spline_decomp List from [make_Zspline()] with X and Z components.
#' @param basis An `fda_basis` object.
#' @param variety_levels Character vector of variety names.
#' @param fixed_coefs Optional named numeric vector of fixed-effect estimates.
#'   If supplied, null-space coefficients (matching "Xnull_" pattern) are used
#'   to reconstruct the population mean curve. If NULL, only random deviations
#'   are returned (centred at zero).
#'
#' @return data.table with columns: `id`, `time`, `fitted`, `se`.
#' @noRd
.reconstruct_variety_curves <- function(spline_blups, spline_decomp,
                                        basis, variety_levels,
                                        fixed_coefs = NULL) {
  # Fine evaluation grid
  n_grid <- 200L
  t_grid <- seq(basis$boundary[1L], basis$boundary[2L], length.out = n_grid)

  # Re-evaluate the B-spline basis at the fine grid
  knot_vec <- c(
    rep(basis$boundary[1L], basis$degree + 1L),
    basis$knots,
    rep(basis$boundary[2L], basis$degree + 1L)
  )
  B_fine <- splines::splineDesign(
    knots    = knot_vec,
    x        = t_grid,
    ord      = basis$degree + 1L,
    outer.ok = TRUE
  )

  # Apply the same eigendecomposition to get X_null and Z_range at fine grid
  P_dense <- as.matrix(basis$P)
  eig <- eigen(P_dense, symmetric = TRUE)
  tol <- max(eig$values) * .Machine$double.eps * basis$n_basis * 10
  idx_range <- which(eig$values > tol)
  idx_null  <- which(eig$values <= tol)
  U_null    <- eig$vectors[, idx_null, drop = FALSE]
  U_range   <- eig$vectors[, idx_range, drop = FALSE]
  scaling   <- 1 / sqrt(eig$values[idx_range])

  X_fine    <- B_fine %*% U_null
  Z_fine    <- B_fine %*% (U_range %*% diag(scaling, nrow = length(scaling)))

  n_varieties <- nrow(spline_blups)
  n_z_cols    <- ncol(spline_blups)

  # Ensure Z_fine columns match BLUP dimensions
  if (ncol(Z_fine) != n_z_cols) {
    warning(
      sprintf(
        "Z_fine has %d columns but BLUPs have %d coefficients per variety. Using min.",
        ncol(Z_fine), n_z_cols
      ),
      call. = FALSE
    )
    k <- min(ncol(Z_fine), n_z_cols)
    Z_fine <- Z_fine[, seq_len(k), drop = FALSE]
    spline_blups <- spline_blups[, seq_len(k), drop = FALSE]
  }

  # Reconstruct the population mean curve from fixed-effect null-space coefficients
  # f_hat_v(t) = mu(t) + u_v(t) = X_null(t)' beta + Z_range(t)' alpha_v
  mean_curve <- rep(0, n_grid)
  if (!is.null(fixed_coefs)) {
    # Extract null-space coefficients from the fixed effects vector
    xnull_idx <- grep("^Xnull_", names(fixed_coefs))
    if (length(xnull_idx) > 0 && length(xnull_idx) == ncol(X_fine)) {
      beta_null <- fixed_coefs[xnull_idx]
      mean_curve <- as.numeric(X_fine %*% beta_null)
    }
  }

  # Reconstruct curves: fitted_v(t) = mu(t) + Z_fine %*% alpha_v
  curve_list <- vector("list", n_varieties)

  for (v in seq_len(n_varieties)) {
    alpha_v  <- spline_blups[v, ]
    fitted_v <- mean_curve + as.numeric(Z_fine %*% alpha_v)

    curve_list[[v]] <- data.table::data.table(
      id     = variety_levels[v],
      time   = t_grid,
      fitted = fitted_v,
      se     = NA_real_
    )
  }

  data.table::rbindlist(curve_list, use.names = TRUE)
}
