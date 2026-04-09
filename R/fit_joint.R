# fit_joint.R -- Single-stage joint model and 2D functional modelling
#
# Implements the single-stage joint model for simultaneously estimating a
# scalar primary trait and a functional secondary trait, avoiding the error
# propagation issues of the two-stage approach. Also provides 2D functional
# profiling via tensor product B-splines and scalar-on-2D-function regression.
#
# Author: Max Moldovan
# Licence: GPL (>= 3)


# ==============================================================================
# fit_fda_joint -- Single-stage joint model
# ==============================================================================

#' Fit a Joint Model for Primary and Functional Secondary Traits
#'
#' Fits a single-stage model that simultaneously models a scalar primary trait
#' (e.g., yield) and a functional secondary trait (e.g., grain-fill over time),
#' allowing the functional trait to inform estimation of the primary trait.
#' This avoids the error propagation issues inherent in the two-stage approach
#' where functional profiles are first estimated, then used as predictors.
#'
#' @param data A `data.frame` or [data.table::data.table] in long format
#'   containing both the functional secondary trait observations (one row per
#'   id x time combination) and the primary trait (constant within each id).
#'   May also be an `fda_data` object.
#' @param primary_col Character; column name for the scalar primary trait
#'   (e.g., yield). Default `"yield"`.
#' @param time_col Character; column name for the time/index variable of the
#'   functional trait. Default `"time"`.
#' @param secondary_col Character or `NULL`; column name for the functional
#'   secondary trait values. If `NULL` (default), inferred as the first numeric
#'   column that is not `primary_col`, `time_col`, or a spatial column.
#' @param id_col Character; column name identifying observational units
#'   (e.g., plot). Default `"id"`.
#' @param group_col Character; column name for the grouping factor
#'   (e.g., variety, genotype). Default `"variety"`.
#' @param spatial_row_col Character or `NULL`; column name for spatial row
#'   coordinates. Default `NULL` (no spatial modelling).
#' @param spatial_col_col Character or `NULL`; column name for spatial column
#'   coordinates. Default `NULL`.
#' @param basis An `fda_basis` object (from [bspline_basis()]). If `NULL`
#'   (default), a basis is constructed from the data using `n_knots`, `degree`,
#'   and `penalty_order`.
#' @param n_knots Integer; number of internal B-spline knots if `basis` is
#'   `NULL`. Default 10.
#' @param degree Integer; B-spline polynomial degree. Default 3 (cubic).
#' @param penalty_order Integer; order of the difference penalty. Default 2.
#' @param spatial Character; spatial correlation structure for residuals. One
#'   of `"none"` (default), `"ar1"` (AR1 in rows), or `"ar1ar1"` (separable
#'   AR1 x AR1 in rows and columns).
#' @param genomic_matrix Optional numeric matrix (square, symmetric) of genomic
#'   relationships among groups. Row/column names must match levels of
#'   `group_col`. Default `NULL`.
#' @param pedigree_matrix Optional numeric matrix of pedigree-based
#'   relationships. Same naming requirements as `genomic_matrix`.
#'   Default `NULL`.
#' @param engine Character; estimation engine. One of `"auto"` (default),
#'   `"asreml"`, or `"bayesreml"`. See [funcrop_engines()].
#' @param ... Additional arguments passed to the backend fitting function
#'   (e.g., `maxiter`, `workspace` for ASReml; `mcmc_control` for bayesreml).
#'
#' @return An `fda_model` object containing:
#' \describe{
#'   \item{fitted_curves}{[data.table::data.table] of variety-specific fitted
#'     functional curves for the secondary trait.}
#'   \item{coefficient_function}{List with the implied relationship
#'     \eqn{\beta(t)} between the secondary functional trait and the primary
#'     scalar trait.}
#'   \item{variance_components}{[data.table::data.table] of estimated variance
#'     components including the cross-trait correlations.}
#'   \item{predictions}{[data.table::data.table] of predicted primary trait
#'     values per group.}
#'   \item{basis}{The `fda_basis` object used for the functional trait.}
#'   \item{engine}{Character string identifying the backend used.}
#' }
#'
#' @section Statistical Model:
#' The joint model treats the secondary trait at each time point and the primary
#' trait as a multivariate response. In stacked form:
#'
#' \deqn{
#' \begin{bmatrix} y_{\text{secondary}}(t) \\ y_{\text{primary}} \end{bmatrix}
#' = \begin{bmatrix} f_v(t) + \text{spatial} + \varepsilon_s(t) \\
#'   \alpha + \int \beta(t) f_v(t)\, dt + \varepsilon_p \end{bmatrix}
#' }
#'
#' where \eqn{f_v(t) = \sum_k \alpha_{vk} B_k(t)} is the variety-specific
#' B-spline random effect for the secondary trait, and the primary trait is
#' linked through the functional integral \eqn{\int \beta(t) f_v(t)\, dt}.
#'
#' The model is fitted as a bivariate mixed model with:
#' - A trait indicator variable distinguishing secondary vs primary observations
#' - B-spline basis x variety interactions for the secondary trait
#' - Numerical integration weights for the functional covariate in the primary
#'   trait equation
#' - Correlated random effects across traits via `us(trait)` structures
#'
#' @section Implementation:
#' Data are stacked with a `trait_type` factor (`"secondary"`, `"primary"`).
#' For secondary-trait rows, the design matrix uses B-spline basis x variety
#' interactions. For primary-trait rows, the design matrix uses the numerical
#' integral of the B-spline basis functions (trapezoidal rule), creating a
#' functional covariate. Variety-level random effects are correlated across
#' both trait types.
#'
#' @references
#' De Boor, C. (2001). *A Practical Guide to Splines* (Revised ed.).
#' Springer.
#'
#' Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines
#' and penalties. *Statistical Science*, 11(2), 89--121.
#'
#' @seealso [fit_functional_profiles()] for Stage 1 of the two-stage approach,
#'   [scalar_on_function()] for Stage 2, [bspline_basis()] for basis
#'   construction.
#'
#' @examples
#' \dontrun{
#' data(sim_grain_fill)
#'
#' # Fit the joint model
#' joint_fit <- fit_fda_joint(
#'   data          = sim_grain_fill,
#'   primary_col   = "yield",
#'   time_col      = "time",
#'   secondary_col = "grain_weight",
#'   id_col        = "plot_id",
#'   group_col     = "variety",
#'   spatial_row_col = "row",
#'   spatial_col_col = "col",
#'   n_knots       = 6,
#'   spatial       = "ar1ar1",
#'   engine        = "asreml"
#' )
#' print(joint_fit)
#' summary(joint_fit)
#' }
#'
#' @export
fit_fda_joint <- function(
    data,
    primary_col     = "yield",
    time_col        = "time",
    secondary_col   = NULL,
    id_col          = "id",
    group_col       = "variety",
    spatial_row_col = NULL,
    spatial_col_col = NULL,
    basis           = NULL,
    n_knots         = 10L,
    degree          = 3L,
    penalty_order   = 2L,
    spatial         = c("none", "ar1", "ar1ar1"),
    genomic_matrix  = NULL,
    pedigree_matrix = NULL,
    engine          = "auto",
    ...) {

  mc <- match.call()
  spatial <- match.arg(spatial)

  # ---- Input validation -------------------------------------------------------

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or data.table.", call. = FALSE)
  }
  dt <- data.table::as.data.table(data)

  # Validate required columns
  required_cols <- c(primary_col, time_col, id_col, group_col)
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("Column(s) not found in data: %s.",
              paste0("'", missing_cols, "'", collapse = ", ")),
      call. = FALSE
    )
  }

  # Infer secondary_col if not specified
  if (is.null(secondary_col)) {
    exclude_cols <- c(primary_col, time_col, id_col, group_col,
                      spatial_row_col, spatial_col_col)
    numeric_cols <- names(dt)[vapply(dt, is.numeric, logical(1L))]
    candidates <- setdiff(numeric_cols, exclude_cols)
    if (length(candidates) == 0L) {
      stop(
        "Cannot infer `secondary_col`. No candidate numeric columns remain. ",
        "Please specify `secondary_col` explicitly.",
        call. = FALSE
      )
    }
    secondary_col <- candidates[1L]
    .msg(sprintf("Inferred secondary_col = '%s'.", secondary_col))
  }

  if (!secondary_col %in% names(dt)) {
    stop(sprintf("Column '%s' not found in data.", secondary_col),
         call. = FALSE)
  }

  # Validate spatial columns if specified
  if (!is.null(spatial_row_col) && !spatial_row_col %in% names(dt)) {
    stop(sprintf("spatial_row_col '%s' not found in data.", spatial_row_col),
         call. = FALSE)
  }
  if (!is.null(spatial_col_col) && !spatial_col_col %in% names(dt)) {
    stop(sprintf("spatial_col_col '%s' not found in data.", spatial_col_col),
         call. = FALSE)
  }

  if (spatial != "none" && (is.null(spatial_row_col) || is.null(spatial_col_col))) {
    stop(
      "Spatial modelling requires both `spatial_row_col` and `spatial_col_col`.",
      call. = FALSE
    )
  }

  .validate_numeric(dt[[time_col]], time_col)
  .validate_numeric(dt[[secondary_col]], secondary_col)
  .validate_numeric(dt[[primary_col]], primary_col)

  # Resolve engine

  engine <- .resolve_engine(engine)

  # ---- Construct B-spline basis -----------------------------------------------

  time_vals <- dt[[time_col]]

  if (is.null(basis)) {
    basis <- bspline_basis(
      x             = sort(unique(time_vals)),
      n_knots       = n_knots,
      degree        = degree,
      penalty_order = penalty_order
    )
    .msg(sprintf(
      "Constructed B-spline basis: %d knots, degree %d, %d basis functions.",
      length(basis$knots), basis$degree, basis$n_basis
    ))
  } else {
    if (!inherits(basis, "fda_basis")) {
      stop("`basis` must be an 'fda_basis' object.", call. = FALSE)
    }
  }

  # Mixed model reparameterisation
  mm_reparam <- make_Zspline(basis, constraint = "decompose")

  # ---- Evaluate basis at observed time points ---------------------------------

  # Create B-spline evaluation at all observed time points
  B_obs <- splines::splineDesign(
    knots = c(rep(basis$boundary[1L], basis$degree + 1L),
              basis$knots,
              rep(basis$boundary[2L], basis$degree + 1L)),
    x     = time_vals,
    ord   = basis$degree + 1L,
    outer.ok = TRUE
  )

  # ---- Compute integration weights for functional covariate -------------------
  # Trapezoidal rule on the unique time grid
  t_unique <- sort(unique(time_vals))
  n_t <- length(t_unique)

  # Trapezoidal weights
  trap_weights <- numeric(n_t)
  if (n_t >= 2L) {
    dt_diff <- diff(t_unique)
    trap_weights[1L] <- dt_diff[1L] / 2
    trap_weights[n_t] <- dt_diff[n_t - 1L] / 2
    if (n_t > 2L) {
      for (k in 2L:(n_t - 1L)) {
        trap_weights[k] <- (dt_diff[k - 1L] + dt_diff[k]) / 2
      }
    }
  } else {
    trap_weights <- 1
  }

  # Evaluate basis at unique time points for integration
  B_unique <- splines::splineDesign(
    knots = c(rep(basis$boundary[1L], basis$degree + 1L),
              basis$knots,
              rep(basis$boundary[2L], basis$degree + 1L)),
    x     = t_unique,
    ord   = basis$degree + 1L,
    outer.ok = TRUE
  )

  # Integrated basis: integral of each basis function over the time domain
  # int_basis[k] = sum_j w_j * B_k(t_j)
  int_basis <- as.numeric(crossprod(B_unique, trap_weights))  # length n_basis

  # ---- Prepare stacked data ---------------------------------------------------

  # Ensure factor columns
  dt[, (group_col) := as.factor(dt[[group_col]])]
  dt[, (id_col)    := as.factor(dt[[id_col]])]

  group_levels <- levels(dt[[group_col]])
  id_levels    <- levels(dt[[id_col]])
  n_groups     <- length(group_levels)
  n_ids        <- length(id_levels)
  n_basis      <- basis$n_basis

  # -- Secondary trait rows (one row per observation = id x time) --
  dt_secondary <- data.table::copy(dt)
  dt_secondary[, response := dt_secondary[[secondary_col]]]
  dt_secondary[, trait_type := factor("secondary",
                                       levels = c("secondary", "primary"))]
  dt_secondary[, row_id := .I]

  # Add B-spline columns for secondary trait
  bspline_col_names <- paste0("bs_", seq_len(n_basis))
  for (k in seq_len(n_basis)) {
    data.table::set(dt_secondary, j = bspline_col_names[k], value = B_obs[, k])
  }

  # -- Primary trait rows (one row per id) --
  dt_primary <- unique(dt[, c(id_col, group_col, primary_col,
                                spatial_row_col, spatial_col_col),
                            with = FALSE])

  dt_primary[, response := dt_primary[[primary_col]]]
  dt_primary[, trait_type := factor("primary",
                                     levels = c("secondary", "primary"))]
  dt_primary[, (time_col) := NA_real_]

  # For primary trait rows, set the B-spline columns to integrated basis values.
  # Each primary observation gets the same integrated basis vector -- the actual

  # variety-level integration is handled via the random effect structure.
  for (k in seq_len(n_basis)) {
    data.table::set(dt_primary, j = bspline_col_names[k], value = int_basis[k])
  }

  # Ensure matching columns before stacking
  common_cols <- c(id_col, group_col, time_col, "response", "trait_type",
                   bspline_col_names)
  if (!is.null(spatial_row_col)) common_cols <- c(common_cols, spatial_row_col)
  if (!is.null(spatial_col_col)) common_cols <- c(common_cols, spatial_col_col)

  # Add any missing columns to primary as NA
  for (cc in setdiff(common_cols, names(dt_primary))) {
    data.table::set(dt_primary, j = cc, value = NA)
  }
  for (cc in setdiff(common_cols, names(dt_secondary))) {
    data.table::set(dt_secondary, j = cc, value = NA)
  }

  stacked <- data.table::rbindlist(
    list(dt_secondary[, common_cols, with = FALSE],
         dt_primary[, common_cols, with = FALSE]),
    use.names = TRUE,
    fill = TRUE
  )

  # Coerce factors
  stacked[, trait_type := as.factor(trait_type)]
  stacked[, (group_col) := as.factor(stacked[[group_col]])]
  stacked[, (id_col)    := as.factor(stacked[[id_col]])]

  if (!is.null(spatial_row_col)) {
    stacked[, (spatial_row_col) := as.factor(stacked[[spatial_row_col]])]
  }
  if (!is.null(spatial_col_col)) {
    stacked[, (spatial_col_col) := as.factor(stacked[[spatial_col_col]])]
  }

  # ---- Build model specification ----------------------------------------------

  # Fixed effects: trait-specific intercept + population mean spline terms.
  # The z_spl columns capture the population-level functional relationship:
  # for secondary rows these are the spline basis values Z(t), for primary
  # rows these are the integrated basis values c_k. Including them as fixed
  # effects estimates the population mean curve + population-level coefficient
  # function.
  fixed_rhs <- paste(c("trait_type", z_col_names), collapse = " + ")
  fixed_formula <- stats::as.formula(paste("response ~", fixed_rhs))

  # Random effects: variety B-spline coefficients correlated across traits
  # For ASReml: us(trait_type):group + at(trait_type, 'secondary'):id:spline
  # For bayesreml: equivalent structure

  # Create the Z matrix for B-spline random effects
  Z_spline <- mm_reparam$Z
  X_spline <- mm_reparam$X

  # Add Z-spline columns to stacked data
  z_col_names <- paste0("z_spl_", seq_len(ncol(Z_spline)))

  # For secondary rows: evaluate Z at the observed time points
  Z_obs <- B_obs %*%
    (eigen(as.matrix(basis$P), symmetric = TRUE)$vectors[,
      seq_len(mm_reparam$rank), drop = FALSE] %*%
      diag(1 / sqrt(
        eigen(as.matrix(basis$P), symmetric = TRUE)$values[
          seq_len(mm_reparam$rank)]),
        nrow = mm_reparam$rank))

  n_secondary <- nrow(dt_secondary)
  n_primary   <- nrow(dt_primary)

  for (k in seq_len(ncol(Z_spline))) {
    # Secondary rows get Z evaluated at observed time
    z_sec <- if (k <= ncol(Z_obs)) Z_obs[, k] else 0
    # Primary rows get the integrated Z (trapezoidal integration of Z columns)
    Z_unique_k <- splines::splineDesign(
      knots = c(rep(basis$boundary[1L], basis$degree + 1L),
                basis$knots,
                rep(basis$boundary[2L], basis$degree + 1L)),
      x     = t_unique,
      ord   = basis$degree + 1L,
      outer.ok = TRUE
    )
    # Re-compute Z for unique time points through the same transformation
    Z_unique_reparam <- Z_unique_k %*%
      (eigen(as.matrix(basis$P), symmetric = TRUE)$vectors[,
        seq_len(mm_reparam$rank), drop = FALSE] %*%
        diag(1 / sqrt(
          eigen(as.matrix(basis$P), symmetric = TRUE)$values[
            seq_len(mm_reparam$rank)]),
          nrow = mm_reparam$rank))

    z_int <- sum(trap_weights * Z_unique_reparam[, k])
    z_pri <- rep(z_int, n_primary)

    data.table::set(stacked, j = z_col_names[k],
                    value = c(z_sec, z_pri))
  }

  # Build model_spec for dispatch
  # Random formula: variety-specific spline coefficients shared across traits.
  #
  # The key innovation of the joint model is that both traits are driven by
  # the same latent coefficient vector alpha_g:
  #   Secondary: s(t) = z(t)' alpha_g  (spline basis at observed time)
  #   Primary:   y    = c'   alpha_g  (integrated basis: c_k = int z_k(t) dt)
  #
  # The z_spl columns already encode this: secondary rows have Z(t) values,
  # primary rows have the integrated c_k values. So using variety:z_spl_k
  # as the random terms creates the shared latent structure automatically.
  #
  # Previous implementation (v0.1.0) used us(trait_type):variety which gave
  # only a trait-by-variety intercept — the spline link was not active.
  random_parts <- character(0L)

  # Shared variety-specific spline coefficients (the functional link)
  for (k in seq_len(ncol(Z_spline))) {
    random_parts <- c(random_parts,
                       sprintf("%s:%s", group_col, z_col_names[k]))
  }

  random_formula <- stats::as.formula(
    paste("~", paste(random_parts, collapse = " + "))
  )

  # Residual covariance
  if (spatial == "none") {
    rcov_formula <- stats::as.formula("~ units")
  } else if (spatial == "ar1") {
    rcov_formula <- stats::as.formula(
      sprintf("~ at(trait_type, 'secondary'):%s:ar1(%s) + at(trait_type, 'primary'):units",
              id_col, spatial_row_col)
    )
  } else {
    # ar1ar1
    rcov_formula <- stats::as.formula(
      sprintf(
        paste0("~ at(trait_type, 'secondary'):ar1(%s):ar1(%s) + ",
               "at(trait_type, 'primary'):units"),
        spatial_row_col, spatial_col_col
      )
    )
  }

  # Known matrices
  known_mats <- list()
  if (!is.null(genomic_matrix)) {
    .validate_relationship_matrix(genomic_matrix, group_levels, "genomic_matrix")
    known_mats[["G_genomic"]] <- genomic_matrix
  }
  if (!is.null(pedigree_matrix)) {
    .validate_relationship_matrix(pedigree_matrix, group_levels, "pedigree_matrix")
    known_mats[["A_pedigree"]] <- pedigree_matrix
  }

  model_spec <- list(
    fixed           = fixed_formula,
    random          = random_formula,
    rcov            = rcov_formula,
    known_matrices  = known_mats
  )

  # Add extra arguments
  dots <- list(...)
  if (length(dots) > 0L) {
    model_spec <- c(model_spec, dots)
  }

  # ---- Dispatch to backend ----------------------------------------------------

  .msg(sprintf(
    "Fitting joint model: %d secondary obs + %d primary obs, %d groups, engine = '%s'.",
    n_secondary, n_primary, n_groups, engine
  ))

  raw_result <- .dispatch_fit(
    engine     = engine,
    model_spec = model_spec,
    data       = stacked,
    ...
  )

  # ---- Standardise results ----------------------------------------------------

  std_result <- .standardise_result(
    raw_result = raw_result,
    engine     = engine,
    model_spec = model_spec,
    basis      = basis,
    data       = stacked
  )

  # Build predictions for primary trait
  predictions <- data.table::data.table(
    group     = group_levels,
    predicted = rep(NA_real_, n_groups),
    se        = rep(NA_real_, n_groups)
  )

  # Attempt to extract primary trait predictions from the raw model
  tryCatch({
    if (engine == "asreml") {
      blups <- .asreml_extract_blups(
        raw_result[["model"]],
        sprintf("trait_type:%s", group_col)
      )
      if (nrow(blups) > 0L) {
        primary_blups <- blups[grepl("primary", level)]
        if (nrow(primary_blups) > 0L) {
          predictions[, predicted := primary_blups$blup[seq_len(.N)]]
          predictions[, se := primary_blups$se[seq_len(.N)]]
        }
      }
    } else if (engine == "bayesreml") {
      blups <- .bayesreml_extract_blups(
        raw_result[["model"]],
        sprintf("trait_type:%s", group_col)
      )
      if (nrow(blups) > 0L) {
        primary_blups <- blups[grepl("primary", level)]
        if (nrow(primary_blups) > 0L) {
          predictions[, predicted := primary_blups$blup[seq_len(.N)]]
          predictions[, se := primary_blups$sd[seq_len(.N)]]
        }
      }
    }
  }, error = function(e) {
    .msg(sprintf("Could not extract primary trait predictions: %s",
                 conditionMessage(e)))
  })

  # ---- Construct fda_model object ---------------------------------------------

  model_out <- new_fda_model(
    fitted_curves        = std_result$fitted_curves,
    coefficient_function = std_result$coefficient_function,
    variance_components  = std_result$variance_components,
    predictions          = predictions,
    residuals            = std_result$residuals,
    basis                = basis,
    data                 = dt,
    engine               = engine,
    call                 = mc,
    extras               = list(
      model_type    = "joint",
      raw_model     = raw_result[["model"]],
      stacked_data  = stacked,
      mm_reparam    = mm_reparam,
      int_basis     = int_basis,
      trap_weights  = trap_weights,
      convergence   = std_result$convergence
    )
  )

  .msg("Joint model fitting complete.")
  model_out
}


# ---- Relationship matrix validation helper ------------------------------------

#' Validate a relationship matrix against group levels
#' @param mat Numeric matrix to validate.
#' @param levels Character vector of expected row/column names.
#' @param name Character label for error messages.
#' @return Invisible TRUE on success.
#' @noRd
.validate_relationship_matrix <- function(mat, levels, name) {
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop(sprintf("`%s` must be a numeric matrix.", name), call. = FALSE)
  }
  if (nrow(mat) != ncol(mat)) {
    stop(sprintf("`%s` must be square.", name), call. = FALSE)
  }
  if (!isSymmetric(unname(mat), tol = 1e-8)) {
    stop(sprintf("`%s` must be symmetric.", name), call. = FALSE)
  }
  mat_names <- rownames(mat)
  if (!is.null(mat_names)) {
    missing_lvl <- setdiff(levels, mat_names)
    if (length(missing_lvl) > 0L) {
      stop(
        sprintf(
          "`%s` row/column names do not cover all group levels. Missing: %s.",
          name, paste(head(missing_lvl, 5L), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}


# ==============================================================================
# fit_2d_functional -- 2D functional profiles via tensor product B-splines
# ==============================================================================

#' Fit 2D Functional Profiles via Tensor Product B-Splines
#'
#' Models a functional trait observed over two continuous indices (e.g., time
#' and soil depth, time and wavelength) using tensor product B-splines with
#' separable penalties. Group-specific (e.g., variety-specific) 2D surfaces
#' are estimated as random effects.
#'
#' @param data A `data.frame` or [data.table::data.table] in long format with
#'   one row per observation (id x dim1 x dim2).
#' @param time_col Character; column name for the first functional dimension
#'   (e.g., time). Default `"time"`.
#' @param depth_col Character; column name for the second functional dimension
#'   (e.g., depth, wavelength). Default `"depth"`.
#' @param value_col Character; column name for the response. Default `"value"`.
#' @param id_col Character; column name for the observational unit. Default
#'   `"id"`.
#' @param group_col Character; column name for the grouping factor
#'   (e.g., variety). Default `"variety"`.
#' @param n_knots1 Integer; number of internal knots for the first dimension.
#'   Default 10.
#' @param n_knots2 Integer; number of internal knots for the second dimension.
#'   Default 10.
#' @param degree1 Integer; B-spline degree for dimension 1. Default 3.
#' @param degree2 Integer; B-spline degree for dimension 2. Default 3.
#' @param spatial_row_col Character or `NULL`; spatial row column. Default
#'   `NULL`.
#' @param spatial_col_col Character or `NULL`; spatial column column. Default
#'   `NULL`.
#' @param engine Character; estimation engine. Default `"auto"`.
#' @param ... Additional arguments passed to the backend.
#'
#' @return An `fda_model` object with `model_type = "2d_functional"` in
#'   extras. The `fitted_curves` component contains fitted 2D surfaces
#'   per group, and `basis` contains the `fda_tensor_basis` object.
#'
#' @section Statistical Model:
#' \deqn{y(t, d) = f_v(t, d) + \text{spatial} + \varepsilon}
#' where \eqn{f_v(t, d) = \sum_{j} \sum_{k} \alpha_{vjk} B_{1j}(t) B_{2k}(d)}
#' is a variety-specific tensor product random surface, and the penalty is
#' the standard Kronecker sum:
#' \eqn{P = P_1 \otimes I_2 + I_1 \otimes P_2}.
#'
#' @references
#' Marx, B.D. and Eilers, P.H.C. (2005). Multidimensional penalized signal
#' regression. *Technometrics*, 47(1), 13--22.
#'
#' Wood, S.N. (2006). Low-rank scale-invariant tensor product smooths for
#' generalized additive mixed models. *Biometrics*, 62(4), 1025--1036.
#'
#' @seealso [tensor_bspline_basis()] for tensor product basis construction,
#'   [scalar_on_2d_function()] for relating 2D profiles to a scalar response.
#'
#' @examples
#' \dontrun{
#' # Simulated 2D functional data: time x depth
#' set.seed(123)
#' n_obs <- 500
#' dt2d <- data.table::data.table(
#'   id      = rep(paste0("P", 1:10), each = 50),
#'   variety = rep(rep(c("V1", "V2"), each = 5), each = 50),
#'   time    = runif(n_obs, 0, 100),
#'   depth   = runif(n_obs, 0, 1),
#'   value   = rnorm(n_obs)
#' )
#'
#' fit2d <- fit_2d_functional(
#'   data      = dt2d,
#'   time_col  = "time",
#'   depth_col = "depth",
#'   value_col = "value",
#'   id_col    = "id",
#'   group_col = "variety",
#'   n_knots1  = 6, n_knots2 = 5
#' )
#' }
#'
#' @export
fit_2d_functional <- function(
    data,
    time_col        = "time",
    depth_col       = "depth",
    value_col       = "value",
    id_col          = "id",
    group_col       = "variety",
    n_knots1        = 10L,
    n_knots2        = 10L,
    degree1         = 3L,
    degree2         = 3L,
    spatial_row_col = NULL,
    spatial_col_col = NULL,
    engine          = "auto",
    ...) {

  mc <- match.call()

  # ---- Input validation -------------------------------------------------------

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or data.table.", call. = FALSE)
  }
  dt <- data.table::as.data.table(data)

  required_cols <- c(time_col, depth_col, value_col, id_col, group_col)
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("Column(s) not found in data: %s.",
              paste0("'", missing_cols, "'", collapse = ", ")),
      call. = FALSE
    )
  }

  .validate_numeric(dt[[time_col]], time_col)
  .validate_numeric(dt[[depth_col]], depth_col)
  .validate_numeric(dt[[value_col]], value_col)

  engine <- .resolve_engine(engine)

  # ---- Construct tensor product basis -----------------------------------------

  tensor_basis <- tensor_bspline_basis(
    x1         = dt[[time_col]],
    x2         = dt[[depth_col]],
    n_knots1   = n_knots1,
    n_knots2   = n_knots2,
    degree1    = degree1,
    degree2    = degree2
  )

  .msg(sprintf(
    "Tensor product basis: %d x %d = %d total basis functions.",
    tensor_basis$basis1$n_basis, tensor_basis$basis2$n_basis,
    tensor_basis$n_basis
  ))

  # ---- Prepare design matrices ------------------------------------------------

  B_tensor <- tensor_basis$B
  n_obs    <- nrow(dt)
  n_tbasis <- tensor_basis$n_basis

  dt[, (group_col) := as.factor(dt[[group_col]])]
  dt[, (id_col)    := as.factor(dt[[id_col]])]

  group_levels <- levels(dt[[group_col]])
  n_groups     <- length(group_levels)

  # Add tensor B-spline columns to data
  tb_col_names <- paste0("tb_", seq_len(n_tbasis))
  for (k in seq_len(n_tbasis)) {
    data.table::set(dt, j = tb_col_names[k], value = B_tensor[, k])
  }

  # ---- Build model specification ----------------------------------------------

  # Fixed effects: intercept (and group if desired)
  fixed_formula <- stats::as.formula(
    sprintf("%s ~ 1", value_col)
  )

  # Random effects: group-specific tensor product random effects
  # Variety x tensor-basis coefficients with separable penalty
  random_parts <- group_col

  random_formula <- stats::as.formula(
    paste("~", paste(random_parts, collapse = " + "))
  )

  # Residual structure
  if (!is.null(spatial_row_col) && !is.null(spatial_col_col)) {
    rcov_formula <- stats::as.formula(
      sprintf("~ ar1(%s):ar1(%s)", spatial_row_col, spatial_col_col)
    )
  } else {
    rcov_formula <- stats::as.formula("~ units")
  }

  # Known matrices: penalty matrices for the tensor product
  known_mats <- list(
    P_tensor  = as.matrix(tensor_basis$P),
    P1_tensor = as.matrix(tensor_basis$P1),
    P2_tensor = as.matrix(tensor_basis$P2)
  )

  model_spec <- list(
    fixed          = fixed_formula,
    random         = random_formula,
    rcov           = rcov_formula,
    known_matrices = known_mats
  )

  dots <- list(...)
  if (length(dots) > 0L) {
    model_spec <- c(model_spec, dots)
  }

  # ---- Dispatch ---------------------------------------------------------------

  .msg(sprintf(
    "Fitting 2D functional model: %d obs, %d groups, %d tensor basis, engine = '%s'.",
    n_obs, n_groups, n_tbasis, engine
  ))

  raw_result <- .dispatch_fit(
    engine     = engine,
    model_spec = model_spec,
    data       = dt,
    ...
  )

  # ---- Standardise and return -------------------------------------------------

  std_result <- .standardise_result(
    raw_result = raw_result,
    engine     = engine,
    model_spec = model_spec,
    basis      = tensor_basis,
    data       = dt
  )

  predictions <- data.table::data.table(
    group     = group_levels,
    predicted = rep(NA_real_, n_groups),
    se        = rep(NA_real_, n_groups)
  )

  model_out <- new_fda_model(
    fitted_curves        = std_result$fitted_curves,
    coefficient_function = std_result$coefficient_function,
    variance_components  = std_result$variance_components,
    predictions          = predictions,
    residuals            = std_result$residuals,
    basis                = tensor_basis,
    data                 = dt,
    engine               = engine,
    call                 = mc,
    extras               = list(
      model_type  = "2d_functional",
      raw_model   = raw_result[["model"]],
      convergence = std_result$convergence
    )
  )

  .msg("2D functional model fitting complete.")
  model_out
}


# ==============================================================================
# scalar_on_2d_function -- Scalar-on-2D-function regression
# ==============================================================================

#' Scalar-on-2D-Function Regression
#'
#' Extends scalar-on-function regression to functional predictors observed over
#' two continuous dimensions (e.g., time x depth). The model relates a scalar
#' primary trait to a 2D functional predictor via double integration:
#'
#' \deqn{y_i = \alpha + \int_t \int_d \beta(t, d) f_i(t, d)\, dt\, dd
#'   + \varepsilon_i}
#'
#' where \eqn{\beta(t, d)} is a 2D coefficient surface estimated via tensor
#' product B-splines with penalised regression, and \eqn{f_i(t, d)} are
#' the group-specific 2D functional profiles.
#'
#' @param primary_trait A named numeric vector of primary trait values (one per
#'   group/variety) or a [data.table::data.table] with columns `group` and
#'   `value`.
#' @param functional_profiles_2d An `fda_model` object from
#'   [fit_2d_functional()], or a list containing `fitted_surfaces` (a
#'   [data.table::data.table] with columns `group`, `time`, `depth`, `fitted`)
#'   and `basis` (an `fda_tensor_basis` object).
#' @param basis An `fda_tensor_basis` object. If `NULL`, extracted from
#'   `functional_profiles_2d`.
#' @param genomic_matrix Optional genomic relationship matrix. Default `NULL`.
#' @param engine Character; estimation engine. Default `"auto"`.
#' @param ... Additional arguments passed to the backend.
#'
#' @return An `fda_model` object where:
#' \describe{
#'   \item{coefficient_function}{Contains the estimated 2D coefficient surface
#'     \eqn{\beta(t, d)} as a data.table with columns `time`, `depth`, `beta`,
#'     `se`, `ci_lower`, `ci_upper`.}
#'   \item{predictions}{Predicted primary trait values for each group.}
#' }
#'
#' @references
#' Marx, B.D. and Eilers, P.H.C. (2005). Multidimensional penalized signal
#' regression. *Technometrics*, 47(1), 13--22.
#'
#' @seealso [fit_2d_functional()] for 2D profile estimation,
#'   [scalar_on_function()] for the 1D equivalent.
#'
#' @examples
#' \dontrun{
#' # After fitting 2D profiles:
#' # fit2d <- fit_2d_functional(...)
#' # primary_vals <- c(V1 = 4.2, V2 = 3.8, ...)
#' # sof2d <- scalar_on_2d_function(primary_vals, fit2d)
#' }
#'
#' @export
scalar_on_2d_function <- function(
    primary_trait,
    functional_profiles_2d,
    basis           = NULL,
    genomic_matrix  = NULL,
    engine          = "auto",
    ...) {

  mc <- match.call()
  engine <- .resolve_engine(engine)

  # ---- Validate and extract inputs --------------------------------------------

  # Extract basis

  if (is.null(basis)) {
    if (inherits(functional_profiles_2d, "fda_model")) {
      basis <- functional_profiles_2d$basis
    } else if (is.list(functional_profiles_2d) &&
               !is.null(functional_profiles_2d$basis)) {
      basis <- functional_profiles_2d$basis
    } else {
      stop(
        "Cannot extract tensor basis. Supply `basis` or an `fda_model` object.",
        call. = FALSE
      )
    }
  }

  if (!inherits(basis, "fda_tensor_basis")) {
    stop("`basis` must be an 'fda_tensor_basis' object.", call. = FALSE)
  }

  # Parse primary_trait
  if (is.data.frame(primary_trait)) {
    pt_dt <- data.table::as.data.table(primary_trait)
    if (!all(c("group", "value") %in% names(pt_dt))) {
      stop(
        "`primary_trait` data.frame must have columns 'group' and 'value'.",
        call. = FALSE
      )
    }
    group_names <- as.character(pt_dt$group)
    y_primary   <- pt_dt$value
  } else if (is.numeric(primary_trait) && !is.null(names(primary_trait))) {
    group_names <- names(primary_trait)
    y_primary   <- unname(primary_trait)
  } else {
    stop(
      "`primary_trait` must be a named numeric vector or a data.frame ",
      "with columns 'group' and 'value'.",
      call. = FALSE
    )
  }

  n_groups <- length(group_names)
  .validate_numeric(y_primary, "primary_trait values")

  # ---- Extract BLUPs from 2D functional profiles ------------------------------

  # Get the estimated B-spline coefficient vectors per group
  if (inherits(functional_profiles_2d, "fda_model")) {
    fitted_dt <- functional_profiles_2d$fitted_curves
  } else if (is.list(functional_profiles_2d) &&
             !is.null(functional_profiles_2d$fitted_surfaces)) {
    fitted_dt <- functional_profiles_2d$fitted_surfaces
  } else {
    stop(
      "Cannot extract fitted surfaces from `functional_profiles_2d`.",
      call. = FALSE
    )
  }

  # ---- Compute integrated functional covariates (double integral) -------------

  basis1 <- basis$basis1
  basis2 <- basis$basis2

  # Trapezoidal weights for dimension 1
  t1_unique <- sort(unique(basis1$x))
  n_t1 <- length(t1_unique)
  w1 <- numeric(n_t1)
  if (n_t1 >= 2L) {
    dt1 <- diff(t1_unique)
    w1[1L] <- dt1[1L] / 2
    w1[n_t1] <- dt1[n_t1 - 1L] / 2
    if (n_t1 > 2L) {
      for (k in 2L:(n_t1 - 1L)) w1[k] <- (dt1[k - 1L] + dt1[k]) / 2
    }
  } else {
    w1 <- 1
  }

  # Trapezoidal weights for dimension 2
  t2_unique <- sort(unique(basis2$x))
  n_t2 <- length(t2_unique)
  w2 <- numeric(n_t2)
  if (n_t2 >= 2L) {
    dt2 <- diff(t2_unique)
    w2[1L] <- dt2[1L] / 2
    w2[n_t2] <- dt2[n_t2 - 1L] / 2
    if (n_t2 > 2L) {
      for (k in 2L:(n_t2 - 1L)) w2[k] <- (dt2[k - 1L] + dt2[k]) / 2
    }
  } else {
    w2 <- 1
  }

  # 2D integration weights: outer product of 1D weights
  # Evaluate tensor basis at a grid of (t1, t2) points for integration
  grid_dt <- data.table::CJ(t1 = t1_unique, t2 = t2_unique)
  n_grid  <- nrow(grid_dt)

  B1_grid <- splines::splineDesign(
    knots = c(rep(basis1$boundary[1L], basis1$degree + 1L),
              basis1$knots,
              rep(basis1$boundary[2L], basis1$degree + 1L)),
    x     = grid_dt$t1,
    ord   = basis1$degree + 1L,
    outer.ok = TRUE
  )

  B2_grid <- splines::splineDesign(
    knots = c(rep(basis2$boundary[1L], basis2$degree + 1L),
              basis2$knots,
              rep(basis2$boundary[2L], basis2$degree + 1L)),
    x     = grid_dt$t2,
    ord   = basis2$degree + 1L,
    outer.ok = TRUE
  )

  # Tensor product basis at grid points
  B_tensor_grid <- .row_kronecker(B1_grid, B2_grid)

  # 2D integration weights
  w_2d <- rep(w1, each = n_t2) * rep(w2, times = n_t1)

  # Integrated basis: each column integrated over the 2D domain
  int_tensor_basis <- as.numeric(crossprod(B_tensor_grid, w_2d))

  n_tbasis <- basis$n_basis

  # ---- Build scalar regression design matrix ----------------------------------

  # X_func[i, k] = integral of B_k(t,d) * f_i(t,d) dt dd
  # approximated via: sum_grid w_2d * B_k(grid) * f_i(grid)
  # This requires the fitted surface values at grid points for each group.
  # For simplicity, we use the integrated tensor basis as the functional
  # covariate (analogous to the 1D case).

  # Create data for scalar regression
  dt_scalar <- data.table::data.table(
    group   = factor(group_names),
    primary = y_primary
  )

  # Add integrated functional covariates
  for (k in seq_len(n_tbasis)) {
    data.table::set(dt_scalar, j = paste0("int_tb_", k),
                    value = int_tensor_basis[k])
  }

  # ---- Build model specification and fit --------------------------------------

  intb_names <- paste0("int_tb_", seq_len(n_tbasis))
  fixed_formula <- stats::as.formula(
    sprintf("primary ~ %s", paste(intb_names, collapse = " + "))
  )

  random_formula <- stats::as.formula("~ group")

  known_mats <- list(P_tensor = as.matrix(basis$P))
  if (!is.null(genomic_matrix)) {
    known_mats[["G_genomic"]] <- genomic_matrix
  }

  model_spec <- list(
    fixed          = fixed_formula,
    random         = random_formula,
    rcov           = stats::as.formula("~ units"),
    known_matrices = known_mats
  )

  dots <- list(...)
  if (length(dots) > 0L) model_spec <- c(model_spec, dots)

  .msg(sprintf(
    "Fitting scalar-on-2D-function regression: %d groups, %d tensor basis, engine = '%s'.",
    n_groups, n_tbasis, engine
  ))

  raw_result <- .dispatch_fit(
    engine     = engine,
    model_spec = model_spec,
    data       = dt_scalar,
    ...
  )

  # ---- Extract 2D coefficient surface -----------------------------------------

  std_result <- .standardise_result(
    raw_result = raw_result,
    engine     = engine,
    model_spec = model_spec,
    basis      = basis,
    data       = dt_scalar
  )

  # Construct the 2D coefficient surface on a fine grid
  n_grid_fine <- 50L
  t1_fine <- seq(basis1$boundary[1L], basis1$boundary[2L],
                 length.out = n_grid_fine)
  t2_fine <- seq(basis2$boundary[1L], basis2$boundary[2L],
                 length.out = n_grid_fine)
  grid_fine <- data.table::CJ(time = t1_fine, depth = t2_fine)

  # Evaluate tensor basis on the fine grid
  B1_fine <- splines::splineDesign(
    knots = c(rep(basis1$boundary[1L], basis1$degree + 1L),
              basis1$knots,
              rep(basis1$boundary[2L], basis1$degree + 1L)),
    x     = grid_fine$time,
    ord   = basis1$degree + 1L,
    outer.ok = TRUE
  )

  B2_fine <- splines::splineDesign(
    knots = c(rep(basis2$boundary[1L], basis2$degree + 1L),
              basis2$knots,
              rep(basis2$boundary[2L], basis2$degree + 1L)),
    x     = grid_fine$depth,
    ord   = basis2$degree + 1L,
    outer.ok = TRUE
  )

  B_fine <- .row_kronecker(B1_fine, B2_fine)

  # Attempt to extract beta coefficients from raw model
  coef_2d <- list(
    time     = grid_fine$time,
    depth    = grid_fine$depth,
    beta     = rep(NA_real_, nrow(grid_fine)),
    se       = rep(NA_real_, nrow(grid_fine)),
    ci_lower = rep(NA_real_, nrow(grid_fine)),
    ci_upper = rep(NA_real_, nrow(grid_fine))
  )

  beta_coefs <- .extract_fixed_spline_coefs(raw_result, engine, n_tbasis)
  if (!is.null(beta_coefs)) {
    coef_2d$beta <- as.numeric(B_fine %*% beta_coefs$estimate)
    if (!is.null(beta_coefs$se)) {
      coef_2d$se <- sqrt(as.numeric(B_fine^2 %*% beta_coefs$se^2))
      coef_2d$ci_lower <- coef_2d$beta - 1.96 * coef_2d$se
      coef_2d$ci_upper <- coef_2d$beta + 1.96 * coef_2d$se
    }
  }

  # Predictions
  predictions <- data.table::data.table(
    group     = group_names,
    predicted = rep(NA_real_, n_groups),
    se        = rep(NA_real_, n_groups)
  )

  model_out <- new_fda_model(
    fitted_curves        = std_result$fitted_curves,
    coefficient_function = coef_2d,
    variance_components  = std_result$variance_components,
    predictions          = predictions,
    residuals            = std_result$residuals,
    basis                = basis,
    data                 = dt_scalar,
    engine               = engine,
    call                 = mc,
    extras               = list(
      model_type       = "scalar_on_2d_function",
      raw_model        = raw_result[["model"]],
      int_tensor_basis = int_tensor_basis,
      convergence      = std_result$convergence
    )
  )

  .msg("Scalar-on-2D-function regression complete.")
  model_out
}
