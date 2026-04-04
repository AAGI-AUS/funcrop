# scalar_on_function.R -- Stage 2: Scalar-on-function regression
#
# Relates a scalar primary trait (e.g., yield) to the functional profiles
# recovered in Stage 1 (fit_functional_profiles). This is the second stage of
# the two-stage FDA approach.
#
# Model: yield_i = alpha + integral( beta(t) * f_hat_i(t) ) dt + u_i + eps_i
# where beta(t) is the coefficient function estimated via B-spline expansion.
#
# Also provides compare_methods() for comparing multiple fda_model fits.
#
# Author: Maksym Gaidashenko
# Licence: GPL (>= 3)


# ---- scalar_on_function ----

#' Scalar-on-Function Regression (Stage 2)
#'
#' Fits a scalar-on-function regression model that relates a scalar primary
#' trait (e.g., yield) to the variety-specific functional profiles estimated
#' in Stage 1 (via [fit_functional_profiles()]). The key output is the
#' coefficient function \eqn{\beta(t)}, which identifies the time periods
#' during which the secondary trait most strongly predicts the primary trait.
#'
#' @param primary_trait One of:
#'   - A named numeric vector with names corresponding to variety identifiers.
#'   - An `fda_data` object containing a `primary_trait` column (extracted
#'     automatically).
#'   - A data.frame / data.table with columns `variety` and `yield` (or
#'     the name specified via the first column).
#' @param functional_profiles One of:
#'   - An `fda_model` object from [fit_functional_profiles()] (preferred).
#'   - A numeric matrix of B-spline coefficients (varieties x basis functions),
#'     in which case `basis` must also be supplied.
#' @param basis An `fda_basis` object (from [bspline_basis()]). Required when
#'   `functional_profiles` is a raw matrix; ignored when it is an `fda_model`
#'   (the basis is extracted from the model object).
#' @param genomic_matrix Optional numeric matrix (symmetric, PSD) of genomic
#'   relationships among varieties. Row/column names must match variety
#'   identifiers. Used as `vm(variety, source = G)` in the mixed model.
#' @param pedigree_matrix Optional numeric matrix (symmetric, PSD) of
#'   pedigree-based relationships. Same naming convention as `genomic_matrix`.
#' @param environment_col Character or `NULL`; name of the environment column
#'   in `primary_trait` if multi-environment data is supplied. Default `NULL`
#'   (single-environment analysis).
#' @param engine Character; estimation engine. One of `"auto"` (default),
#'   `"asreml"`, or `"bayesreml"`. See [funcrop_engines()].
#' @param ... Additional arguments passed to the backend fitting function.
#'
#' @return An `fda_model` object containing:
#' \describe{
#'   \item{fitted_curves}{Empty data.table (not applicable for Stage 2).}
#'   \item{coefficient_function}{List with `time`, `beta`, `se`, `ci_lower`,
#'     `ci_upper` -- the estimated coefficient function \eqn{\hat{\beta}(t)}
#'     evaluated on a fine grid.}
#'   \item{variance_components}{data.table of estimated variance components.}
#'   \item{predictions}{data.table of predicted primary trait values per
#'     variety.}
#'   \item{residuals}{Model residuals.}
#'   \item{basis}{The `fda_basis` object used.}
#'   \item{data}{The input `primary_trait` data.}
#'   \item{engine}{Character string identifying the engine used.}
#'   \item{call}{The matched function call.}
#'   \item{extras}{List with `beta_coefs` (estimated B-spline coefficients for
#'     \eqn{\beta(t)}), `J_matrix` (inner product matrix of B-spline basis),
#'     `C_matrix` (functional covariate matrix), `raw_model`.}
#' }
#'
#' @details
#' ## Model formulation
#'
#' \deqn{y_i = \alpha + \int \beta(t) \hat{f}_i(t) \, dt + u_i + \varepsilon_i}
#'
#' where:
#' - \eqn{\hat{f}_i(t)} is the Stage 1 estimated functional profile for
#'   variety \eqn{i}.
#' - \eqn{\beta(t) = \sum_k b_k B_k(t)} is the coefficient function,
#'   expanded in the same B-spline basis.
#' - The integral is approximated using the inner product matrix:
#'   \deqn{J_{jk} = \int B_j(t) B_k(t) \, dt}
#'   computed via numerical quadrature (composite Simpson's rule).
#' - The functional covariate matrix is:
#'   \deqn{C = \hat{A} J}
#'   where \eqn{\hat{A}} is the variety x basis coefficient matrix from
#'   Stage 1. Each row of \eqn{C} is the functional covariate for one
#'   variety.
#' - \eqn{u_i} is an optional genomic/pedigree random effect:
#'   \eqn{u \sim N(0, \sigma_g^2 G)} where \eqn{G} is the relationship
#'   matrix.
#' - \eqn{\varepsilon_i \sim N(0, \sigma^2)} is the residual.
#'
#' ## Smoothness of beta(t)
#'
#' The coefficient function is regularised via a P-spline penalty on the
#' B-spline coefficients \eqn{b}: the penalty matrix from the basis is
#' incorporated as a random effect precision structure (mixed model
#' representation of P-splines).
#'
#' @references
#' Ramsay, J.O. and Silverman, B.W. (2005). *Functional Data Analysis*
#' (2nd ed.). Springer.
#'
#' Reiss, P.T. and Ogden, R.T. (2007). Functional principal component
#' regression and functional partial least squares. *Journal of the American
#' Statistical Association*, 102(479), 984--996.
#'
#' Marx, B.D. and Eilers, P.H.C. (1999). Generalized linear regression on
#' sampled signals and curves: a P-spline approach. *Technometrics*, 41(1),
#' 1--13.
#'
#' @seealso [fit_functional_profiles()] for Stage 1, [compare_methods()] for
#'   model comparison, [bspline_basis()] for basis construction.
#'
#' @examples
#' \dontrun{
#' # Assuming `profiles` is an fda_model from fit_functional_profiles()
#' data(sim_grain_fill)
#'
#' # Named yield vector
#' yield_vec <- setNames(
#'   sim_grain_fill[, .(yield = mean(yield)), by = variety]$yield,
#'   sim_grain_fill[, .(yield = mean(yield)), by = variety]$variety
#' )
#'
#' # Stage 2: scalar-on-function regression
#' sof_fit <- scalar_on_function(
#'   primary_trait        = yield_vec,
#'   functional_profiles  = profiles,
#'   engine               = "auto"
#' )
#'
#' # Inspect the coefficient function
#' plot(sof_fit, which = "coef")
#' coef(sof_fit)
#' }
#'
#' @export
scalar_on_function <- function(
    primary_trait,
    functional_profiles,
    basis              = NULL,
    genomic_matrix     = NULL,
    pedigree_matrix    = NULL,
    environment_col    = NULL,
    engine             = "auto",
    ...
) {

  call <- match.call()

  # ===========================================================================
  # Input validation
  # ===========================================================================

  if (missing(primary_trait) || is.null(primary_trait)) {
    stop("`primary_trait` must be supplied.", call. = FALSE)
  }
  if (missing(functional_profiles) || is.null(functional_profiles)) {
    stop("`functional_profiles` must be supplied.", call. = FALSE)
  }

  # --- Resolve engine ---
  engine <- .resolve_engine(engine)

  # ===========================================================================
  # Parse primary_trait
  # ===========================================================================

  if (inherits(primary_trait, "fda_data")) {
    meta <- attr(primary_trait, "fda_meta")
    if (!meta$has_primary_trait) {
      stop(
        "`primary_trait` is an fda_data object but does not contain a ",
        "primary_trait column.", call. = FALSE
      )
    }
    pt_dt <- as.data.table.fda_data(primary_trait)
    # Extract one value per group (variety)
    group_col <- if ("group" %in% names(pt_dt)) "group" else NULL
    if (is.null(group_col)) {
      stop(
        "fda_data object must have a `group` column (variety identifier) ",
        "for scalar-on-function regression.", call. = FALSE
      )
    }
    yield_dt <- unique(pt_dt[, .(variety = as.character(get(group_col)),
                                  yield = primary_trait)])
    # Check uniqueness
    if (anyDuplicated(yield_dt$variety)) {
      # Take the mean per variety
      yield_dt <- yield_dt[, .(yield = mean(yield, na.rm = TRUE)),
                            by = variety]
    }
    yield_vec <- stats::setNames(yield_dt$yield, yield_dt$variety)
  } else if (is.numeric(primary_trait) && !is.null(names(primary_trait))) {
    # Named numeric vector
    yield_vec <- primary_trait
    if (anyNA(yield_vec)) {
      n_na <- sum(is.na(yield_vec))
      warning(
        sprintf("Removing %d NA values from `primary_trait`.", n_na),
        call. = FALSE
      )
      yield_vec <- yield_vec[!is.na(yield_vec)]
    }
    if (length(yield_vec) < 3L) {
      stop("`primary_trait` must have at least 3 non-NA values.", call. = FALSE)
    }
  } else if (is.data.frame(primary_trait)) {
    pt_dt <- data.table::as.data.table(primary_trait)
    # Expect columns: first = variety, second = yield (or named)
    if (ncol(pt_dt) < 2L) {
      stop(
        "`primary_trait` data.frame must have at least 2 columns ",
        "(variety identifier and yield).", call. = FALSE
      )
    }
    variety_cn <- names(pt_dt)[1L]
    yield_cn   <- names(pt_dt)[2L]
    yield_vec  <- stats::setNames(pt_dt[[yield_cn]],
                                   as.character(pt_dt[[variety_cn]]))
  } else {
    stop(
      "`primary_trait` must be a named numeric vector, fda_data object, ",
      "or data.frame.", call. = FALSE
    )
  }

  # ===========================================================================
  # Parse functional_profiles
  # ===========================================================================

  if (inherits(functional_profiles, "fda_model")) {
    # Extract from Stage 1 fda_model
    if (is.null(functional_profiles$extras$spline_blups)) {
      stop(
        "`functional_profiles` fda_model does not contain spline_blups in ",
        "extras. Was it produced by fit_functional_profiles()?", call. = FALSE
      )
    }
    alpha_matrix <- functional_profiles$extras$spline_blups
    basis <- functional_profiles$basis

    if (!inherits(basis, "fda_basis")) {
      stop(
        "Could not extract a valid fda_basis from the functional_profiles model.",
        call. = FALSE
      )
    }
  } else if (is.matrix(functional_profiles) && is.numeric(functional_profiles)) {
    alpha_matrix <- functional_profiles
    if (is.null(basis)) {
      stop(
        "`basis` must be supplied when `functional_profiles` is a raw matrix.",
        call. = FALSE
      )
    }
    if (!inherits(basis, "fda_basis")) {
      stop("`basis` must be an 'fda_basis' object from bspline_basis().",
           call. = FALSE)
    }
  } else {
    stop(
      "`functional_profiles` must be an fda_model from ",
      "fit_functional_profiles() or a numeric matrix of B-spline coefficients.",
      call. = FALSE
    )
  }

  # Validate dimensions
  if (is.null(rownames(alpha_matrix))) {
    warning(
      "`alpha_matrix` has no rownames; assuming rows match `primary_trait` names.",
      call. = FALSE
    )
  }

  # --- Align varieties between yield and alpha_matrix ---
  alpha_varieties <- rownames(alpha_matrix)
  yield_varieties <- names(yield_vec)

  if (!is.null(alpha_varieties) && !is.null(yield_varieties)) {
    common <- intersect(alpha_varieties, yield_varieties)
    if (length(common) == 0L) {
      stop(
        "No common variety identifiers between `primary_trait` and ",
        "`functional_profiles`.", call. = FALSE
      )
    }
    if (length(common) < length(alpha_varieties) ||
        length(common) < length(yield_varieties)) {
      n_dropped <- max(length(alpha_varieties), length(yield_varieties)) -
        length(common)
      .msg("Aligning varieties: ", length(common), " in common, ",
           n_dropped, " dropped.")
    }
    alpha_matrix <- alpha_matrix[common, , drop = FALSE]
    yield_vec    <- yield_vec[common]
  } else if (nrow(alpha_matrix) != length(yield_vec)) {
    stop(
      sprintf(
        "Number of varieties in `functional_profiles` (%d) does not match ",
        "`primary_trait` (%d) and no names available for alignment.",
        nrow(alpha_matrix), length(yield_vec)
      ),
      call. = FALSE
    )
  }

  n_var   <- nrow(alpha_matrix)
  n_coefs <- ncol(alpha_matrix)

  if (n_var < 3L) {
    stop("At least 3 varieties are required for scalar-on-function regression.",
         call. = FALSE)
  }

  # --- Validate optional relationship matrices ---
  if (!is.null(genomic_matrix)) {
    .validate_relationship_matrix(genomic_matrix, names(yield_vec), "genomic_matrix")
  }
  if (!is.null(pedigree_matrix)) {
    .validate_relationship_matrix(pedigree_matrix, names(yield_vec), "pedigree_matrix")
  }

  # ===========================================================================
  # Compute the inner product matrix J
  # ===========================================================================

  # J_{jk} = integral( B_j(t) * B_k(t) ) dt over the basis domain.
  # Computed via composite Simpson's rule on a fine grid.
  J <- .compute_bspline_inner_product(basis)

  # ===========================================================================
  # Compute the functional covariate matrix C
  # ===========================================================================

  # Note: alpha_matrix from Stage 1 is in the transformed (Z) space if
  # make_Zspline was used. We need to map back to the full basis or use the
  # appropriate inner product.
  #
  # If alpha_matrix columns = n_z_cols (random effect space), then the
  # "basis functions" are not the raw B-splines but the transformed ones.
  # We compute C using the transformed basis inner products.

  if (n_coefs == basis$n_basis) {
    # Alpha is in original B-spline space: C = Alpha %*% J
    C <- alpha_matrix %*% J
  } else {
    # Alpha is in the Z-space (from make_Zspline decomposition).
    # The Z-space basis functions are Z_fine = B %*% U_range %*% diag(1/sqrt(lambda)).
    # Inner product in Z-space: J_z = integral( Z_j(t) * Z_k(t) ) dt
    J_z <- .compute_z_space_inner_product(basis)
    if (ncol(J_z) != n_coefs) {
      # Dimension mismatch -- fall back to using what we have
      warning(
        sprintf(
          "Z-space inner product matrix has %d columns but alpha_matrix has %d. ",
          ncol(J_z), n_coefs
        ),
        "Using truncated/padded version.",
        call. = FALSE
      )
      k <- min(ncol(J_z), n_coefs)
      J_z <- J_z[seq_len(k), seq_len(k), drop = FALSE]
      alpha_matrix <- alpha_matrix[, seq_len(k), drop = FALSE]
    }
    C <- alpha_matrix %*% J_z
  }

  # C is now n_var x n_basis (or n_var x n_z_cols) -- the functional covariate.
  # Each column of C corresponds to a B-spline coefficient b_k in beta(t).

  # ===========================================================================
  # Build model data for Stage 2
  # ===========================================================================

  variety_ids <- names(yield_vec) %||% rownames(alpha_matrix) %||%
    paste0("V", seq_len(n_var))

  model_dt <- data.table::data.table(
    variety = factor(variety_ids, levels = variety_ids),
    yield   = as.numeric(yield_vec)
  )

  # Add C matrix columns as covariates
  n_c_cols <- ncol(C)
  # Column names: use letters to avoid ASReml parsing digits as indices
  c_letters <- c(letters, paste0(rep(letters, each = 26), letters))
  c_names <- paste0("fc", c_letters[seq_len(n_c_cols)])
  for (k in seq_len(n_c_cols)) {
    data.table::set(model_dt, j = c_names[k], value = C[, k])
  }

  # ===========================================================================
  # Build model specification
  # ===========================================================================

  # Fixed effects: intercept + C columns (with penalty for smoothness)
  # In the mixed model representation of P-splines, the smooth beta(t) is

  # split into fixed (null space) and random (range space) components.
  # For Stage 2, we use the C columns as fixed effects with a penalty on
  # the b coefficients -- implemented via a random effect with known precision.

  # Apply the same decomposition to the beta(t) coefficients:
  # beta(t) = X_null %*% alpha_fixed + Z_range %*% u_beta
  # where u_beta ~ N(0, sigma_beta^2 * I)

  # For simplicity: fit C columns as fixed effects first, then add penalty
  # as random effect structure if the backend supports it.
  fixed_rhs <- paste(c_names, collapse = " + ")
  fixed_formula <- stats::as.formula(paste("yield ~", fixed_rhs))

  # Random effects
  random_terms <- character(0L)
  known_mats   <- list()

  # Genomic random effect
  if (!is.null(genomic_matrix)) {
    # Subset and reorder to match variety_ids
    G_aligned <- genomic_matrix[variety_ids, variety_ids]
    known_mats[["G_genomic"]] <- G_aligned
    random_terms <- c(random_terms,
                      "vm(variety, source = G_genomic)")
  }

  # Pedigree random effect
  if (!is.null(pedigree_matrix)) {
    A_aligned <- pedigree_matrix[variety_ids, variety_ids]
    known_mats[["A_pedigree"]] <- A_aligned
    random_terms <- c(random_terms,
                      "vm(variety, source = A_pedigree)")
  }

  # Smoothness penalty on beta(t) via random effect representation
  # The penalty matrix P operates on the b coefficients.
  # In the mixed model: some b's are fixed (null space), rest are random.
  # For the Stage 2 model with C as covariates, the penalty translates to:
  # a ridge-type penalty on the C-column regression coefficients.
  # This is implemented as: the C columns corresponding to the range space

  # of the penalty are modelled as random effects with identity covariance.
  beta_decomp <- make_Zspline(basis, constraint = "decompose")
  n_null_beta <- beta_decomp$null_dim
  n_rand_beta <- beta_decomp$rank

  # Split the C matrix into null-space (fixed) and range-space (random) parts
  if (n_c_cols == basis$n_basis) {
    # C is in original basis space -- apply the decomposition to columns
    P_dense <- as.matrix(basis$P)
    eig <- eigen(P_dense, symmetric = TRUE)
    tol <- max(eig$values) * .Machine$double.eps * basis$n_basis * 10
    idx_range <- which(eig$values > tol)
    idx_null  <- which(eig$values <= tol)
    U_null  <- eig$vectors[, idx_null, drop = FALSE]
    U_range <- eig$vectors[, idx_range, drop = FALSE]
    scaling <- 1 / sqrt(eig$values[idx_range])

    # C_fixed = C %*% U_null  (null space: unpenalised)
    # C_random = C %*% U_range %*% diag(scaling)  (range space: penalised)
    C_fixed  <- C %*% U_null
    C_random <- C %*% (U_range %*% diag(scaling, nrow = length(scaling)))

    # Replace columns in model_dt
    # Remove old C columns
    for (cn in c_names) {
      data.table::set(model_dt, j = cn, value = NULL)
    }

    # Add null-space fixed columns
    fn_letters <- c(letters, paste0(rep(letters, each = 26), letters))
    fixed_c_names <- paste0("fn", fn_letters[seq_len(ncol(C_fixed))])
    for (k in seq_len(ncol(C_fixed))) {
      data.table::set(model_dt, j = fixed_c_names[k], value = C_fixed[, k])
    }

    # Add range-space random columns
    rn_letters <- c(letters, paste0(rep(letters, each = 26), letters))
    rand_c_names <- paste0("rn", rn_letters[seq_len(ncol(C_random))])
    for (k in seq_len(ncol(C_random))) {
      data.table::set(model_dt, j = rand_c_names[k], value = C_random[, k])
    }

    fixed_rhs <- paste(fixed_c_names, collapse = " + ")
    fixed_formula <- stats::as.formula(paste("yield ~", fixed_rhs))

    # The random C columns represent the penalised part of beta(t)
    # Each Cr_k column is a continuous covariate treated as a random regression
    # coefficient. In both ASReml and bayesreml, these are independent random
    # effects with identity covariance (iid shrinkage on the wiggly part of
    # beta(t)). This is the same approach as Model 4 in the illustration.
    random_terms <- c(random_terms, rand_c_names)
  } else {
    # C is in Z-space -- treat all C columns as fixed (penalty already absorbed
    # by the Z-space transformation from Stage 1)
    .msg("Functional covariates in Z-space: treating as fixed effects.")
  }

  # Build random formula
  if (length(random_terms) > 0L) {
    random_formula <- stats::as.formula(paste("~", paste(random_terms,
                                                          collapse = " + ")))
  } else {
    random_formula <- NULL
  }

  model_spec <- list(
    fixed          = fixed_formula,
    random         = random_formula,
    rcov           = NULL,
    known_matrices = known_mats
  )

  # ===========================================================================
  # Dispatch to backend
  # ===========================================================================

  .msg("Fitting scalar-on-function model via '", engine, "' engine...")
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

  # -- Reconstruct beta(t) coefficient function --
  coef_fn <- .reconstruct_beta_function(
    raw_result  = raw_result,
    engine      = engine,
    basis       = basis,
    n_c_cols    = n_c_cols,
    fixed_c_names = if (exists("fixed_c_names")) fixed_c_names else c_names,
    rand_c_names  = if (exists("rand_c_names")) rand_c_names else NULL,
    eig           = if (exists("eig")) eig else NULL,
    idx_null      = if (exists("idx_null")) idx_null else NULL,
    idx_range     = if (exists("idx_range")) idx_range else NULL,
    scaling       = if (exists("scaling")) scaling else NULL
  )

  # -- Predictions for primary trait --
  predictions <- data.table::data.table(
    variety   = variety_ids,
    observed  = as.numeric(yield_vec),
    predicted = NA_real_
  )
  # Try to extract fitted values
  fitted_vals <- tryCatch(
    stats::fitted(raw_result[["model"]]),
    error = function(e) NULL
  )
  if (!is.null(fitted_vals) && length(fitted_vals) == n_var) {
    predictions[, predicted := fitted_vals]
  }

  # ===========================================================================
  # Construct and return fda_model
  # ===========================================================================

  new_fda_model(
    fitted_curves          = data.table::data.table(),
    coefficient_function   = coef_fn,
    variance_components    = std_result$variance_components,
    predictions            = predictions,
    residuals              = std_result$residuals,
    basis                  = basis,
    data                   = model_dt,
    engine                 = engine,
    call                   = call,
    extras                 = list(
      beta_coefs    = coef_fn$beta_coefs,
      J_matrix      = if (exists("J")) J else NULL,
      C_matrix      = C,
      alpha_matrix  = alpha_matrix,
      raw_model     = std_result$raw_model,
      convergence   = std_result$convergence,
      model_spec    = model_spec
    )
  )
}


# ---- compare_methods ----

#' Compare Multiple FDA Model Fits
#'
#' Takes multiple `fda_model` objects and computes comparison metrics including
#' information criteria (AIC, BIC, WAIC where available), prediction accuracy
#' (RMSE, R-squared), and coverage of confidence/credible intervals.
#'
#' @param ... Two or more `fda_model` objects to compare.
#' @param labels Optional character vector of model labels. If `NULL`,
#'   labels are generated from the call or assigned sequentially.
#'
#' @return An `fda_comparison` object (see `new_fda_comparison()`) containing:
#' \describe{
#'   \item{models}{Named list of the input `fda_model` objects.}
#'   \item{metrics}{data.table of comparison metrics with one row per model
#'     and columns: `model`, `engine`, `aic`, `bic`, `waic`, `rmse`,
#'     `r_squared`, `coverage_95`.}
#'   \item{label}{Character description of the comparison.}
#' }
#'
#' @details
#' Metrics computed:
#' \describe{
#'   \item{AIC, BIC}{Extracted from the fitted model where available (ASReml).
#'     `NA` for Bayesian models (use WAIC instead).}
#'   \item{WAIC}{Widely applicable information criterion. Available for
#'     bayesreml models only.}
#'   \item{RMSE}{Root mean squared error of predictions vs observed primary
#'     trait values.}
#'   \item{R-squared}{Coefficient of determination for the primary trait
#'     predictions.}
#'   \item{Coverage}{Proportion of observed values falling within the 95\%
#'     confidence/credible intervals of the fitted curves (where available).}
#' }
#'
#' @seealso [fit_functional_profiles()], [scalar_on_function()],
#'   `new_fda_comparison()`.
#'
#' @examples
#' \dontrun{
#' # Compare two models with different basis sizes
#' comp <- compare_methods(model_10k, model_20k,
#'                         labels = c("10 knots", "20 knots"))
#' print(comp)
#' summary(comp)
#' plot(comp)
#' }
#'
#' @export
compare_methods <- function(..., labels = NULL) {

  models <- list(...)

  # ===========================================================================
  # Input validation
  # ===========================================================================

  if (length(models) < 2L) {
    stop("At least 2 fda_model objects must be supplied for comparison.",
         call. = FALSE)
  }

  # Check all are fda_model
  is_model <- vapply(models, inherits, logical(1L), "fda_model")
  if (!all(is_model)) {
    bad_idx <- which(!is_model)
    stop(
      sprintf(
        "All arguments must be fda_model objects. Arguments at positions %s are not.",
        paste(bad_idx, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # --- Labels ---
  n_models <- length(models)
  if (is.null(labels)) {
    # Try to extract from call
    cl <- match.call(expand.dots = FALSE)
    dot_args <- cl[["..."]]
    if (!is.null(dot_args)) {
      labels <- vapply(dot_args, deparse, character(1L))
    } else {
      labels <- paste0("model_", seq_len(n_models))
    }
  }

  if (!is.character(labels) || length(labels) != n_models) {
    stop(
      sprintf("`labels` must be a character vector of length %d.", n_models),
      call. = FALSE
    )
  }

  if (anyDuplicated(labels)) {
    labels <- make.unique(labels, sep = "_")
  }

  names(models) <- labels

  # ===========================================================================
  # Compute metrics
  # ===========================================================================

  metrics_list <- vector("list", n_models)

  for (i in seq_len(n_models)) {
    m <- models[[i]]
    lbl <- labels[i]

    # --- Information criteria ---
    aic_val  <- .safe_extract_ic(m, "aic")
    bic_val  <- .safe_extract_ic(m, "bic")
    waic_val <- .safe_extract_ic(m, "waic")

    # --- Prediction accuracy ---
    rmse_val <- NA_real_
    r2_val   <- NA_real_

    preds <- m$predictions
    if (is.data.table(preds) && nrow(preds) > 0L) {
      obs_col  <- intersect(names(preds), c("observed", "yield", "primary_trait"))
      pred_col <- intersect(names(preds), c("predicted", "fitted", "pred"))

      if (length(obs_col) > 0L && length(pred_col) > 0L) {
        obs_vals  <- preds[[obs_col[1L]]]
        pred_vals <- preds[[pred_col[1L]]]
        valid     <- !is.na(obs_vals) & !is.na(pred_vals)

        if (sum(valid) >= 3L) {
          resid    <- obs_vals[valid] - pred_vals[valid]
          rmse_val <- sqrt(mean(resid^2))
          ss_res   <- sum(resid^2)
          ss_tot   <- sum((obs_vals[valid] - mean(obs_vals[valid]))^2)
          r2_val   <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
        }
      }
    }

    # --- Coverage (95% CI) ---
    coverage_val <- .compute_coverage(m)

    metrics_list[[i]] <- data.table::data.table(
      model       = lbl,
      engine      = m$engine,
      aic         = aic_val,
      bic         = bic_val,
      waic        = waic_val,
      rmse        = rmse_val,
      r_squared   = r2_val,
      coverage_95 = coverage_val
    )
  }

  metrics <- data.table::rbindlist(metrics_list, use.names = TRUE)

  # ===========================================================================
  # Construct and return fda_comparison
  # ===========================================================================

  comp_label <- paste0(
    "Comparison of ", n_models, " models (",
    paste(unique(metrics$engine), collapse = " vs "), ")"
  )

  new_fda_comparison(
    models  = models,
    metrics = metrics,
    label   = comp_label
  )
}


# ==============================================================================
# Internal helper functions
# ==============================================================================

#' Compute inner product matrix of B-spline basis functions
#'
#' Computes \eqn{J_{jk} = \int B_j(t) B_k(t) \, dt} via composite Simpson's
#' rule on a fine grid. This matrix is needed for the scalar-on-function
#' integral approximation.
#'
#' @param basis An `fda_basis` object.
#' @param n_quad Integer; number of quadrature points. Default 500.
#'
#' @return Symmetric numeric matrix (n_basis x n_basis).
#' @noRd
.compute_bspline_inner_product <- function(basis, n_quad = 500L) {
  if (!inherits(basis, "fda_basis")) {
    stop("`basis` must be an fda_basis object.", call. = FALSE)
  }

  # Fine quadrature grid
  t_quad <- seq(basis$boundary[1L], basis$boundary[2L], length.out = n_quad)
  h <- diff(t_quad[1:2])

  # Evaluate basis at quadrature points

  knot_vec <- c(
    rep(basis$boundary[1L], basis$degree + 1L),
    basis$knots,
    rep(basis$boundary[2L], basis$degree + 1L)
  )
  B_quad <- splines::splineDesign(
    knots    = knot_vec,
    x        = t_quad,
    ord      = basis$degree + 1L,
    outer.ok = TRUE
  )

  # Simpson's rule weights: 1, 4, 2, 4, 2, ..., 4, 1 (for odd n_quad)
  # Fall back to trapezoidal if n_quad is even
  if (n_quad %% 2L == 1L) {
    w <- rep(c(2, 4), length.out = n_quad)
    w[1L] <- 1
    w[n_quad] <- 1
    w <- w * h / 3
  } else {
    # Trapezoidal rule
    w <- rep(h, n_quad)
    w[1L] <- h / 2
    w[n_quad] <- h / 2
  }

  # J = B' %*% diag(w) %*% B = (B * sqrt(w))' %*% (B * sqrt(w))
  # More efficient: J = crossprod(B_quad * sqrt(w))
  B_weighted <- B_quad * sqrt(w)
  J <- crossprod(B_weighted)

  J
}


#' Compute inner product matrix in the Z-space (reparameterised basis)
#'
#' For the reparameterised basis Z = B %*% U_range %*% diag(1/sqrt(lambda)),
#' computes the inner product matrix in this transformed space.
#'
#' @param basis An `fda_basis` object.
#' @param n_quad Integer; number of quadrature points. Default 500.
#'
#' @return Symmetric numeric matrix (n_z_cols x n_z_cols).
#' @noRd
.compute_z_space_inner_product <- function(basis, n_quad = 500L) {
  if (!inherits(basis, "fda_basis")) {
    stop("`basis` must be an fda_basis object.", call. = FALSE)
  }

  # Quadrature grid
  t_quad <- seq(basis$boundary[1L], basis$boundary[2L], length.out = n_quad)
  h <- diff(t_quad[1:2])

  # Evaluate basis at quadrature points
  knot_vec <- c(
    rep(basis$boundary[1L], basis$degree + 1L),
    basis$knots,
    rep(basis$boundary[2L], basis$degree + 1L)
  )
  B_quad <- splines::splineDesign(
    knots    = knot_vec,
    x        = t_quad,
    ord      = basis$degree + 1L,
    outer.ok = TRUE
  )

  # Eigendecomposition for Z-space transformation
  P_dense <- as.matrix(basis$P)
  eig <- eigen(P_dense, symmetric = TRUE)
  tol <- max(eig$values) * .Machine$double.eps * basis$n_basis * 10
  idx_range <- which(eig$values > tol)
  U_range <- eig$vectors[, idx_range, drop = FALSE]
  scaling <- 1 / sqrt(eig$values[idx_range])

  Z_quad <- B_quad %*% (U_range %*% diag(scaling, nrow = length(scaling)))

  # Simpson's / trapezoidal weights
  if (n_quad %% 2L == 1L) {
    w <- rep(c(2, 4), length.out = n_quad)
    w[1L] <- 1
    w[n_quad] <- 1
    w <- w * h / 3
  } else {
    w <- rep(h, n_quad)
    w[1L] <- h / 2
    w[n_quad] <- h / 2
  }

  Z_weighted <- Z_quad * sqrt(w)
  J_z <- crossprod(Z_weighted)

  J_z
}


#' Validate a relationship matrix (genomic or pedigree)
#'
#' Checks that the matrix is square, symmetric, has matching row/column names,
#' and covers the required variety identifiers.
#'
#' @param mat Numeric matrix to validate.
#' @param variety_ids Character vector of required variety identifiers.
#' @param name Name of the argument for error messages.
#' @noRd
.validate_relationship_matrix <- function(mat, variety_ids, name) {
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop(sprintf("`%s` must be a numeric matrix.", name), call. = FALSE)
  }
  if (nrow(mat) != ncol(mat)) {
    stop(sprintf("`%s` must be a square matrix.", name), call. = FALSE)
  }
  if (!isSymmetric(mat, tol = 1e-6)) {
    stop(sprintf("`%s` must be symmetric.", name), call. = FALSE)
  }
  rn <- rownames(mat)
  if (is.null(rn)) {
    stop(sprintf("`%s` must have row and column names.", name), call. = FALSE)
  }
  missing_ids <- setdiff(variety_ids, rn)
  if (length(missing_ids) > 0L) {
    stop(
      sprintf(
        "`%s` is missing entries for varieties: %s.",
        name, paste(head(missing_ids, 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


#' Reconstruct the coefficient function beta(t) from Stage 2 model
#'
#' Combines the fixed and random effect estimates for the C columns to
#' reconstruct the coefficient function on a fine time grid.
#'
#' @param raw_result Raw backend result.
#' @param engine Character; engine name.
#' @param basis fda_basis object.
#' @param n_c_cols Number of C columns.
#' @param fixed_c_names Names of fixed C columns.
#' @param rand_c_names Names of random C columns (NULL if none).
#' @param eig Eigendecomposition result (NULL if C is in original space).
#' @param idx_null Null-space eigenvalue indices.
#' @param idx_range Range-space eigenvalue indices.
#' @param scaling Eigenvalue scaling factors.
#'
#' @return List with `time`, `beta`, `se`, `ci_lower`, `ci_upper`,
#'   `beta_coefs`.
#' @noRd
.reconstruct_beta_function <- function(raw_result, engine, basis, n_c_cols,
                                       fixed_c_names, rand_c_names, eig,
                                       idx_null, idx_range, scaling) {
  model <- raw_result[["model"]]
  n_basis <- basis$n_basis

  # Fine evaluation grid
  n_grid <- 200L
  t_grid <- seq(basis$boundary[1L], basis$boundary[2L], length.out = n_grid)

  knot_vec <- c(
    rep(basis$boundary[1L], basis$degree + 1L),
    basis$knots,
    rep(basis$boundary[2L], basis$degree + 1L)
  )
  B_grid <- splines::splineDesign(
    knots    = knot_vec,
    x        = t_grid,
    ord      = basis$degree + 1L,
    outer.ok = TRUE
  )

  # Try to extract fixed effect coefficients
  beta_raw <- tryCatch({
    coefs <- summary(model, coef = TRUE)
    if (is.list(coefs) && !is.null(coefs$coef.fixed)) {
      coefs$coef.fixed
    } else if (is.data.frame(coefs)) {
      coefs
    } else {
      stats::coef(model)
    }
  }, error = function(e) {
    tryCatch(stats::coef(model), error = function(e2) NULL)
  })

  # Default output
  beta_coefs <- rep(0, n_basis)
  se_coefs   <- rep(NA_real_, n_basis)

  if (!is.null(beta_raw)) {
    # Extract coefficients corresponding to the C columns
    if (is.matrix(beta_raw) || is.data.frame(beta_raw)) {
      coef_names <- rownames(beta_raw)
      coef_vals  <- beta_raw[, 1L]
      coef_se    <- if (ncol(beta_raw) >= 2L) beta_raw[, 2L] else NULL
    } else if (is.numeric(beta_raw)) {
      coef_names <- names(beta_raw)
      coef_vals  <- beta_raw
      coef_se    <- NULL
    } else {
      coef_names <- NULL
      coef_vals  <- NULL
      coef_se    <- NULL
    }

    if (!is.null(coef_names) && !is.null(coef_vals)) {
      # Match fixed C coefficients
      fixed_idx <- match(fixed_c_names, coef_names)
      fixed_idx <- fixed_idx[!is.na(fixed_idx)]

      if (length(fixed_idx) > 0L && !is.null(eig) && !is.null(idx_null)) {
        # Reconstruct in original basis space
        b_null <- coef_vals[fixed_idx]
        if (length(b_null) == length(idx_null)) {
          # Map back: beta_orig = U_null %*% b_null + U_range %*% diag(scaling) %*% b_range
          U_null  <- eig$vectors[, idx_null, drop = FALSE]
          U_range <- eig$vectors[, idx_range, drop = FALSE]
          beta_coefs <- as.numeric(U_null %*% b_null)

          # TODO: extract random part (b_range) from BLUPs if available
          # For now, set range component to zero (conservative)
        }
      } else if (length(fixed_idx) > 0L) {
        # C is in original space -- coefficients are beta directly
        beta_coefs <- coef_vals[fixed_idx]
        if (!is.null(coef_se)) {
          se_coefs <- coef_se[fixed_idx]
        }
      }
    }
  }

  # Evaluate beta(t) = B_grid %*% beta_coefs
  if (length(beta_coefs) != n_basis) {
    # Pad or truncate
    beta_coefs <- c(beta_coefs, rep(0, max(0L, n_basis - length(beta_coefs))))
    beta_coefs <- beta_coefs[seq_len(n_basis)]
  }

  beta_vals <- as.numeric(B_grid %*% beta_coefs)

  # Standard error propagation
  if (!all(is.na(se_coefs)) && length(se_coefs) == n_basis) {
    se_vals <- sqrt(as.numeric(B_grid^2 %*% se_coefs^2))
  } else {
    se_vals <- rep(NA_real_, n_grid)
  }

  # Confidence intervals
  if (!all(is.na(se_vals))) {
    ci_lower <- beta_vals - 1.96 * se_vals
    ci_upper <- beta_vals + 1.96 * se_vals
  } else {
    ci_lower <- rep(NA_real_, n_grid)
    ci_upper <- rep(NA_real_, n_grid)
  }

  list(
    time       = t_grid,
    beta       = beta_vals,
    se         = se_vals,
    ci_lower   = ci_lower,
    ci_upper   = ci_upper,
    beta_coefs = beta_coefs
  )
}


#' Safely extract information criterion from an fda_model
#'
#' @param model An `fda_model` object.
#' @param criterion Character: `"aic"`, `"bic"`, or `"waic"`.
#' @return Numeric value or `NA_real_`.
#' @noRd
.safe_extract_ic <- function(model, criterion = c("aic", "bic", "waic")) {
  criterion <- match.arg(criterion)
  raw_model <- model$extras$raw_model

  if (is.null(raw_model)) return(NA_real_)

  tryCatch({
    switch(criterion,
      aic = {
        if (model$engine == "asreml") {
          # ASReml: summary has loglik; AIC = -2*loglik + 2*npar
          smry <- summary(raw_model)
          loglik <- smry$loglik %||% raw_model$loglik
          npar <- length(raw_model$vparameters)
          if (!is.null(loglik)) -2 * loglik + 2 * npar else NA_real_
        } else {
          stats::AIC(raw_model)
        }
      },
      bic = {
        if (model$engine == "asreml") {
          smry <- summary(raw_model)
          loglik <- smry$loglik %||% raw_model$loglik
          npar <- length(raw_model$vparameters)
          nobs <- nrow(model$data)
          if (!is.null(loglik)) -2 * loglik + log(nobs) * npar else NA_real_
        } else {
          stats::BIC(raw_model)
        }
      },
      waic = {
        if (model$engine == "bayesreml") {
          # bayesreml v0.1.0: WAIC stored in fit$extras$summary or similar
          w <- tryCatch({
            waic_val <- raw_model$extras$summary$waic
            if (is.null(waic_val)) NA_real_ else waic_val
          }, error = function(e) NA_real_)
          if (is.list(w)) w$waic else w
        } else {
          NA_real_
        }
      }
    )
  }, error = function(e) NA_real_)
}


#' Compute 95% coverage for fitted curves
#'
#' For models with CI columns in fitted_curves, computes the proportion of
#' observed data points falling within the 95% CI.
#'
#' @param model An `fda_model` object.
#' @return Numeric coverage proportion, or `NA_real_` if CIs are not available.
#' @noRd
.compute_coverage <- function(model) {
  fc <- model$fitted_curves

  if (!is.data.table(fc) || nrow(fc) == 0L) return(NA_real_)

  if (!all(c("ci_lower", "ci_upper") %in% names(fc))) return(NA_real_)
  if (all(is.na(fc$ci_lower)) || all(is.na(fc$ci_upper))) return(NA_real_)

  # Match fitted curves to observed data
  obs_data <- model$data
  if (!is.data.frame(obs_data)) return(NA_real_)

  # Find value column
  val_col <- intersect(names(obs_data), c("value", "grain_weight", "ndvi"))
  if (length(val_col) == 0L) return(NA_real_)

  # Merge on id + time (approximate match for fine-grid curves)
  if (!all(c("id", "time") %in% names(obs_data))) return(NA_real_)

  # Use nearest-neighbour matching on the time grid
  obs_dt <- data.table::as.data.table(obs_data)
  fc_valid <- fc[!is.na(ci_lower) & !is.na(ci_upper)]

  if (nrow(fc_valid) == 0L) return(NA_real_)

  # For each observed (id, time), find the nearest fitted curve point
  n_covered <- 0L
  n_total   <- 0L

  for (uid in unique(obs_dt$id)) {
    obs_sub <- obs_dt[id == uid]
    fc_sub  <- fc_valid[id == uid]
    if (nrow(fc_sub) == 0L) next

    for (r in seq_len(nrow(obs_sub))) {
      t_obs <- obs_sub$time[r]
      y_obs <- obs_sub[[val_col[1L]]][r]
      if (is.na(y_obs)) next

      # Nearest time in fitted curve
      nearest_idx <- which.min(abs(fc_sub$time - t_obs))
      ci_lo <- fc_sub$ci_lower[nearest_idx]
      ci_hi <- fc_sub$ci_upper[nearest_idx]

      n_total <- n_total + 1L
      if (y_obs >= ci_lo && y_obs <= ci_hi) {
        n_covered <- n_covered + 1L
      }
    }
  }

  if (n_total == 0L) return(NA_real_)
  n_covered / n_total
}
