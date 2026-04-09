# funcrop: Shared result standardisation across backends
#
# Internal functions that take the raw output from either the ASReml or
# bayesreml backend and produce a standardised result structure. This
# ensures all downstream code (S3 methods, plotting, inference) operates
# on a common representation regardless of the estimation engine used.

# ---- Main standardisation function -------------------------------------------

#' Standardise raw backend result into funcrop's internal format
#'
#' Takes the raw output from `.asreml_fit()` or `.bayesreml_fit()` and creates
#' a standardised named list containing all quantities needed to construct an
#' `fda_model` S3 object.
#'
#' @param raw_result Named list returned by `.asreml_fit()` or
#'   `.bayesreml_fit()`.
#' @param engine Character string: `"asreml"` or `"bayesreml"`.
#' @param model_spec Named list with the original model specification (formulas,
#'   known matrices, etc.).
#' @param basis The funcrop basis object used in the model (contains B-spline
#'   evaluation matrix, knots, penalty).
#' @param data The original data (data.frame or data.table) passed to the
#'   fitting function.
#'
#' @return Named list suitable for `new_fda_model()`:
#'   \describe{
#'     \item{`fitted_curves`}{data.table of fitted functional curves per group.}
#'     \item{`coefficient_function`}{data.table of the estimated coefficient
#'       function beta(t) on a fine grid.}
#'     \item{`variance_components`}{data.table of estimated variance
#'       components.}
#'     \item{`blups`}{data.table of BLUPs/posterior means for random effects.}
#'     \item{`residuals`}{Numeric vector of model residuals.}
#'     \item{`convergence`}{List with convergence diagnostics.}
#'     \item{`engine`}{Character string identifying the engine used.}
#'     \item{`raw_model`}{The original fitted model object.}
#'   }
#' @noRd
.standardise_result <- function(raw_result, engine, model_spec, basis, data) {
  # --- Variance components ---
  vc <- switch(
    engine,
    asreml    = .asreml_extract_vc(raw_result[["model"]]),
    bayesreml = .bayesreml_extract_vc(raw_result[["model"]]),
    lme4      = .lme4_extract_vc(raw_result[["model"]]),
    mgcv      = .mgcv_extract_vc(raw_result[["model"]]),
    stop(sprintf("Unknown engine '%s'.", engine), call. = FALSE)
  )

  # --- BLUPs for random terms ---
  random_terms <- .extract_random_term_names(model_spec[["random"]])
  blups <- switch(
    engine,
    asreml    = .asreml_extract_blups(raw_result[["model"]], random_terms),
    bayesreml = .bayesreml_extract_blups(raw_result[["model"]], random_terms),
    lme4      = .lme4_extract_blups(raw_result[["model"]], random_terms),
    mgcv      = .mgcv_extract_blups(raw_result[["model"]], random_terms)
  )

  # --- Fitted curves (variety-specific functional profiles) ---
  fitted_curves <- .extract_fitted_curves(raw_result, engine, basis, data)

  # --- Coefficient function beta(t) ---
  coef_function <- .extract_coefficient_function(raw_result, engine, basis)

  # --- Residuals ---
  resid_vec <- tryCatch(
    stats::residuals(raw_result[["model"]]),
    error = function(e) {
      warning("Could not extract residuals: ", conditionMessage(e),
              call. = FALSE)
      rep(NA_real_, nrow(data))
    }
  )

  # --- Convergence diagnostics ---
  convergence <- list(
    converged = raw_result[["converged"]],
    engine    = engine
  )
  if (engine == "asreml") {
    convergence[["log_lik"]] <- raw_result[["log_lik"]]
    convergence[["n_iter"]]  <- raw_result[["n_iter"]]
  }
  if (engine == "bayesreml") {
    convergence[["summary"]] <- raw_result[["summary"]]
  }

  list(
    fitted_curves         = fitted_curves,
    coefficient_function  = coef_function,
    variance_components   = vc,
    blups                 = blups,
    residuals             = resid_vec,
    convergence           = convergence,
    engine                = engine,
    raw_model             = raw_result[["model"]]
  )
}


# ---- Fitted curves extraction ------------------------------------------------

#' Reconstruct fitted variety curves from BLUPs and B-spline basis
#'
#' For each variety (or group), computes the fitted functional curve by
#' combining the estimated spline coefficients (BLUPs or posterior means)
#' with the B-spline basis functions:
#' \deqn{\hat{u}_v(t) = \sum_{k=1}^{K} \hat{\alpha}_{vk} B_k(t)}
#'
#' The curve is evaluated on a fine grid of time points for smooth
#' visualisation.
#'
#' @param raw_result Named list from the backend.
#' @param engine Character: `"asreml"` or `"bayesreml"`.
#' @param basis funcrop basis object containing `Bmatrix` and `knots`.
#' @param data Original data (used to identify group/variety structure).
#'
#' @return data.table with columns:
#'   \describe{
#'     \item{`id`}{Group/variety identifier.}
#'     \item{`time`}{Time point (fine grid).}
#'     \item{`fitted`}{Fitted curve value at that time point.}
#'     \item{`se`}{Standard error (ASReml) or posterior SD (bayesreml).}
#'     \item{`ci_lower`}{Lower 95\% CI (bayesreml only; NA for ASReml).}
#'     \item{`ci_upper`}{Upper 95\% CI (bayesreml only; NA for ASReml).}
#'   }
#' @noRd
.extract_fitted_curves <- function(raw_result, engine, basis, data) {
  # Determine evaluation grid: use fine grid spanning the basis knot range
  if (is.list(basis) && !is.null(basis[["knots"]])) {
    knot_range <- range(basis[["knots"]])
  } else {
    # Fallback: use data range (assumes a 'time' column)
    if ("time" %in% names(data)) {
      knot_range <- range(data[["time"]], na.rm = TRUE)
    } else {
      warning("Cannot determine time range for fitted curves.", call. = FALSE)
      return(data.table::data.table())
    }
  }

  # Fine evaluation grid (200 points for smooth curves)
  n_grid <- 200L
  t_grid <- seq(knot_range[1L], knot_range[2L], length.out = n_grid)

  # Evaluate B-spline basis at the fine grid
  if (is.list(basis) && !is.null(basis[["degree"]])) {
    # Use splines::bs() or splineDesign() to evaluate at new points
    B_grid <- splines::bs(
      x       = t_grid,
      knots   = basis[["internal_knots"]],
      degree  = basis[["degree"]],
      Boundary.knots = knot_range
    )
  } else if (is.list(basis) && !is.null(basis[["Bmatrix"]])) {
    # Attempt to re-evaluate using stored basis info
    # If basis evaluation function is stored, use it; otherwise interpolate
    warning(
      "Basis re-evaluation at fine grid requires knot information. ",
      "Using original evaluation points.",
      call. = FALSE
    )
    B_grid <- basis[["Bmatrix"]]
    t_grid <- if (!is.null(basis[["time_points"]])) {
      basis[["time_points"]]
    } else {
      seq_len(nrow(B_grid))
    }
    n_grid <- length(t_grid)
  } else {
    warning("Cannot evaluate basis at fine grid.", call. = FALSE)
    return(data.table::data.table())
  }

  # Extract spline coefficient BLUPs per variety
  # The BLUP data.table has columns: term, level, blup, se (and possibly
  # ci_lower, ci_upper for bayesreml)
  spline_terms <- .find_spline_terms(raw_result, engine)

  if (length(spline_terms) == 0L) {
    .msg("No spline random effect terms found; skipping fitted curves.")
    return(data.table::data.table())
  }

  blup_dt <- switch(
    engine,
    asreml    = .asreml_extract_blups(raw_result[["model"]], spline_terms),
    bayesreml = .bayesreml_extract_blups(raw_result[["model"]], spline_terms),
    lme4      = .lme4_extract_blups(raw_result[["model"]], spline_terms),
    mgcv      = .mgcv_extract_blups(raw_result[["model"]], spline_terms),
    data.table::data.table()
  )

  if (is.null(blup_dt) || !is.data.frame(blup_dt) || nrow(blup_dt) == 0L) {
    return(data.table::data.table())
  }

  # Parse variety:spline_coef levels into variety ID and coefficient index.
  # Convention: level format is "variety_name:coef_k" or just sequential
  # within each variety.
  n_basis <- ncol(B_grid)
  variety_ids <- unique(blup_dt[["level"]])

  # If levels encode variety:coefficient, parse them
  if (any(grepl(":", blup_dt[["level"]]))) {
    blup_dt[, c("variety", "coef_idx") := data.table::tstrsplit(
      level, ":", fixed = TRUE
    )]
    blup_dt[, coef_idx := as.integer(coef_idx)]
    variety_ids <- unique(blup_dt[["variety"]])
  } else {
    # Assume levels are ordered: first n_basis belong to variety 1, etc.
    n_levels <- nrow(blup_dt)
    if (n_levels %% n_basis != 0L) {
      warning(
        "Number of BLUPs is not a multiple of the number of basis functions. ",
        "Cannot reconstruct fitted curves.",
        call. = FALSE
      )
      return(data.table::data.table())
    }
    n_varieties <- n_levels %/% n_basis
    blup_dt[, variety  := rep(seq_len(n_varieties), each = n_basis)]
    blup_dt[, coef_idx := rep(seq_len(n_basis), times = n_varieties)]
    variety_ids <- unique(blup_dt[["variety"]])
  }

  # Reconstruct fitted curves: fitted_v(t) = B_grid %*% alpha_v
  curve_list <- vector("list", length(variety_ids))

  for (j in seq_along(variety_ids)) {
    vid <- variety_ids[j]
    alpha_v <- blup_dt[variety == vid, blup]

    if (length(alpha_v) != n_basis) {
      warning(
        sprintf(
          "Variety '%s' has %d coefficients but basis has %d functions. Skipping.",
          as.character(vid), length(alpha_v), n_basis
        ),
        call. = FALSE
      )
      next
    }

    fitted_v <- as.numeric(B_grid %*% alpha_v)

    # Standard error / uncertainty propagation
    # For ASReml: SE of linear combination via delta method (approximate)
    # For bayesreml: propagate posterior SD through basis (approximate)
    se_v <- blup_dt[variety == vid, ]
    if ("se" %in% names(se_v)) {
      # Approximate SE of curve: SE(fitted_v(t)) ~= sqrt(B_grid^2 %*% se_v^2)
      # This assumes independence of coefficient estimates (conservative approx)
      se_alpha <- se_v[["se"]]
      se_fitted <- sqrt(as.numeric(B_grid^2 %*% se_alpha^2))
    } else {
      se_fitted <- rep(NA_real_, n_grid)
    }

    # Credible intervals (bayesreml only)
    ci_lower <- rep(NA_real_, n_grid)
    ci_upper <- rep(NA_real_, n_grid)
    if (engine == "bayesreml" && all(c("ci_lower", "ci_upper") %in%
                                     names(se_v))) {
      # Approximate propagation for CIs (same basis-coefficient approach)
      ci_lo_alpha <- se_v[["ci_lower"]]
      ci_hi_alpha <- se_v[["ci_upper"]]
      # Note: this is an approximation. Proper CI propagation requires
      # posterior draws of the full curve. For now, use linear propagation.
      ci_lower <- as.numeric(B_grid %*% ci_lo_alpha)
      ci_upper <- as.numeric(B_grid %*% ci_hi_alpha)
    }

    curve_list[[j]] <- data.table::data.table(
      id       = as.character(vid),
      time     = t_grid,
      fitted   = fitted_v,
      se       = se_fitted,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    )
  }

  data.table::rbindlist(curve_list, use.names = TRUE)
}


# ---- Coefficient function extraction -----------------------------------------

#' Extract the coefficient function beta(t)
#'
#' Reconstructs the estimated coefficient function from the fixed effect
#' B-spline coefficients:
#' \deqn{\hat{\beta}(t) = \sum_{k=1}^{K} \hat{\beta}_k B_k(t)}
#'
#' This is the key output for scalar-on-function regression, representing
#' the time-varying effect of the functional predictor on the scalar response.
#'
#' @param raw_result Named list from the backend.
#' @param engine Character: `"asreml"` or `"bayesreml"`.
#' @param basis funcrop basis object.
#'
#' @return data.table with columns:
#'   \describe{
#'     \item{`time`}{Evaluation time points (fine grid).}
#'     \item{`beta`}{Estimated coefficient function value.}
#'     \item{`se`}{Standard error (ASReml) or posterior SD (bayesreml).}
#'     \item{`ci_lower`}{Lower 95\% CI (bayesreml; derived from SE for ASReml).}
#'     \item{`ci_upper`}{Upper 95\% CI (bayesreml; derived from SE for ASReml).}
#'   }
#' @noRd
.extract_coefficient_function <- function(raw_result, engine, basis) {
  # Determine evaluation grid
  if (is.list(basis) && !is.null(basis[["knots"]])) {
    knot_range <- range(basis[["knots"]])
  } else {
    warning(
      "Cannot determine time range for coefficient function evaluation.",
      call. = FALSE
    )
    return(data.table::data.table())
  }

  n_grid <- 200L
  t_grid <- seq(knot_range[1L], knot_range[2L], length.out = n_grid)

  # Evaluate B-spline basis at the fine grid
  if (is.list(basis) && !is.null(basis[["degree"]]) &&
      !is.null(basis[["internal_knots"]])) {
    B_grid <- splines::bs(
      x              = t_grid,
      knots          = basis[["internal_knots"]],
      degree         = basis[["degree"]],
      Boundary.knots = knot_range
    )
  } else {
    warning(
      "Insufficient basis information for coefficient function evaluation.",
      call. = FALSE
    )
    return(data.table::data.table())
  }

  n_basis <- ncol(B_grid)

  # Extract fixed effect spline coefficients
  beta_coefs <- .extract_fixed_spline_coefs(raw_result, engine, n_basis)

  if (is.null(beta_coefs)) {
    .msg("No fixed-effect spline coefficients found; ",
         "skipping coefficient function.")
    return(data.table::data.table())
  }

  # Reconstruct beta(t) = B_grid %*% beta_k
  beta_vals <- as.numeric(B_grid %*% beta_coefs[["estimate"]])

  # Standard error propagation
  if (!is.null(beta_coefs[["se"]])) {
    se_vals <- sqrt(as.numeric(B_grid^2 %*% beta_coefs[["se"]]^2))
  } else {
    se_vals <- rep(NA_real_, n_grid)
  }

  # Confidence / credible intervals
  if (engine == "bayesreml" && !is.null(beta_coefs[["ci_lower"]])) {
    ci_lower <- as.numeric(B_grid %*% beta_coefs[["ci_lower"]])
    ci_upper <- as.numeric(B_grid %*% beta_coefs[["ci_upper"]])
  } else if (!all(is.na(se_vals))) {
    # Wald-type 95% CI for ASReml
    ci_lower <- beta_vals - 1.96 * se_vals
    ci_upper <- beta_vals + 1.96 * se_vals
  } else {
    ci_lower <- rep(NA_real_, n_grid)
    ci_upper <- rep(NA_real_, n_grid)
  }

  data.table::data.table(
    time     = t_grid,
    beta     = beta_vals,
    se       = se_vals,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}


# ---- Helper functions --------------------------------------------------------

#' Extract random term names from a formula
#'
#' Parses a random effects formula (one-sided, e.g., `~ variety + variety:time`)
#' and returns the term labels as a character vector.
#'
#' @param random_formula A one-sided formula for random effects.
#' @return Character vector of term names.
#' @noRd
.extract_random_term_names <- function(random_formula) {
  if (is.null(random_formula)) return(character(0L))

  if (is.character(random_formula)) {
    random_formula <- stats::as.formula(random_formula)
  }

  # Use terms() to parse the formula; labels() gives the term names
  tryCatch({
    tt <- stats::terms(random_formula)
    labels(tt)
  }, error = function(e) {
    # Fallback: simple text parsing for complex ASReml-style formulas
    # that terms() may not handle (e.g., vm(), str())
    ftext <- deparse(random_formula, width.cutoff = 500L)
    # Remove the leading "~"
    ftext <- sub("^\\s*~\\s*", "", ftext)
    # Split on "+"
    parts <- trimws(strsplit(ftext, "\\+")[[1L]])
    parts[nzchar(parts)]
  })
}


#' Find spline-related random effect terms in a fitted model
#'
#' Identifies which random terms correspond to B-spline coefficients by
#' checking for known naming conventions (spline_coef, basis, bspline, etc.).
#'
#' @param raw_result Named list from the backend.
#' @param engine Character: `"asreml"` or `"bayesreml"`.
#' @return Character vector of spline-related term names.
#' @noRd
.find_spline_terms <- function(raw_result, engine) {
  model <- raw_result[["model"]]

  # Get all random term names from the model
  all_terms <- tryCatch({
    if (engine == "asreml" && .has_asreml()) {
      names(model$G.param)
    } else if (engine == "bayesreml" && .has_bayesreml()) {
      names(model$ranef)
    } else if (engine == "lme4" && .has_lme4()) {
      names(lme4::ranef(model))
    } else if (engine == "mgcv" && .has_mgcv()) {
      # mgcv stores smooth terms; extract labels
      vapply(model$smooth, function(s) s$label, character(1L))
    } else {
      character(0L)
    }
  }, error = function(e) character(0L))

  if (length(all_terms) == 0L) return(character(0L))

  # Filter for spline-related terms using naming conventions
  # Includes Zrange_ (v0.2.0 decomposed basis), Zspline_ (MET single-stage),
  # and legacy Bsp_ / spline_coef patterns for backward compatibility.
  spline_patterns <- c("Zrange_", "Zspline_", "spline", "basis", "bspline",
                        "Bsp", "bs\\(", "str\\(", "spline_coef")
  pattern <- paste(spline_patterns, collapse = "|")
  spline_terms <- all_terms[grepl(pattern, all_terms, ignore.case = TRUE)]

  spline_terms
}


#' Extract fixed effect spline coefficients from fitted model
#'
#' Retrieves the estimated B-spline coefficients from the fixed effects of
#' the model. These are the beta_k coefficients for the coefficient function
#' beta(t) = sum_k beta_k * B_k(t).
#'
#' @param raw_result Named list from the backend.
#' @param engine Character: `"asreml"` or `"bayesreml"`.
#' @param n_basis Integer: expected number of basis functions. Used to identify
#'   which fixed effects correspond to spline coefficients.
#'
#' @return Named list with elements `estimate`, `se`, and optionally
#'   `ci_lower`, `ci_upper`. Returns NULL if no spline coefficients found.
#' @noRd
.extract_fixed_spline_coefs <- function(raw_result, engine, n_basis) {
  model <- raw_result[["model"]]

  if (engine == "asreml" && .has_asreml()) {
    # ASReml: fixed effect coefficients via coef() or summary
    coefs <- tryCatch(
      summary(model, coef = TRUE)$coef.fixed,
      error = function(e) NULL
    )

    if (is.null(coefs)) return(NULL)

    # Identify spline-related fixed effects by name pattern or count
    coef_names <- rownames(coefs)
    spline_idx <- grep("spline|basis|bs\\(", coef_names, ignore.case = TRUE)

    # Fallback: if no pattern match, check if last n_basis coefficients
    # correspond to the spline (common convention)
    if (length(spline_idx) == 0L && nrow(coefs) >= n_basis) {
      # Heuristic: assume last n_basis fixed effects are spline coefficients
      spline_idx <- seq(nrow(coefs) - n_basis + 1L, nrow(coefs))
    }

    if (length(spline_idx) != n_basis) return(NULL)

    return(list(
      estimate = coefs[spline_idx, 1L],
      se       = coefs[spline_idx, 2L]
    ))
  }

  if (engine == "bayesreml" && .has_bayesreml()) {
    # bayesreml: extract from summary
    fit_summary <- raw_result[["summary"]]
    if (is.null(fit_summary)) return(NULL)

    if (is.data.frame(fit_summary)) {
      param_names <- rownames(fit_summary)
      if (is.null(param_names)) param_names <- fit_summary[["parameter"]]

      spline_idx <- grep("spline|basis|bs\\(", param_names,
                         ignore.case = TRUE)

      # Filter to fixed effects only
      fixed_idx <- grep("^beta|^fixed", param_names, ignore.case = TRUE)
      spline_fixed_idx <- intersect(spline_idx, fixed_idx)

      if (length(spline_fixed_idx) == 0L && length(fixed_idx) >= n_basis) {
        spline_fixed_idx <- tail(fixed_idx, n_basis)
      }

      if (length(spline_fixed_idx) != n_basis) return(NULL)

      # Extract posterior summaries
      result <- list(estimate = fit_summary[spline_fixed_idx, "mean"])

      if ("sd" %in% names(fit_summary)) {
        result[["se"]] <- fit_summary[spline_fixed_idx, "sd"]
      }
      ci_cols <- intersect(names(fit_summary), c("2.5%", "q025", "ci_lower"))
      if (length(ci_cols) > 0L) {
        result[["ci_lower"]] <- fit_summary[spline_fixed_idx, ci_cols[1L]]
      }
      ci_cols_u <- intersect(names(fit_summary), c("97.5%", "q975", "ci_upper"))
      if (length(ci_cols_u) > 0L) {
        result[["ci_upper"]] <- fit_summary[spline_fixed_idx, ci_cols_u[1L]]
      }

      return(result)
    }
  }

  NULL
}
