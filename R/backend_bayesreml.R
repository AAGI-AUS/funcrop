# funcrop: bayesreml/greta backend
#
# All functions in this file are internal (not exported). They implement the
# interface between funcrop's model specification and bayesreml's Bayesian MCMC
# estimation engine (built on greta/TensorFlow).
#
# bayesreml uses ASReml-style formula syntax, so much of the formula translation
# is pass-through. The key difference is in result extraction: bayesreml returns
# full posterior distributions rather than point estimates.
#
# IMPORTANT: Every call to bayesreml:: functions is guarded by
# requireNamespace() so the package loads cleanly without bayesreml installed.

# ---- Main fitting function ---------------------------------------------------

#' Fit a model via bayesreml
#'
#' Constructs and executes a `bayesreml::bayesreml()` call from funcrop's
#' internal model specification. Handles formula construction, MCMC
#' configuration, convergence diagnostics (Rhat), and result extraction.
#'
#' @param model_spec Named list with model specification. Expected elements:
#'   \describe{
#'     \item{`fixed`}{Formula for fixed effects.}
#'     \item{`random`}{Formula for random effects.}
#'     \item{`rcov`}{Residual covariance formula (default `~ units`).}
#'     \item{`known_matrices`}{Named list of matrices (e.g., G, Z_spline).}
#'   }
#' @param data data.frame or data.table of observations.
#' @param mcmc_control Named list of MCMC parameters:
#'   \describe{
#'     \item{`n_samples`}{Total posterior draws per chain (default 2000).}
#'     \item{`warmup`}{Warmup/burn-in samples per chain (default 1000).}
#'     \item{`chains`}{Number of MCMC chains (default 4).}
#'     \item{`prior_fixed_sd`}{SD for weakly informative normal prior on
#'       fixed effects (default 10).}
#'     \item{`prior_vc_sd`}{SD for half-normal prior on variance component
#'       standard deviations (default 1).}
#'   }
#'
#' @return Named list:
#'   \describe{
#'     \item{`model`}{The fitted bayesreml object.}
#'     \item{`converged`}{Logical: TRUE if all Rhat values < 1.05.}
#'     \item{`draws`}{Posterior draws as an mcmc.list object.}
#'     \item{`summary`}{Summary data.frame with posterior statistics.}
#'   }
#' @noRd
.bayesreml_fit <- function(model_spec, data, mcmc_control = list()) {
  if (!.has_bayesreml()) {
    stop(
      "bayesreml is not installed. Cannot use the 'bayesreml' engine.\n",
      "Install via: install.packages('bayesreml')",
      call. = FALSE
    )
  }

  # Build bayesreml-compatible formulas
  formulas <- .bayesreml_build_formulas(model_spec)

  # MCMC configuration with sensible defaults
  n_samples     <- mcmc_control[["n_samples"]]     %||% 2000L
  warmup        <- mcmc_control[["warmup"]]         %||% 1000L
  chains        <- mcmc_control[["chains"]]         %||% 4L
  prior_fixed_sd <- mcmc_control[["prior_fixed_sd"]] %||% 10
  prior_vc_sd   <- mcmc_control[["prior_vc_sd"]]    %||% 1

  # Construct the bayesreml call arguments
  bayesreml_args <- list(
    fixed    = formulas[["fixed"]],
    random   = formulas[["random"]],
    rcov     = formulas[["rcov"]],
    data     = if (data.table::is.data.table(data)) data.table::copy(data) else data,
    n_samples = as.integer(n_samples),
    warmup   = as.integer(warmup),
    chains   = as.integer(chains)
  )

  # Register known matrices if supplied.
  # bayesreml uses the same vm() syntax as ASReml for known relationship
  # matrices. We attach them to the formula environment.
  known_mats <- model_spec[["known_matrices"]]
  if (length(known_mats) > 0L) {
    mat_env <- new.env(parent = environment(formulas[["fixed"]]))
    for (nm in names(known_mats)) {
      assign(nm, known_mats[[nm]], envir = mat_env)
    }
    environment(bayesreml_args[["fixed"]])  <- mat_env
    environment(bayesreml_args[["random"]]) <- mat_env
    if (!is.null(bayesreml_args[["rcov"]])) {
      environment(bayesreml_args[["rcov"]]) <- mat_env
    }
  }

  # Prior specification (bayesreml-specific arguments)
  bayesreml_args[["prior_fixed_sd"]] <- prior_fixed_sd
  bayesreml_args[["prior_vc_sd"]]    <- prior_vc_sd

  # Execute the bayesreml call with error trapping
  .msg(sprintf(
    "Running bayesreml MCMC: %d samples, %d warmup, %d chains...",
    n_samples, warmup, chains
  ))

  fit <- tryCatch(
    do.call(bayesreml::bayesreml, bayesreml_args),
    error = function(e) {
      stop(
        "bayesreml fitting failed with error:\n  ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  # Extract posterior draws and summary
  # bayesreml v0.1.0 stores draws in fit$greta$draws
  draws <- tryCatch(
    fit$greta$draws,
    error = function(e) {
      warning(
        "Failed to extract posterior draws: ", conditionMessage(e),
        call. = FALSE
      )
      NULL
    }
  )

  fit_summary <- tryCatch(
    summary(fit),
    error = function(e) {
      warning(
        "Failed to generate bayesreml summary: ", conditionMessage(e),
        call. = FALSE
      )
      NULL
    }
  )

  # Convergence diagnostic: check that all Rhat values are below threshold.
  # Rhat > 1.05 indicates poor mixing / non-convergence (Vehtari et al., 2021).
  converged <- .bayesreml_check_convergence(fit_summary, threshold = 1.05)

  if (isFALSE(converged)) {
    warning(
      "bayesreml MCMC may not have converged. ",
      "Some Rhat values exceed 1.05. Consider increasing `n_samples` ",
      "or `warmup`, or checking model specification.",
      call. = FALSE
    )
  }

  list(
    model     = fit,
    converged = converged,
    draws     = draws,
    summary   = fit_summary
  )
}


# ---- Formula construction ----------------------------------------------------

#' Build bayesreml formulas from funcrop model specification
#'
#' bayesreml uses the same formula syntax as ASReml-R (it was designed as a
#' Bayesian drop-in replacement), so formula translation is largely
#' pass-through. B-spline Z matrices are passed as known matrices via `vm()`.
#'
#' @param model_spec Named list with `fixed`, `random`, `rcov` elements.
#' @return Named list of formulas: `fixed`, `random`, `rcov`.
#' @noRd
.bayesreml_build_formulas <- function(model_spec) {
  # Fixed effects -- direct pass-through
  fixed <- model_spec[["fixed"]]
  if (is.character(fixed)) {
    fixed <- stats::as.formula(fixed)
  }

  # Random effects -- direct pass-through (bayesreml uses ASReml syntax)
  random <- model_spec[["random"]]
  if (is.character(random)) {
    random <- stats::as.formula(random)
  }

  # Residual covariance
  rcov <- model_spec[["rcov"]]
  if (is.null(rcov)) {
    rcov <- stats::as.formula("~ units")
  } else if (is.character(rcov)) {
    rcov <- stats::as.formula(rcov)
  }

  list(
    fixed  = fixed,
    random = random,
    rcov   = rcov
  )
}


# ---- BLUP extraction (posterior mean) ----------------------------------------

#' Extract posterior mean BLUPs from a bayesreml model
#'
#' Computes the posterior mean, standard deviation, and 95% credible intervals
#' for specified random effect terms from the bayesreml posterior draws.
#'
#' @param bayesreml_model A fitted bayesreml object.
#' @param terms Character vector of random terms to extract (e.g., `"variety"`,
#'   `"variety:spline_coef"`).
#'
#' @return A data.table with columns: `term`, `level`, `blup` (posterior mean),
#'   `sd` (posterior sd), `ci_lower` (2.5th percentile), `ci_upper` (97.5th
#'   percentile).
#' @noRd
.bayesreml_extract_blups <- function(bayesreml_model, terms) {
  if (!.has_bayesreml()) {
    stop("bayesreml is required for BLUP extraction.", call. = FALSE)
  }

  results <- vector("list", length(terms))

  for (i in seq_along(terms)) {
    tm <- terms[i]

    # bayesreml v0.1.0 stores BLUPs in fit$extras$blups (named list by term).
    draws_mat <- tryCatch({
      blups_list <- bayesreml_model$extras$blups
      if (!is.null(blups_list) && tm %in% names(blups_list)) {
        blups_list[[tm]]
      } else {
        NULL
      }
    }, error = function(e) {
      warning(
        sprintf("Failed to extract BLUPs for term '%s': %s", tm,
                conditionMessage(e)),
        call. = FALSE
      )
      NULL
    })

    if (is.null(draws_mat)) next

    # draws_mat: matrix with rows = posterior draws, cols = levels
    level_names <- colnames(draws_mat)
    if (is.null(level_names)) {
      level_names <- paste0(tm, "_", seq_len(ncol(draws_mat)))
    }

    # Compute posterior summaries for each level
    post_mean  <- colMeans(draws_mat)
    post_sd    <- apply(draws_mat, 2L, stats::sd)
    post_lower <- apply(draws_mat, 2L, stats::quantile, probs = 0.025)
    post_upper <- apply(draws_mat, 2L, stats::quantile, probs = 0.975)

    results[[i]] <- data.table::data.table(
      term     = tm,
      level    = level_names,
      blup     = post_mean,
      sd       = post_sd,
      ci_lower = post_lower,
      ci_upper = post_upper
    )
  }

  data.table::rbindlist(results, use.names = TRUE)
}


# ---- Variance component extraction -------------------------------------------

#' Extract variance components from a bayesreml model
#'
#' Returns posterior mean, standard deviation, and 95% credible intervals for
#' all variance components in the model.
#'
#' @param bayesreml_model A fitted bayesreml object.
#' @return A data.table with columns: `component`, `estimate` (posterior mean),
#'   `sd`, `ci_lower` (2.5%), `ci_upper` (97.5%).
#' @noRd
.bayesreml_extract_vc <- function(bayesreml_model) {
  if (!.has_bayesreml()) {
    stop("bayesreml is required for variance component extraction.",
         call. = FALSE)
  }

  # bayesreml v0.1.0 stores variance components in fit$extras$variance_comps
  vc_summary <- tryCatch(
    bayesreml_model$extras$variance_comps,
    error = function(e) {
      stop(
        "Failed to extract variance components: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  # vc_summary is expected to be a data.frame or matrix with rows = components
  if (is.data.frame(vc_summary)) {
    dt <- data.table::as.data.table(vc_summary)
    # Standardise column names
    expected_cols <- c("component", "estimate", "sd", "ci_lower", "ci_upper")
    # Attempt name mapping if columns differ
    if (!"component" %in% names(dt) && !is.null(rownames(vc_summary))) {
      dt[, component := rownames(vc_summary)]
    }
    return(dt)
  }

  # Fallback: extract from posterior draws manually
  # bayesreml v0.1.0: draws available via fit$greta$draws
  vc_draws <- tryCatch({
    all_draws <- bayesreml_model$greta$draws
    if (is.null(all_draws)) {
      stop("No posterior draws available in bayesreml fit object.")
    }
    # Extract variance-related parameters from draws
    if (inherits(all_draws, "mcmc.list")) {
      all_names <- colnames(all_draws[[1L]])
    } else if (is.matrix(all_draws)) {
      all_names <- colnames(all_draws)
    } else {
      stop("Unexpected draws format.")
    }
    vc_idx <- grep("^sigma|^tau|^variance", all_names)
    if (length(vc_idx) == 0L) {
      stop("No variance component parameters found in draws.")
    }
    if (is.matrix(all_draws)) {
      all_draws[, vc_idx, drop = FALSE]
    } else {
      # Combine chains for mcmc.list
      do.call(rbind, lapply(all_draws, function(ch) ch[, vc_idx, drop = FALSE]))
    }
  }, error = function(e) {
    stop(
      "Failed to extract VC draws: ", conditionMessage(e),
      call. = FALSE
    )
  })

  # vc_draws: matrix (n_draws x n_components)
  comp_names <- colnames(vc_draws)
  if (is.null(comp_names)) {
    comp_names <- paste0("vc_", seq_len(ncol(vc_draws)))
  }

  data.table::data.table(
    component = comp_names,
    estimate  = colMeans(vc_draws),
    sd        = apply(vc_draws, 2L, stats::sd),
    ci_lower  = apply(vc_draws, 2L, stats::quantile, probs = 0.025),
    ci_upper  = apply(vc_draws, 2L, stats::quantile, probs = 0.975)
  )
}


# ---- Full posterior extraction -----------------------------------------------

#' Extract full posterior draws for specific parameters
#'
#' Returns the complete MCMC posterior draws for a specified parameter or set
#' of parameters. Useful for custom posterior analysis, posterior predictive
#' checks, or computing derived quantities.
#'
#' @param bayesreml_model A fitted bayesreml object.
#' @param parameter Character string or vector specifying which parameters to
#'   extract. Use `"fixed"` for all fixed effects, `"random"` for all random
#'   effects, `"vc"` for variance components, or a specific parameter name.
#'
#' @return A matrix (n_draws x n_params) of posterior draws, or an mcmc.list
#'   if multiple chains are preserved.
#' @noRd
.bayesreml_extract_posterior <- function(bayesreml_model, parameter) {
  if (!.has_bayesreml()) {
    stop("bayesreml is required for posterior extraction.", call. = FALSE)
  }

  # bayesreml v0.1.0 stores draws in fit$greta$draws
  draws <- tryCatch(
    bayesreml_model$greta$draws,
    error = function(e) {
      stop(
        "Failed to extract posterior draws: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  # `draws` is expected to be an mcmc.list or matrix.
  # Filter to requested parameters.
  if (is.character(parameter) && length(parameter) == 1L) {
    # Special keywords
    if (parameter == "fixed") {
      param_regex <- "^beta\\[|^fixed\\["
    } else if (parameter == "random") {
      param_regex <- "^u\\[|^ranef\\["
    } else if (parameter == "vc") {
      param_regex <- "^sigma|^tau|^variance"
    } else {
      # Treat as a literal parameter name or regex
      param_regex <- parameter
    }

    # Extract matching columns from draws
    if (inherits(draws, "mcmc.list")) {
      all_names <- colnames(draws[[1L]])
      idx <- grep(param_regex, all_names)
      if (length(idx) == 0L) {
        stop(
          sprintf("No parameters matching '%s' found in posterior draws.",
                  parameter),
          call. = FALSE
        )
      }
      # Return as mcmc.list (preserving chain structure)
      return(lapply(draws, function(ch) ch[, idx, drop = FALSE]))
    }

    if (is.matrix(draws)) {
      all_names <- colnames(draws)
      idx <- grep(param_regex, all_names)
      if (length(idx) == 0L) {
        stop(
          sprintf("No parameters matching '%s' found in posterior draws.",
                  parameter),
          call. = FALSE
        )
      }
      return(draws[, idx, drop = FALSE])
    }
  }

  # Multiple parameter names: select columns directly
  if (is.character(parameter) && length(parameter) > 1L) {
    if (is.matrix(draws)) {
      missing_params <- setdiff(parameter, colnames(draws))
      if (length(missing_params) > 0L) {
        stop(
          sprintf("Parameters not found in posterior draws: %s",
                  paste(missing_params, collapse = ", ")),
          call. = FALSE
        )
      }
      return(draws[, parameter, drop = FALSE])
    }
  }

  # Fallback: return all draws
  draws
}


# ---- B-spline structure for bayesreml ----------------------------------------

#' Build bayesreml known_matrices and formula for B-spline random effects
#'
#' Constructs the Z matrix and associated formula components for fitting
#' variety-specific B-spline random curves in bayesreml. The approach mirrors
#' the ASReml backend but uses bayesreml's matrix registration mechanism.
#'
#' @section Statistical background:
#' Same Kronecker product structure as `.asreml_build_spline_str()`:
#' \deqn{Z = I_V \otimes Z_{spline}}
#' bayesreml handles the penalty structure through its Bayesian prior
#' specification rather than ASReml's `str()` interface. The penalty matrix
#' is incorporated as a precision matrix prior on the spline coefficients.
#'
#' @param basis A funcrop basis object or numeric matrix of B-spline
#'   evaluations (n_time x K).
#' @param variety_levels Character vector of variety factor levels.
#' @param var_structure Character: `"penalised"`, `"iid"`, or `"us"`.
#'   Default `"penalised"`.
#' @param penalty_matrix Optional penalty matrix (K x K). Required for
#'   `"penalised"` structure if not available from `basis`.
#'
#' @return Named list:
#'   \describe{
#'     \item{`Z`}{Sparse design matrix (n_obs x VK).}
#'     \item{`formula_term`}{Character string for the random formula.}
#'     \item{`known_matrices`}{Named list of matrices to register (Z matrix
#'       and optionally the penalty precision matrix).}
#'   }
#' @noRd
.bayesreml_build_spline_str <- function(basis,
                                        variety_levels,
                                        var_structure = c("penalised", "iid",
                                                          "us"),
                                        penalty_matrix = NULL) {
  var_structure <- match.arg(var_structure)

  # Extract B-spline evaluation matrix
  if (is.list(basis) && !is.null(basis[["Bmatrix"]])) {
    Z_spline <- basis[["Bmatrix"]]
    if (is.null(penalty_matrix) && !is.null(basis[["penalty"]])) {
      penalty_matrix <- basis[["penalty"]]
    }
  } else if (is.matrix(basis)) {
    Z_spline <- basis
  } else {
    stop(
      "`basis` must be a funcrop basis object or a numeric matrix.",
      call. = FALSE
    )
  }

  n_time  <- nrow(Z_spline)
  n_basis <- ncol(Z_spline)
  n_var   <- length(variety_levels)

  if (n_var < 1L) {
    stop("`variety_levels` must have at least one element.", call. = FALSE)
  }

  # Construct combined design matrix: Z = I_variety (x) Z_spline
  I_var <- Matrix::Diagonal(n_var)
  Z_combined <- Matrix::kronecker(I_var, Z_spline)

  # Prepare known matrices for bayesreml registration
  known_mats <- list()
  known_mats[["Z_spline_combined"]] <- as.matrix(Z_combined)

  # For the penalised structure, bayesreml uses the penalty matrix as a
  # precision prior on the spline coefficients:
  #   alpha ~ N(0, sigma^2_s * P^{-1})
  # which is equivalent to specifying precision = (1/sigma^2_s) * P.
  if (var_structure == "penalised") {
    if (is.null(penalty_matrix)) {
      stop(
        "Penalty matrix required for 'penalised' variance structure.\n",
        "Supply via `basis$penalty` or `penalty_matrix` argument.",
        call. = FALSE
      )
    }
    known_mats[["P_spline"]] <- as.matrix(penalty_matrix)

    # For variety-specific penalised curves, the full precision structure is:
    # Prec(alpha) = (1/sigma^2_s) * (I_V (x) P)
    P_full <- Matrix::kronecker(I_var, penalty_matrix)
    known_mats[["P_spline_full"]] <- as.matrix(P_full)
  }

  # bayesreml formula term: use vm() to reference the combined Z matrix
  # The exact syntax depends on bayesreml's API, but follows ASReml convention.
  formula_term <- switch(var_structure,
    iid = "vm(spline_coef, source = Z_spline_combined)",
    penalised = paste0(
      "vm(spline_coef, source = Z_spline_combined, ",
      "precision = P_spline_full)"
    ),
    us = paste0(
      "us(variety):vm(spline_coef, source = Z_spline_combined)"
    )
  )

  list(
    Z              = Z_combined,
    formula_term   = formula_term,
    var_structure  = var_structure,
    known_matrices = known_mats,
    n_basis        = n_basis,
    n_varieties    = n_var
  )
}


# ---- Convergence diagnostics -------------------------------------------------

#' Check MCMC convergence via Rhat
#'
#' Examines the Rhat (potential scale reduction factor) values from the
#' bayesreml summary. Convergence criterion: all Rhat < threshold (default
#' 1.05, following Vehtari et al., 2021, Bayesian Analysis).
#'
#' @param fit_summary Summary object from bayesreml (data.frame or list).
#' @param threshold Numeric: Rhat threshold for convergence (default 1.05).
#' @return Logical scalar: TRUE if converged (all Rhat below threshold).
#' @noRd
.bayesreml_check_convergence <- function(fit_summary, threshold = 1.05) {
  if (is.null(fit_summary)) {
    return(FALSE)
  }

  # Extract Rhat values -- bayesreml may store these under different names
  rhat_vals <- NULL

  if (is.data.frame(fit_summary)) {
    # Try common column names
    rhat_col <- intersect(
      names(fit_summary),
      c("Rhat", "rhat", "R_hat", "psrf")
    )
    if (length(rhat_col) > 0L) {
      rhat_vals <- fit_summary[[rhat_col[1L]]]
    }
  } else if (is.list(fit_summary) && !is.null(fit_summary[["Rhat"]])) {
    rhat_vals <- fit_summary[["Rhat"]]
  }

  if (is.null(rhat_vals) || length(rhat_vals) == 0L) {
    # Cannot assess convergence without Rhat values
    warning(
      "Could not extract Rhat values from bayesreml summary. ",
      "Convergence status unknown.",
      call. = FALSE
    )
    return(NA)
  }

  # Remove any NA Rhat values (some parameters may not have valid diagnostics)
  rhat_vals <- rhat_vals[!is.na(rhat_vals)]

  if (length(rhat_vals) == 0L) {
    return(NA)
  }

  all(rhat_vals < threshold)
}
