# funcrop: ASReml-R v4.2 backend
#
# All functions in this file are internal (not exported). They implement the
# interface between funcrop's model specification and ASReml-R's formula-based
# API for REML estimation of linear mixed models with B-spline random effects.
#
# IMPORTANT: Every call to asreml:: functions is guarded by requireNamespace()
# so the package loads cleanly even without ASReml installed.

# ---- Main fitting function ---------------------------------------------------

#' Fit a model via ASReml-R
#'
#' Constructs and executes an `asreml::asreml()` call from funcrop's internal
#' model specification. Handles formula construction, relationship matrix
#' registration, convergence checking, and result extraction.
#'
#' @param model_spec Named list with model specification. Expected elements:
#'   \describe{
#'     \item{`fixed`}{Formula for fixed effects.}
#'     \item{`random`}{Formula for random effects.}
#'     \item{`rcov`}{Residual covariance formula (default `~ units`).}
#'     \item{`known_matrices`}{Named list of matrices (e.g., G, Z_spline).}
#'     \item{`start_values`}{Optional named vector of starting values.}
#'     \item{`maxiter`}{Maximum REML iterations (default 50).}
#'     \item{`workspace`}{ASReml workspace in bytes (default 128e6).}
#'   }
#' @param data data.frame or data.table of observations.
#' @param control Named list of additional ASReml control parameters passed to
#'   `asreml::asreml.options()`.
#'
#' @return Named list:
#'   \describe{
#'     \item{`model`}{The fitted asreml object.}
#'     \item{`converged`}{Logical: did REML converge?}
#'     \item{`log_lik`}{Log-likelihood at convergence.}
#'     \item{`n_iter`}{Number of REML iterations.}
#'   }
#' @noRd
.asreml_fit <- function(model_spec, data, control = list()) {
  if (!.has_asreml()) {
    stop(
      "ASReml-R is not installed. Cannot use the 'asreml' engine.\n",
      "See https://vsni.co.uk/software/asreml-r for installation.",
      call. = FALSE
    )
  }

  # Build ASReml-compatible formulas from the model specification
  formulas <- .asreml_build_formulas(model_spec)

  # Prepare arguments for asreml() call
  maxiter   <- model_spec[["maxiter"]]   %||% 50L
  workspace <- model_spec[["workspace"]] %||% 128e6

  asreml_args <- list(
    fixed     = formulas[["fixed"]],
    random    = formulas[["random"]],
    residual  = formulas[["rcov"]],
    data      = data,
    maxiter   = as.integer(maxiter),
    workspace = workspace,
    trace     = isTRUE(getOption("funcrop.verbose", default = FALSE))
  )

  # Register known matrices (relationship / design matrices) if supplied.
  # ASReml-R v4.2 uses vm(term, source = mat) -- matrices must be assigned

  # into the data environment or registered via asreml's attr mechanism.
  known_mats <- model_spec[["known_matrices"]]
  if (length(known_mats) > 0L) {
    # Assign matrices into a local environment attached to the data
    mat_env <- new.env(parent = environment(formulas[["fixed"]]))
    for (nm in names(known_mats)) {
      assign(nm, known_mats[[nm]], envir = mat_env)
    }
    # Ensure formulas use this environment
    environment(asreml_args[["fixed"]])    <- mat_env
    environment(asreml_args[["random"]])   <- mat_env
    if (!is.null(asreml_args[["residual"]])) {
      environment(asreml_args[["residual"]]) <- mat_env
    }
  }

  # Starting values for variance components (optional)
  if (!is.null(model_spec[["start_values"]])) {
    asreml_args[["start.values"]] <- TRUE
  }

  # Apply any user-level control overrides
  if (length(control) > 0L) {
    asreml_args <- utils::modifyList(asreml_args, control)
  }

  # Execute the ASReml call with error trapping
  fit <- tryCatch(
    do.call(asreml::asreml, asreml_args),
    error = function(e) {
      stop(
        "ASReml fitting failed with error:\n  ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  # If start.values were requested, re-fit with user-supplied values
  if (!is.null(model_spec[["start_values"]])) {
    sv <- fit  # start.values = TRUE returns initial value table, not a model
    sv_tab <- sv$vparameters.table
    # Overwrite with user-supplied values where names match
    user_sv <- model_spec[["start_values"]]
    matched <- intersect(names(user_sv), sv_tab$Component)
    if (length(matched) > 0L) {
      idx <- match(matched, sv_tab$Component)
      sv_tab$Value[idx] <- user_sv[matched]
    }
    asreml_args[["start.values"]] <- FALSE
    asreml_args[["G.param"]]      <- sv_tab
    asreml_args[["R.param"]]      <- sv_tab

    fit <- tryCatch(
      do.call(asreml::asreml, asreml_args),
      error = function(e) {
        stop(
          "ASReml fitting (with start values) failed:\n  ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }

  # Extract convergence and summary information
  converged <- fit$converge
  log_lik   <- fit$loglik
  n_iter    <- fit$nedf  # number of iterations is in the trace

  # Attempt to get actual iteration count from the model
  n_iter <- tryCatch(
    length(fit$monitor$loglik),
    error = function(e) NA_integer_
  )

  list(
    model     = fit,
    converged = isTRUE(converged),
    log_lik   = log_lik,
    n_iter    = n_iter
  )
}


# ---- Formula construction ----------------------------------------------------

#' Build ASReml formulas from funcrop model specification
#'
#' Translates the internal model_spec into ASReml-compatible formulas.
#' Handles special terms: `vm()` for relationship matrices, `at()` for
#' heterogeneous variances, `us()`/`fa()` for structured covariance,
#' `ar1()` for spatial correlation, and `str()` for custom B-spline
#' structures.
#'
#' @param model_spec Named list with `fixed`, `random`, `rcov` elements.
#'   These may be formulas or character strings.
#' @return Named list of formulas: `fixed`, `random`, `rcov`.
#' @noRd
.asreml_build_formulas <- function(model_spec) {
  # Fixed effects formula -- pass through (ASReml uses standard R formula syntax)
  fixed <- model_spec[["fixed"]]
  if (is.character(fixed)) {
    fixed <- stats::as.formula(fixed)
  }

  # Random effects formula
  random <- model_spec[["random"]]
  if (is.character(random)) {
    random <- stats::as.formula(random)
  }

  # Residual covariance formula (default: ~ units for iid residuals)
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


# ---- BLUP extraction --------------------------------------------------------

#' Extract BLUPs from an ASReml model
#'
#' Uses `predict()` from ASReml to obtain best linear unbiased predictions
#' (BLUPs) for specified random terms.
#'
#' @param asreml_model A fitted asreml object.
#' @param terms Character vector of random terms to predict (e.g., `"variety"`,
#'   `"variety:spline_basis"`).
#'
#' @return A data.table with columns: `term`, `level`, `blup`, `se`.
#' @noRd
.asreml_extract_blups <- function(asreml_model, terms) {
  if (!.has_asreml()) {
    stop("ASReml-R is required for BLUP extraction.", call. = FALSE)
  }

  results <- vector("list", length(terms))

  for (i in seq_along(terms)) {
    tm <- terms[i]

    pred <- tryCatch(
      asreml::predict.asreml(
        asreml_model,
        classify = tm,
        present  = all.vars(asreml_model$call$random)
      ),
      error = function(e) {
        warning(
          sprintf("Failed to extract BLUPs for term '%s': %s", tm,
                  conditionMessage(e)),
          call. = FALSE
        )
        NULL
      }
    )

    if (is.null(pred)) next

    # ASReml predict returns a list with $pvals (predictions data.frame)
    pvals <- pred$pvals

    # Construct data.table from predictions
    # The classify columns form the 'level' identifier
    classify_cols <- setdiff(names(pvals), c("predicted.value", "std.error",
                                              "status"))
    if (length(classify_cols) == 1L) {
      level_vec <- as.character(pvals[[classify_cols]])
    } else {
      # Multiple classify columns -- paste together
      level_vec <- do.call(
        paste, c(lapply(classify_cols, function(cc) pvals[[cc]]), sep = ":")
      )
    }

    results[[i]] <- data.table::data.table(
      term  = tm,
      level = level_vec,
      blup  = pvals[["predicted.value"]],
      se    = pvals[["std.error"]]
    )
  }

  data.table::rbindlist(results, use.names = TRUE)
}


# ---- Variance component extraction -------------------------------------------

#' Extract variance components from an ASReml model
#'
#' Retrieves the variance component estimates, standard errors, z-ratios,
#' and boundary status from an ASReml model summary.
#'
#' @param asreml_model A fitted asreml object.
#' @return A data.table with columns: `component`, `estimate`, `se`,
#'   `z_ratio`, `bound`.
#' @noRd
.asreml_extract_vc <- function(asreml_model) {
  if (!.has_asreml()) {
    stop("ASReml-R is required for variance component extraction.", call. = FALSE)
  }

  # summary.asreml() returns a list with $varcomp data.frame
  vc_summary <- tryCatch(
    summary(asreml_model)$varcomp,
    error = function(e) {
      stop(
        "Failed to extract variance components: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  data.table::data.table(
    component = rownames(vc_summary),
    estimate  = vc_summary[["component"]],
    se        = vc_summary[["std.error"]],
    z_ratio   = vc_summary[["z.ratio"]],
    bound     = vc_summary[["bound"]]
  )
}


# ---- B-spline structure for ASReml -------------------------------------------

#' Build ASReml str() structure for B-spline random effects
#'
#' Constructs the design matrix Z and variance structure needed to fit
#' variety-specific B-spline random curves in ASReml. This is the core function
#' bridging funcrop's B-spline basis system with ASReml's `str()` interface.
#'
#' @section Statistical background:
#' For \eqn{V} varieties each with a smooth random curve over time, we model:
#' \deqn{u_v(t) = \sum_{k=1}^{K} \alpha_{vk} B_k(t)}
#' where \eqn{B_k(t)} are the B-spline basis functions and \eqn{\alpha_{vk}}
#' are variety-specific random spline coefficients.
#'
#' The combined design matrix is the Kronecker product:
#' \deqn{Z = I_V \otimes Z_{spline}}
#' where \eqn{Z_{spline}} is the \eqn{n_{time} \times K} B-spline evaluation
#' matrix from [make_Zspline()] (or equivalent basis evaluation), and
#' \eqn{I_V} is the \eqn{V \times V} identity matrix selecting the variety.
#'
#' The associated variance structure for the spline coefficients can be:
#' - `"iid"`: \eqn{Var(\alpha) = \sigma^2_s I_{VK}} -- independent coefficients
#' - `"penalised"`: \eqn{Var(\alpha) = \sigma^2_s (I_V \otimes P^{-})} where
#'   \eqn{P} is the P-spline penalty matrix
#' - `"us"`: \eqn{Var(\alpha) = G_s \otimes I_K} -- unstructured covariance
#'   among varieties for each spline coefficient
#'
#' @param basis A funcrop basis object (containing `Bmatrix`, `penalty`, and
#'   `knots`), or a numeric matrix of B-spline evaluations (n_time x K).
#' @param variety_levels Character vector of variety factor levels. Used to
#'   construct the Kronecker product design matrix.
#' @param var_structure Character string: `"iid"`, `"penalised"`, or `"us"`.
#'   Default `"penalised"`.
#' @param penalty_matrix Optional penalty matrix (K x K). If NULL and
#'   var_structure is `"penalised"`, attempts to extract from `basis`.
#'
#' @return Named list suitable for ASReml's `str()` interface:
#'   \describe{
#'     \item{`Z`}{Sparse design matrix (n_obs x VK).}
#'     \item{`formula_term`}{Character string for inclusion in the random
#'       formula.}
#'     \item{`variance_model`}{List describing the variance structure for
#'       ASReml's `str()`.}
#'     \item{`known_matrices`}{Named list of matrices that need to be
#'       registered (e.g., the inverse penalty matrix for `vm()`).}
#'   }
#' @noRd
.asreml_build_spline_str <- function(basis,
                                     variety_levels,
                                     var_structure = c("penalised", "iid", "us"),
                                     penalty_matrix = NULL) {
  var_structure <- match.arg(var_structure)

  # Extract the B-spline evaluation matrix
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

  n_time <- nrow(Z_spline)
  n_basis <- ncol(Z_spline)
  n_var   <- length(variety_levels)

  if (n_var < 1L) {
    stop("`variety_levels` must have at least one element.", call. = FALSE)
  }

  # Construct the combined design matrix: Z = I_variety (x) Z_spline
  # This is block-diagonal: each variety gets its own copy of Z_spline.
  # Using sparse representation for efficiency.
  I_var <- Matrix::Diagonal(n_var)
  Z_combined <- Matrix::kronecker(I_var, Z_spline)

  # Build the variance structure specification
  known_mats <- list()
  formula_term <- NULL

  switch(var_structure,
    iid = {
      # Simple iid variance on all VK coefficients
      # In ASReml: str(~ variety:basis_coef, ~id(variety):id(basis_coef))
      formula_term <- sprintf(
        "str(~ variety:spline_coef, ~ id(variety):id(spline_coef))"
      )
    },
    penalised = {
      # P-spline penalty: Var(alpha_v) = sigma^2_s * P^{-}
      # where P^{-} is the generalised inverse of the penalty matrix.
      # In ASReml, pass the penalty inverse as a known matrix via vm().
      if (is.null(penalty_matrix)) {
        stop(
          "Penalty matrix required for 'penalised' variance structure.\n",
          "Supply via `basis$penalty` or `penalty_matrix` argument.",
          call. = FALSE
        )
      }
      # Compute generalised inverse of penalty matrix for use as vm() source.
      # Add small ridge for numerical stability if penalty is singular
      # (which it is for difference penalties: rank = K - order).
      P_inv <- tryCatch(
        solve(penalty_matrix + Matrix::Diagonal(n_basis) * 1e-8),
        error = function(e) {
          # Fallback: Moore-Penrose pseudoinverse via SVD
          svd_P <- svd(as.matrix(penalty_matrix))
          tol <- max(dim(penalty_matrix)) * max(svd_P$d) * .Machine$double.eps
          pos <- svd_P$d > tol
          if (sum(pos) == 0L) {
            return(Matrix::Diagonal(n_basis))
          }
          Matrix::Matrix(
            svd_P$v[, pos, drop = FALSE] %*%
              diag(1 / svd_P$d[pos], nrow = sum(pos)) %*%
              t(svd_P$u[, pos, drop = FALSE])
          )
        }
      )

      known_mats[["P_inv_spline"]] <- as.matrix(P_inv)

      # ASReml formula term: variety crossed with spline coefficients,
      # spline coefficients penalised via vm() with the penalty inverse
      formula_term <- paste0(
        "str(~ variety:spline_coef, ",
        "~ id(variety):vm(spline_coef, source = P_inv_spline))"
      )
    },
    us = {
      # Unstructured variety covariance for each spline coefficient
      # Var(alpha) = G_var (x) I_K
      # In ASReml: str(~ variety:spline_coef, ~ us(variety):id(spline_coef))
      formula_term <- sprintf(
        "str(~ variety:spline_coef, ~ us(variety):id(spline_coef))"
      )
    }
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

# Note: %||% (null-coalescing) is available from base R (>= 4.0.0).
# funcrop requires R >= 4.1.0, so no custom definition is needed.
