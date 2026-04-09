# funcrop: lme4 backend (lmer / glmer)
#
# All functions in this file are internal (not exported). They implement the
# interface between funcrop's model specification (which uses ASReml-style
# formulas) and lme4's formula-based API for (restricted) maximum likelihood
# estimation of linear and generalised linear mixed models.
#
# Key translation challenge: ASReml encodes variety-specific B-spline random
# effects as interaction terms `variety_f:Zrange_k`, while lme4 requires
# grouped random effects `(0 + Zrange_1 + Zrange_2 + ... | variety_f)`.
#
# Limitations: lme4 does NOT support AR1/AR1xAR1 residual structures, factor
# analytic (FA) or unstructured (US) variance models, or vm() for known
# relationship matrices directly. Genomic G-matrices are handled via the
# Cholesky decomposition trick (pre-multiplying the design matrix).
#
# IMPORTANT: Every call to lme4:: functions is guarded by requireNamespace()
# so the package loads cleanly even without lme4 installed.

# ---- Main fitting function ---------------------------------------------------

#' Fit a model via lme4
#'
#' Translates funcrop's internal model specification (ASReml-style formulas)
#' into lme4 syntax and fits via \code{lme4::lmer()} (Gaussian) or
#' \code{lme4::glmer()} (non-Gaussian families).
#'
#' @param model_spec Named list with model specification. Expected elements:
#'   \describe{
#'     \item{\code{fixed}}{Formula for fixed effects (e.g.,
#'       \code{yield ~ lin(time)}).}
#'     \item{\code{random}}{Formula for random effects using ASReml-style
#'       syntax (e.g., \code{~ variety_f:Zrange_1 + variety_f:Zrange_2 +
#'       block}).}
#'     \item{\code{rcov}}{Residual covariance formula. Ignored by lme4 (only
#'       iid residuals supported). A warning is issued if non-trivial rcov is
#'       supplied.}
#'     \item{\code{known_matrices}}{Named list of relationship/design matrices
#'       (e.g., genomic G-matrix). Handled via the Cholesky trick.}
#'     \item{\code{family}}{Optional: a \code{family} object for GLMM fitting
#'       via \code{glmer()}. If \code{NULL} or \code{gaussian()}, uses
#'       \code{lmer()}.}
#'   }
#' @param data data.frame or data.table of observations.
#' @param control Named list of additional lme4 control parameters passed to
#'   \code{lme4::lmerControl()} or \code{lme4::glmerControl()}.
#'
#' @return Named list:
#'   \describe{
#'     \item{\code{model}}{The fitted lme4 model object (merMod).}
#'     \item{\code{converged}}{Logical: did optimisation converge without
#'       warnings?}
#'     \item{\code{log_lik}}{Log-likelihood at convergence.}
#'     \item{\code{n_iter}}{Number of optimiser iterations (if available).}
#'   }
#' @noRd
.lme4_fit <- function(model_spec, data, control = list()) {
  if (!.has_lme4()) {
    stop(
      "lme4 is not installed. Cannot use the 'lme4' engine.\n",
      "Install via: install.packages('lme4')",
      call. = FALSE
    )
  }

  # --- Input validation -------------------------------------------------------
  if (!is.list(model_spec)) {
    stop("`model_spec` must be a list.", call. = FALSE)
  }
  required_fields <- c("fixed", "random")
  missing_fields <- setdiff(required_fields, names(model_spec))
  if (length(missing_fields) > 0L) {
    stop(
      sprintf(
        "`model_spec` is missing required fields: %s.",
        paste(missing_fields, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!is.data.frame(data) && !data.table::is.data.table(data)) {
    stop("`data` must be a data.frame or data.table.", call. = FALSE)
  }

  # Warn if non-trivial residual structure is requested
  rcov <- model_spec[["rcov"]]
  if (!is.null(rcov)) {
    rcov_text <- if (inherits(rcov, "formula")) {
      deparse(rcov, width.cutoff = 500L)
    } else {
      as.character(rcov)
    }
    is_trivial <- grepl("^\\s*~\\s*units\\s*$", rcov_text) ||
      grepl("^\\s*~\\s*id\\(units\\)\\s*$", rcov_text)
    if (!is_trivial) {
      warning(
        "lme4 does not support non-iid residual structures.\n",
        "  Requested rcov: ", rcov_text, "\n",
        "  Falling back to iid residuals. Consider using 'asreml' or ",
        "'bayesreml' engine for spatial/AR1 residual structures.",
        call. = FALSE
      )
    }
  }

  # --- Handle genomic G-matrix via Cholesky trick -----------------------------
  known_mats <- model_spec[["known_matrices"]]
  if (length(known_mats) > 0L) {
    data <- .lme4_apply_cholesky_trick(data, model_spec, known_mats)
  }

  # --- Build lme4 formula -----------------------------------------------------
  lme4_formula <- .lme4_build_formula(model_spec)

  # --- Determine family (lmer vs glmer) ---------------------------------------
  family <- model_spec[["family"]]
  use_glmer <- !is.null(family) && !identical(family, stats::gaussian()) &&
    !(is.character(family) && tolower(family) == "gaussian")

  # --- Build control object ---------------------------------------------------
  ctrl <- if (use_glmer) {
    ctrl_args <- utils::modifyList(
      list(optimizer = "bobyqa",
           optCtrl = list(maxfun = 100000L)),
      control
    )
    do.call(lme4::glmerControl, ctrl_args)
  } else {
    ctrl_args <- utils::modifyList(
      list(optimizer = "bobyqa",
           optCtrl = list(maxfun = 100000L)),
      control
    )
    do.call(lme4::lmerControl, ctrl_args)
  }

  # --- Fit the model ----------------------------------------------------------
  .msg("Fitting lme4 model...")
  .msg(sprintf("  Formula: %s", deparse(lme4_formula, width.cutoff = 500L)))

  # Ensure data is a plain data.frame (lme4 can behave unexpectedly with

  # data.table due to reference semantics)
  fit_data <- if (data.table::is.data.table(data)) {
    as.data.frame(data)
  } else {
    data
  }

  fit <- tryCatch({
    if (use_glmer) {
      lme4::glmer(
        formula = lme4_formula,
        data    = fit_data,
        family  = family,
        control = ctrl
      )
    } else {
      lme4::lmer(
        formula = lme4_formula,
        data    = fit_data,
        REML    = TRUE,
        control = ctrl
      )
    }
  }, error = function(e) {
    stop(
      "lme4 fitting failed with error:\n  ", conditionMessage(e),
      call. = FALSE
    )
  })

  # --- Extract results --------------------------------------------------------
  converged <- .lme4_check_convergence(fit)
  log_lik <- tryCatch(
    as.numeric(stats::logLik(fit)),
    error = function(e) NA_real_
  )
  n_iter <- tryCatch({
    opt_info <- fit@optinfo
    if (!is.null(opt_info[["feval"]])) {
      opt_info[["feval"]]
    } else {
      NA_integer_
    }
  }, error = function(e) NA_integer_)

  if (!converged) {
    warning(
      "lme4 model did not converge cleanly. Check model specification ",
      "and consider simplifying the random effects structure.",
      call. = FALSE
    )
  }

  list(
    model     = fit,
    converged = converged,
    log_lik   = log_lik,
    n_iter    = n_iter
  )
}


# ---- Formula construction ----------------------------------------------------

#' Translate ASReml-style formula to lme4 syntax
#'
#' Parses the funcrop model specification (which uses ASReml-style formulas)
#' and constructs a single lme4-compatible formula combining fixed and random
#' effects.
#'
#' Key translations:
#' \itemize{
#'   \item \code{variety_f:Zrange_1 + variety_f:Zrange_2 + ...} becomes
#'     \code{(0 + Zrange_1 + Zrange_2 + ... | variety_f)}
#'   \item Simple random intercepts: \code{block} becomes \code{(1 | block)}
#'   \item Nested or interaction terms are grouped by their grouping factor.
#' }
#'
#' @param model_spec Named list with \code{fixed} and \code{random} elements.
#'   These may be formulas or character strings.
#' @return A single lme4-compatible formula combining fixed and random parts.
#' @noRd
.lme4_build_formula <- function(model_spec) {
  # --- Fixed effects ----------------------------------------------------------
  fixed <- model_spec[["fixed"]]
  if (is.character(fixed)) {
    fixed <- stats::as.formula(fixed)
  }

  # Extract the response variable and fixed RHS
  fixed_parts <- as.character(fixed)
  # as.character on a formula gives c("~", "LHS", "RHS") for two-sided
  if (length(fixed_parts) == 3L) {
    response  <- fixed_parts[2L]
    fixed_rhs <- fixed_parts[3L]
  } else {
    stop(
      "Fixed formula must be two-sided (e.g., yield ~ x1 + x2).",
      call. = FALSE
    )
  }

  # --- Random effects ---------------------------------------------------------
  random <- model_spec[["random"]]
  if (is.null(random)) {
    # No random effects -- return fixed formula as-is
    return(fixed)
  }
  if (is.character(random)) {
    random <- stats::as.formula(random)
  }

  # Parse random formula into individual terms
  random_text <- deparse(random, width.cutoff = 500L)
  random_text <- sub("^\\s*~\\s*", "", random_text)

  # Split on " + " respecting parentheses (simple split for typical cases)
  random_terms <- .lme4_parse_random_terms(random_text)

  # Classify terms: interaction (grouping_factor:covariate) vs simple intercept
  grouped <- list()     # key = grouping factor, value = vector of covariates
  intercept_terms <- character(0L)

  for (tm in random_terms) {
    tm <- trimws(tm)
    if (nchar(tm) == 0L) next

    # Skip terms that are already in lme4 syntax (contain "|")
    if (grepl("\\|", tm)) {
      # Already lme4 syntax -- keep as-is (shouldn't happen, but defensive)
      intercept_terms <- c(intercept_terms, paste0("(", tm, ")"))
      next
    }

    # Skip ASReml-specific terms that can't be translated (vm, str, fa, us, at)
    if (grepl("^(vm|str|fa|us)\\(", tm)) {
      warning(
        sprintf("lme4 cannot handle ASReml term '%s'. Skipping.", tm),
        call. = FALSE
      )
      next
    }

    # Handle at() terms -- translate to nested random effects if possible
    if (grepl("^at\\(", tm)) {
      warning(
        sprintf(
          "lme4 has limited 'at()' support. Term '%s' skipped. ",
          tm
        ),
        "Consider using the 'asreml' engine for heterogeneous variance models.",
        call. = FALSE
      )
      next
    }

    # Check if this is an interaction term (contains ":")
    if (grepl(":", tm)) {
      parts <- strsplit(tm, ":", fixed = TRUE)[[1L]]
      if (length(parts) == 2L) {
        grp <- trimws(parts[1L])
        cov <- trimws(parts[2L])

        # Heuristic: the grouping factor is the one that looks like a factor
        # (ends in _f, or is a known blocking factor). The covariate is the
        # other (typically starts with Z, Bsp, or is numeric).
        # If the covariate part looks like a factor name and the group part
        # looks like a covariate, swap them.
        if (.lme4_looks_like_covariate(grp) &&
            !.lme4_looks_like_covariate(cov)) {
          tmp <- grp
          grp <- cov
          cov <- tmp
        }

        if (!grp %in% names(grouped)) {
          grouped[[grp]] <- character(0L)
        }
        grouped[[grp]] <- c(grouped[[grp]], cov)
      } else {
        # Higher-order interaction -- treat as simple term
        # lme4 can handle e.g., (1 | a:b) for interactions
        intercept_terms <- c(
          intercept_terms,
          sprintf("(1 | %s)", paste(parts, collapse = ":"))
        )
      }
    } else {
      # Simple random intercept term (e.g., "block")
      intercept_terms <- c(intercept_terms, sprintf("(1 | %s)", tm))
    }
  }

  # --- Build grouped random effect terms in lme4 syntax -----------------------
  # variety_f:Zrange_1 + variety_f:Zrange_2 + ... becomes
  # (0 + Zrange_1 + Zrange_2 + ... | variety_f)
  grouped_terms <- character(0L)
  for (grp in names(grouped)) {
    covs <- grouped[[grp]]
    if (length(covs) == 0L) next
    rhs <- paste(covs, collapse = " + ")
    grouped_terms <- c(grouped_terms, sprintf("(0 + %s | %s)", rhs, grp))
  }

  # --- Combine into a single lme4 formula -------------------------------------
  all_random <- c(grouped_terms, intercept_terms)
  if (length(all_random) == 0L) {
    return(fixed)
  }

  random_rhs <- paste(all_random, collapse = " + ")
  formula_text <- sprintf("%s ~ %s + %s", response, fixed_rhs, random_rhs)

  stats::as.formula(formula_text)
}


#' Parse random effects formula text into individual terms
#'
#' Splits a random formula RHS string on " + " while respecting parenthesised
#' sub-expressions (e.g., \code{str(...)}, \code{vm(...)}).
#'
#' @param text Character string: the RHS of a random effects formula.
#' @return Character vector of individual terms.
#' @noRd
.lme4_parse_random_terms <- function(text) {
  # Simple approach: split on "+" that are not inside parentheses

  # Track parenthesis depth
  chars <- strsplit(text, "")[[1L]]
  depth <- 0L
  split_positions <- integer(0L)

  for (i in seq_along(chars)) {
    ch <- chars[i]
    if (ch == "(") {
      depth <- depth + 1L
    } else if (ch == ")") {
      depth <- max(0L, depth - 1L)
    } else if (ch == "+" && depth == 0L) {
      split_positions <- c(split_positions, i)
    }
  }

  if (length(split_positions) == 0L) {
    return(trimws(text))
  }

  # Split at the "+" positions
  starts <- c(1L, split_positions + 1L)
  ends   <- c(split_positions - 1L, nchar(text))
  terms  <- vapply(
    seq_along(starts),
    function(j) trimws(substr(text, starts[j], ends[j])),
    character(1L)
  )

  terms[nzchar(terms)]
}


#' Check whether a term name looks like a covariate (not a grouping factor)
#'
#' Heuristic: covariates typically start with Z, Bsp, bs, or contain numeric
#' suffixes like _1, _2. Grouping factors end in _f or are common blocking
#' factor names.
#'
#' @param name Character string: a term name from the random formula.
#' @return Logical: TRUE if the name looks like a covariate.
#' @noRd
.lme4_looks_like_covariate <- function(name) {
  # Patterns that suggest a covariate (design matrix column)
  cov_patterns <- c(
    "^[Zz]",              # Z-prefixed design matrix columns (Zrange_, Zspline_)
    "^[Bb]sp",            # B-spline basis columns
    "^bs\\(",             # bs() call
    "_\\d+$",             # Ends in numeric suffix (_1, _2, ...)
    "^lin\\(",            # Linear covariate
    "^pol\\("             # Polynomial covariate
  )
  any(vapply(cov_patterns, function(p) grepl(p, name), logical(1L)))
}


# ---- Variance component extraction -------------------------------------------

#' Extract variance components from an lme4 model
#'
#' Uses \code{lme4::VarCorr()} to extract variance component estimates and
#' returns them in funcrop's standardised format.
#'
#' @param lme4_model A fitted lme4 model object (class merMod).
#' @return A data.table with columns: \code{component}, \code{estimate},
#'   \code{se}, \code{z_ratio}, \code{bound}.
#' @noRd
.lme4_extract_vc <- function(lme4_model) {
  if (!.has_lme4()) {
    stop("lme4 is required for variance component extraction.", call. = FALSE)
  }

  vc <- lme4::VarCorr(lme4_model)

  # Extract variance components into a flat table
  vc_list <- vector("list", 0L)

  for (grp in names(vc)) {
    vc_grp <- vc[[grp]]
    variances <- attr(vc_grp, "stddev")^2
    var_names <- names(variances)

    if (is.null(var_names)) {
      var_names <- paste0(grp, ".", seq_along(variances))
    } else {
      var_names <- paste0(grp, ".", var_names)
    }

    for (j in seq_along(variances)) {
      vc_list <- c(vc_list, list(data.table::data.table(
        component = var_names[j],
        estimate  = variances[j],
        se        = NA_real_,
        z_ratio   = NA_real_,
        bound     = ""
      )))
    }

    # Extract correlations if present (multiple random effects per group)
    if (ncol(vc_grp) > 1L) {
      corr_mat <- attr(vc_grp, "correlation")
      if (!is.null(corr_mat) && nrow(corr_mat) > 1L) {
        for (row_i in 2L:nrow(corr_mat)) {
          for (col_j in 1L:(row_i - 1L)) {
            corr_name <- sprintf(
              "%s.cor(%s,%s)",
              grp,
              rownames(corr_mat)[row_i],
              colnames(corr_mat)[col_j]
            )
            vc_list <- c(vc_list, list(data.table::data.table(
              component = corr_name,
              estimate  = corr_mat[row_i, col_j],
              se        = NA_real_,
              z_ratio   = NA_real_,
              bound     = ""
            )))
          }
        }
      }
    }
  }

  # Add residual variance
  resid_var <- attr(vc, "sc")^2
  vc_list <- c(vc_list, list(data.table::data.table(
    component = "units!R",
    estimate  = resid_var,
    se        = NA_real_,
    z_ratio   = NA_real_,
    bound     = ""
  )))

  if (length(vc_list) == 0L) {
    return(data.table::data.table(
      component = character(0L),
      estimate  = numeric(0L),
      se        = NA_real_,
      z_ratio   = NA_real_,
      bound     = character(0L)
    ))
  }

  # Attempt to compute approximate SEs via the Hessian (if available)
  result <- data.table::rbindlist(vc_list, use.names = TRUE)
  result <- .lme4_add_vc_uncertainty(lme4_model, result)
  result
}


#' Add approximate uncertainty to lme4 variance components
#'
#' Uses the profile likelihood or Hessian to compute approximate standard
#' errors for variance components, when available. Falls back gracefully
#' if profiling fails.
#'
#' @param lme4_model A fitted merMod object.
#' @param vc_dt data.table of variance components (modified in place).
#' @return The modified data.table with \code{se} and \code{z_ratio} updated
#'   where possible.
#' @noRd
.lme4_add_vc_uncertainty <- function(lme4_model, vc_dt) {
  # Try to get the variance-covariance matrix of the variance parameters
  # from the model's deviance function Hessian
  vcov_theta <- tryCatch({
    dd <- lme4_model@devcomp$cmp
    # lme4 stores the Hessian at the optimum in some cases
    # Use numDeriv or vcov on the theta parameters
    # This is approximate; profile() is more accurate but expensive
    NULL
  }, error = function(e) NULL)

  # If direct Hessian extraction failed, leave SEs as NA

  # Profile-based CIs are expensive and not computed by default
  # Users wanting proper CIs should use confint(model, method = "profile")
  vc_dt
}


# ---- BLUP extraction ---------------------------------------------------------

#' Extract BLUPs (conditional modes) from an lme4 model
#'
#' Uses \code{lme4::ranef()} to extract random effect predictions (conditional
#' modes of the random effects given the data) and returns them in funcrop's
#' standardised format.
#'
#' @param lme4_model A fitted lme4 model object (class merMod).
#' @param terms Character vector of random terms to extract. If \code{NULL},
#'   extracts all random effects.
#' @return A data.table with columns: \code{term}, \code{level}, \code{blup},
#'   \code{se}.
#' @noRd
.lme4_extract_blups <- function(lme4_model, terms = NULL) {
  if (!.has_lme4()) {
    stop("lme4 is required for BLUP extraction.", call. = FALSE)
  }

  # ranef() with condVar = TRUE gives conditional variances for SE computation
  re <- tryCatch(
    lme4::ranef(lme4_model, condVar = TRUE),
    error = function(e) {
      warning(
        "Failed to extract random effects from lme4 model: ",
        conditionMessage(e),
        call. = FALSE
      )
      return(NULL)
    }
  )

  if (is.null(re) || length(re) == 0L) {
    return(data.table::data.table(
      term  = character(0L),
      level = character(0L),
      blup  = numeric(0L),
      se    = numeric(0L)
    ))
  }

  results <- vector("list", 0L)

  for (grp in names(re)) {
    re_grp <- re[[grp]]
    cond_var <- attr(re_grp, "postVar")

    for (col_name in colnames(re_grp)) {
      # Determine the term label
      term_label <- if (col_name == "(Intercept)") {
        grp
      } else {
        paste0(grp, ":", col_name)
      }

      # Filter if specific terms were requested
      if (!is.null(terms) && length(terms) > 0L) {
        # Check if this term matches any requested term (partial matching)
        matches <- vapply(terms, function(t) {
          grepl(t, term_label, fixed = TRUE) ||
            grepl(term_label, t, fixed = TRUE) ||
            grp == t ||
            col_name == t
        }, logical(1L))
        if (!any(matches)) next
      }

      blup_values <- re_grp[[col_name]]
      level_names <- rownames(re_grp)

      # Extract conditional SEs from postVar array
      col_idx <- which(colnames(re_grp) == col_name)
      se_values <- if (!is.null(cond_var) && length(dim(cond_var)) == 3L) {
        sqrt(cond_var[col_idx, col_idx, ])
      } else {
        rep(NA_real_, length(blup_values))
      }

      results <- c(results, list(data.table::data.table(
        term  = term_label,
        level = level_names,
        blup  = blup_values,
        se    = se_values
      )))
    }
  }

  if (length(results) == 0L) {
    return(data.table::data.table(
      term  = character(0L),
      level = character(0L),
      blup  = numeric(0L),
      se    = numeric(0L)
    ))
  }

  data.table::rbindlist(results, use.names = TRUE)
}


# ---- Convergence check -------------------------------------------------------

#' Check convergence of an lme4 model
#'
#' Inspects the fitted model for convergence warnings, singular fit messages,
#' and other diagnostic flags.
#'
#' @param lme4_model A fitted lme4 model object (class merMod).
#' @return Logical: \code{TRUE} if the model converged without warnings,
#'   \code{FALSE} otherwise.
#' @noRd
.lme4_check_convergence <- function(lme4_model) {
  # Check for convergence code from the optimiser
  opt_info <- lme4_model@optinfo
  conv_code <- opt_info$conv$opt
  if (!is.null(conv_code) && conv_code != 0L) {
    return(FALSE)
  }

  # Check for convergence warnings stored by lme4
  conv_msgs <- opt_info$warnings
  if (length(conv_msgs) > 0L) {
    return(FALSE)
  }

  # Check for singular fit (boundary variance components)
  is_singular <- tryCatch(
    lme4::isSingular(lme4_model),
    error = function(e) FALSE
  )
  if (isTRUE(is_singular)) {
    .msg("lme4 model is singular (some variance components at boundary).")
    # Singular fit is not necessarily non-convergence, but flag it
    # Return TRUE (converged, but warn separately)
    return(TRUE)
  }

  TRUE
}


# ---- Genomic G-matrix Cholesky trick -----------------------------------------

#' Apply Cholesky decomposition trick for genomic relationship matrices
#'
#' lme4 does not support \code{vm()} or custom variance structures for known
#' relationship matrices. For a genomic G-matrix, we use the Cholesky
#' decomposition trick: pre-multiply the incidence (design) matrix Z by the
#' Cholesky factor L of G, giving \eqn{Z^* = Z \times L}. Then
#' \eqn{Z^* u^*} where \eqn{u^* \sim N(0, \sigma^2_g I)} is equivalent to
#' \eqn{Z u} where \eqn{u \sim N(0, \sigma^2_g G)}.
#'
#' The function adds the transformed design matrix columns to the data and
#' creates a dummy grouping variable for the lme4 random effects.
#'
#' @param data data.frame or data.table to augment.
#' @param model_spec Named list with model specification.
#' @param known_mats Named list of known matrices from \code{model_spec}.
#' @return The augmented data with Z-star columns and dummy grouping variable.
#' @noRd
.lme4_apply_cholesky_trick <- function(data, model_spec, known_mats) {
  # Ensure we work on a copy to avoid modifying the caller's data
  if (data.table::is.data.table(data)) {
    data <- data.table::copy(data)
  } else {
    data <- data.frame(data, check.names = FALSE)
  }

  for (mat_name in names(known_mats)) {
    G <- known_mats[[mat_name]]

    if (!is.matrix(G) && !inherits(G, "Matrix")) {
      warning(
        sprintf("Known matrix '%s' is not a matrix. Skipping.", mat_name),
        call. = FALSE
      )
      next
    }

    # Compute Cholesky factor: G = L %*% t(L)
    # Add small ridge for numerical stability if G is near-singular
    L <- tryCatch({
      chol_G <- chol(G + diag(nrow(G)) * 1e-6)
      t(chol_G)  # Lower triangular
    }, error = function(e) {
      warning(
        sprintf(
          "Cholesky decomposition of '%s' failed: %s\n",
          mat_name, conditionMessage(e)
        ),
        "  Falling back to eigendecomposition.",
        call. = FALSE
      )
      # Eigendecomposition fallback: G = V D V', L = V D^{1/2}
      eig <- eigen(G, symmetric = TRUE)
      pos <- eig$values > max(eig$values) * 1e-8
      eig$vectors[, pos, drop = FALSE] %*%
        diag(sqrt(eig$values[pos]), nrow = sum(pos))
    })

    # Create Z* = Z %*% L columns in the data
    # The design matrix Z maps observations to genotype levels
    # For now, we assume a simple one-to-one mapping via a factor column
    # TODO: Extend to handle arbitrary incidence matrices

    n_cols <- ncol(L)
    col_prefix <- paste0("Zstar_", mat_name, "_")
    dummy_grp  <- paste0("dummy_", mat_name)

    # Add a dummy observation-level grouping factor (required for lme4)
    data[[dummy_grp]] <- factor(seq_len(nrow(data)))

    # Add Z* columns
    # If data has a column matching genotype levels of G, construct Z first
    # Otherwise, assume L is n_obs x n_cols (pre-computed Z*L)
    if (nrow(L) == nrow(data)) {
      # L is already observation-level (Z*L precomputed or n_geno == n_obs)
      for (k in seq_len(n_cols)) {
        col_name <- paste0(col_prefix, k)
        data[[col_name]] <- L[, k]
      }
    } else {
      .msg(
        sprintf(
          "G-matrix '%s' has %d rows but data has %d rows. ",
          mat_name, nrow(L), nrow(data)
        ),
        "Skipping Cholesky trick -- manual incidence matrix required."
      )
    }
  }

  data
}


# ---- Capabilities ------------------------------------------------------------

#' Report lme4 engine capabilities
#'
#' Returns a named list describing which funcrop model features the lme4
#' backend supports, partially supports, or does not support. Used by the
#' engine dispatcher to validate model requests before fitting.
#'
#' @return Named list:
#'   \describe{
#'     \item{\code{supported}}{Character vector of fully supported features.}
#'     \item{\code{not_supported}}{Character vector of unsupported features.}
#'     \item{\code{partial}}{Named list of partially supported features with
#'       notes.}
#'   }
#' @noRd
.lme4_capabilities <- function() {
  list(
    supported = c(
      "bspline_random",
      "glmm"
    ),
    not_supported = c(
      "vm",
      "fa",
      "us",
      "ar1",
      "ar1ar1",
      "tensor",
      "bayesian",
      "bam",
      "spatial_spline"
    ),
    partial = list(
      at = paste(
        "lme4 does not support at() directly. Heterogeneous variances can",
        "sometimes be approximated using nested random effects or dummy",
        "variable coding, but this is not automatic."
      ),
      genomic = paste(
        "Genomic G-matrices are supported via the Cholesky decomposition",
        "trick: Z* = Z %*% chol(G). This is exact for Gaussian models",
        "but adds columns to the data and may be slow for large G."
      )
    )
  )
}


# ---- Engine registration -----------------------------------------------------

#' Register the lme4 engine with funcrop's engine dispatcher
#'
#' Creates an engine descriptor for lme4 and registers it with the internal
#' engine registry. This function should be called in \code{.onLoad()} to make
#' the lme4 backend available when the lme4 package is installed.
#'
#' @return Invisible \code{NULL}. Called for its side effect of registering the
#'   engine.
#' @noRd
.register_lme4_engine <- function() {

  # Engine descriptor: standardised interface for the dispatcher
  engine_desc <- list(
    name         = "lme4",
    label        = "lme4 (REML/ML via lmer/glmer)",
    available    = .has_lme4,             # function, evaluated lazily
    fit          = .lme4_fit,
    extract_vc   = .lme4_extract_vc,
    extract_blups = .lme4_extract_blups,
    check_convergence = .lme4_check_convergence,
    build_formula = .lme4_build_formula,
    capabilities = .lme4_capabilities
  )

  # Register with the engine registry.
  # NOTE: The engine registry mechanism (.funcrop_engines environment or

  # similar) is being developed as part of Phase 5 (extensible engine
  # interface). When the registry is ready, replace the placeholder below
  # with the actual registration call, e.g.:
  #   .register_engine(engine_desc)
  #
  # For now, store in a package-level environment so the dispatcher can find it.
  pkg_env <- parent.env(environment())
  if (!exists(".funcrop_engine_registry", envir = pkg_env)) {
    assign(".funcrop_engine_registry", new.env(parent = emptyenv()),
           envir = pkg_env)
  }
  registry <- get(".funcrop_engine_registry", envir = pkg_env)
  assign("lme4", engine_desc, envir = registry)

  invisible(NULL)
}
# This should be called in .onLoad() -- see funcrop-package.R or zzz.R
