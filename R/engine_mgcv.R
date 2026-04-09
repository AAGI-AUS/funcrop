# funcrop: mgcv backend (gam / bam / gamm)
#
# All functions in this file are internal (not exported). They implement the
# interface between funcrop's model specification and mgcv's generalised
# additive model framework for REML estimation.
#
# mgcv provides three fitting routes:
#   - gam()  : standard GAM with REML/ML (default)
#   - bam()  : memory-efficient GAM for large datasets (>50k rows)
#   - gamm() : GAM with nlme correlation structures (e.g., AR1 residuals)
#
# IMPORTANT: Every call to mgcv:: or nlme:: functions is guarded by
# requireNamespace() so the package loads cleanly without mgcv installed.

# ---- Main fitting function ---------------------------------------------------

#' Fit a model via mgcv
#'
#' Translates funcrop's internal model specification into mgcv syntax and
#' dispatches to `mgcv::gam()`, `mgcv::bam()`, or `mgcv::gamm()` depending
#' on data size and model features requested (e.g., AR1 residuals).
#'
#' @param model_spec Named list with model specification. Expected elements:
#'   \describe{
#'     \item{`fixed`}{Formula for fixed effects (e.g., `yield ~ lin(time)`).}
#'     \item{`random`}{Formula for random effects (e.g.,
#'       `~ variety_f:Zrange_1 + variety_f:Zrange_2 + block`).}
#'     \item{`rcov`}{Optional residual covariance specification. If it
#'       contains AR1, gamm() with nlme correlation is used.}
#'     \item{`known_matrices`}{Named list of relationship/design matrices.}
#'   }
#' @param data data.frame or data.table of observations.
#' @param control Named list of additional control parameters:
#'   \describe{
#'     \item{`method`}{Smoothing parameter estimation method (default
#'       `"REML"`).}
#'     \item{`sub_engine`}{Force sub-engine: `"gam"`, `"bam"`, or `"gamm"`.
#'       If NULL (default), auto-selected.}
#'     \item{`bam_threshold`}{Row count threshold for switching to bam()
#'       (default 50000).}
#'     \item{`bam_discrete`}{Use discrete = TRUE in bam() (default TRUE).}
#'     \item{`ar1_form`}{One-sided formula for AR1 correlation grouping
#'       (e.g., `~ row | col`). Only used with gamm().}
#'     \item{`ar1_value`}{Starting value for AR1 parameter (default 0.3).}
#'   }
#'
#' @return Named list:
#'   \describe{
#'     \item{`model`}{The fitted mgcv model object (gam, bam, or gamm).}
#'     \item{`converged`}{Logical: did the optimisation converge?}
#'     \item{`log_lik`}{Log-likelihood at convergence.}
#'     \item{`n_iter`}{Number of iterations (where available).}
#'   }
#' @noRd
.mgcv_fit <- function(model_spec, data, control = list()) {
  if (!.has_mgcv()) {
    stop(
      "mgcv is not installed. Cannot use the 'mgcv' engine.\n",
      "Install via: install.packages('mgcv')",
      call. = FALSE
    )
  }

  # --- Validate inputs --------------------------------------------------------
  if (!is.list(model_spec)) {
    stop("`model_spec` must be a list.", call. = FALSE)
  }
  if (!"fixed" %in% names(model_spec)) {
    stop("`model_spec` must contain at least a `fixed` element.", call. = FALSE)
  }
  if (!is.data.frame(data) && !data.table::is.data.table(data)) {
    stop("`data` must be a data.frame or data.table.", call. = FALSE)
  }

  # Coerce data.table to data.frame for mgcv compatibility (mgcv sometimes

  # has issues with data.table class)
  if (data.table::is.data.table(data)) {
    data <- as.data.frame(data)
  }

  # --- Build the mgcv formula -------------------------------------------------
  mgcv_formula <- .mgcv_build_formula(model_spec, data)

  # --- Determine sub-engine ---------------------------------------------------
  method        <- control[["method"]]        %||% "REML"
  sub_engine    <- control[["sub_engine"]]     %||% NULL
  bam_threshold <- control[["bam_threshold"]] %||% 50000L
  bam_discrete  <- control[["bam_discrete"]]  %||% TRUE

  # Detect AR1 request from rcov specification
  has_ar1 <- .mgcv_detect_ar1(model_spec)

  if (is.null(sub_engine)) {
    if (has_ar1) {
      sub_engine <- "gamm"
    } else if (nrow(data) > bam_threshold) {
      sub_engine <- "bam"
    } else {
      sub_engine <- "gam"
    }
  }

  # Validate sub_engine choice
  sub_engine <- match.arg(sub_engine, c("gam", "bam", "gamm"))

  if (sub_engine == "gamm") {
    if (!requireNamespace("nlme", quietly = TRUE)) {
      stop(
        "nlme is required for gamm() fitting (AR1 residual correlation).\n",
        "Install via: install.packages('nlme')",
        call. = FALSE
      )
    }
  }

  .msg(sprintf("Fitting model using mgcv::%s() with method = '%s'...",
               sub_engine, method))

  # --- Fit the model ----------------------------------------------------------
  fit <- tryCatch(
    switch(sub_engine,
      gam = {
        mgcv::gam(
          formula  = mgcv_formula,
          data     = data,
          method   = method
        )
      },
      bam = {
        mgcv::bam(
          formula  = mgcv_formula,
          data     = data,
          method   = method,
          discrete = bam_discrete
        )
      },
      gamm = {
        # Build correlation structure for AR1
        cor_struct <- NULL
        if (has_ar1) {
          ar1_form  <- control[["ar1_form"]]  %||% NULL
          ar1_value <- control[["ar1_value"]] %||% 0.3

          if (is.null(ar1_form)) {
            # Attempt to parse from rcov specification
            ar1_form <- .mgcv_parse_ar1_form(model_spec[["rcov"]])
          }
          if (!is.null(ar1_form)) {
            cor_struct <- nlme::corAR1(value = ar1_value, form = ar1_form)
          }
        }

        mgcv::gamm(
          formula     = mgcv_formula,
          data        = data,
          method      = method,
          correlation = cor_struct
        )
      },
      stop(sprintf("Unknown mgcv sub-engine '%s'.", sub_engine), call. = FALSE)
    ),
    error = function(e) {
      stop(
        sprintf("mgcv::%s() fitting failed with error:\n  %s",
                sub_engine, conditionMessage(e)),
        call. = FALSE
      )
    }
  )

  # --- Extract convergence and summary info -----------------------------------
  converged <- .mgcv_check_convergence(fit)

  # Log-likelihood extraction
  log_lik <- tryCatch({
    if (sub_engine == "gamm") {
      stats::logLik(fit$lme)
    } else {
      stats::logLik(fit)
    }
  }, error = function(e) NA_real_)
  log_lik <- as.numeric(log_lik)


  # Iteration count
  n_iter <- tryCatch({
    if (sub_engine == "gamm") {
      # gamm uses lme underneath; extract iteration count
      fit$lme$numIter
    } else {
      # gam/bam: outer iteration count from mgcv
      if (!is.null(fit$outer.info)) {
        length(fit$outer.info$grad)
      } else {
        fit$iter
      }
    }
  }, error = function(e) NA_integer_)

  list(
    model      = fit,
    converged  = isTRUE(converged),
    log_lik    = log_lik,
    n_iter     = n_iter
  )
}


# ---- Formula construction ----------------------------------------------------

#' Translate funcrop model_spec to mgcv formula
#'
#' Parses the fixed and random formulas from funcrop's model specification
#' and constructs a single mgcv-compatible formula. Random effects are
#' translated into smooth terms:
#' \itemize{
#'   \item `variety_f:Zrange_k` becomes `s(variety_f, by = Zrange_k, bs = "re")`
#'   \item `block` becomes `s(block, bs = "re")`
#'   \item Tensor product terms become `te()` or `ti()` terms
#' }
#'
#' @param model_spec Named list with `fixed`, `random`, and optionally
#'   `known_matrices` elements.
#' @param data data.frame used to identify column types and validate terms.
#'
#' @return A formula suitable for `mgcv::gam()` / `mgcv::bam()` /
#'   `mgcv::gamm()`.
#' @noRd
.mgcv_build_formula <- function(model_spec, data) {
  # --- Fixed effects (LHS ~ RHS) ---
  fixed <- model_spec[["fixed"]]
  if (is.character(fixed)) {
    fixed <- stats::as.formula(fixed)
  }

  # Extract response and fixed predictors
  fixed_char <- deparse(fixed, width.cutoff = 500L)
  # Split into LHS and RHS
  parts <- strsplit(fixed_char, "~", fixed = TRUE)[[1L]]
  if (length(parts) != 2L) {
    stop("Fixed formula must be two-sided (response ~ predictors).", call. = FALSE)
  }
  response_str <- trimws(parts[1L])
  fixed_rhs    <- trimws(parts[2L])

  # --- Random effects -> mgcv smooth terms ---
  random <- model_spec[["random"]]
  smooth_terms <- character(0L)

  if (!is.null(random)) {
    if (is.character(random)) {
      random <- stats::as.formula(random)
    }
    # Parse the one-sided formula into individual terms
    random_char <- deparse(random, width.cutoff = 500L)
    random_char <- sub("^\\s*~\\s*", "", random_char)
    raw_terms   <- trimws(strsplit(random_char, "\\+")[[1L]])
    raw_terms   <- raw_terms[nzchar(raw_terms)]

    for (tm in raw_terms) {
      smooth_terms <- c(smooth_terms, .mgcv_translate_term(tm, data))
    }
  }

  # --- Combine into a single formula ---
  rhs_combined <- fixed_rhs
  if (length(smooth_terms) > 0L) {
    rhs_combined <- paste(
      c(rhs_combined, smooth_terms),
      collapse = " + "
    )
  }

  formula_str <- paste(response_str, "~", rhs_combined)
  stats::as.formula(formula_str, env = environment(fixed))
}


#' Translate a single random effect term to mgcv syntax
#'
#' Handles the following patterns:
#' \itemize{
#'   \item `variety_f:Zrange_k` -> `s(variety_f, by = Zrange_k, bs = "re")`
#'   \item `block` (factor) -> `s(block, bs = "re")`
#'   \item `te(x, z)` or `ti(x, z)` -> pass through (already mgcv syntax)
#'   \item `str(...)` or `vm(...)` -> warning + skip (ASReml-specific)
#' }
#'
#' @param term Character string: a single term from the random formula.
#' @param data data.frame for column type checking.
#' @return Character string: the mgcv-compatible smooth term.
#' @noRd
.mgcv_translate_term <- function(term, data) {
  term <- trimws(term)

  # --- Pass-through: already mgcv smooth syntax ---
  if (grepl("^s\\(|^te\\(|^ti\\(|^t2\\(", term)) {
    return(term)
  }

  # --- ASReml-specific terms: warn and skip ---
  if (grepl("^str\\(|^vm\\(|^us\\(|^fa\\(|^corh?\\(", term)) {
    warning(
      sprintf(
        "mgcv backend does not support ASReml-specific term '%s'. Skipping.",
        term
      ),
      call. = FALSE
    )
    return(NULL)
  }

  # --- Interaction term: factor:continuous -> s(factor, by = continuous, bs = "re") ---
  if (grepl(":", term, fixed = TRUE)) {
    parts <- trimws(strsplit(term, ":", fixed = TRUE)[[1L]])

    if (length(parts) == 2L) {
      # Determine which part is the factor and which is the continuous covariate.
      # Convention: Zrange_k / Zspline_k / Bsp_k are continuous (numeric columns
      # in data or known_matrices); variety_f / block / env are factors.
      factor_part <- NULL
      cont_part   <- NULL

      for (p in parts) {
        if (p %in% names(data)) {
          if (is.factor(data[[p]]) || is.character(data[[p]])) {
            factor_part <- p
          } else if (is.numeric(data[[p]])) {
            cont_part <- p
          }
        } else {
          # Column not in data — use naming heuristics
          if (grepl("^Z(range|spline)_|^Bsp_|^bs\\(", p, ignore.case = TRUE)) {
            cont_part <- p
          } else {
            # Default: assume it is a factor grouping variable
            factor_part <- p
          }
        }
      }

      if (!is.null(factor_part) && !is.null(cont_part)) {
        return(sprintf('s(%s, by = %s, bs = "re")', factor_part, cont_part))
      }

      # Fallback for two-factor interaction: s(f1, f2, bs = "re")
      if (is.null(cont_part)) {
        return(sprintf('s(%s, %s, bs = "re")', parts[1L], parts[2L]))
      }
    }

    # Three-way or higher interaction: attempt generic translation
    warning(
      sprintf(
        "mgcv translation for multi-way interaction '%s' is approximate. ",
        term
      ),
      call. = FALSE
    )
    all_parts <- paste(parts, collapse = ", ")
    return(sprintf('s(%s, bs = "re")', all_parts))
  }

  # --- Simple term (single variable) -> s(term, bs = "re") for factors ---
  # For a single factor (e.g., block, variety_f), use random effect smooth.
  if (term %in% names(data)) {
    if (is.factor(data[[term]]) || is.character(data[[term]])) {
      return(sprintf('s(%s, bs = "re")', term))
    }
    # Numeric: could be a random slope — use bs = "re" as well
    return(sprintf('s(%s, bs = "re")', term))
  }

  # Term not found in data — attempt as-is (may be constructed downstream)
  sprintf('s(%s, bs = "re")', term)
}


# ---- Variance component extraction -------------------------------------------

#' Extract variance components from an mgcv model
#'
#' For gam/bam: extracts variance components from `mgcv::gam.vcomp()` or
#' from the smoothing parameter / scale parameter.
#' For gamm: extracts from the nlme `lme` component via `nlme::VarCorr()`.
#'
#' @param mgcv_model A fitted mgcv model object (gam, bam, or gamm list).
#' @return A data.table with columns: `component`, `estimate`, `se`,
#'   `z_ratio`, `bound`.
#' @noRd
.mgcv_extract_vc <- function(mgcv_model) {
  if (!.has_mgcv()) {
    stop("mgcv is required for variance component extraction.", call. = FALSE)
  }

  # --- gamm returns a list with $gam and $lme ---
  if (is.list(mgcv_model) && !is.null(mgcv_model[["lme"]])) {
    return(.mgcv_extract_vc_gamm(mgcv_model))
  }

  # --- gam / bam ---
  # Try gam.vcomp() first (available in mgcv >= 1.8)
  vc <- tryCatch({
    vc_raw <- mgcv::gam.vcomp(mgcv_model, rescale = FALSE)
    # gam.vcomp returns a matrix with rows = components
    if (is.matrix(vc_raw)) {
      data.table::data.table(
        component = rownames(vc_raw),
        estimate  = vc_raw[, 1L],
        se        = if (ncol(vc_raw) >= 2L) vc_raw[, 2L] else NA_real_,
        z_ratio   = NA_real_,
        bound     = ""
      )
    } else if (is.numeric(vc_raw)) {
      # Named vector (older mgcv versions)
      data.table::data.table(
        component = names(vc_raw),
        estimate  = unname(vc_raw),
        se        = NA_real_,
        z_ratio   = NA_real_,
        bound     = ""
      )
    } else {
      NULL
    }
  }, error = function(e) NULL)

  if (!is.null(vc)) return(vc)

  # Fallback: extract from smoothing parameters and scale
  sp_vals <- mgcv_model$sp
  sig2    <- mgcv_model$sig2

  components <- character(0L)
  estimates  <- numeric(0L)

  if (length(sp_vals) > 0L) {
    # Smoothing parameters are inverse variance components: sp = 1/sigma^2
    # Variance component = scale / sp (for penalised terms)
    vc_from_sp <- sig2 / sp_vals
    components <- c(components, names(sp_vals))
    estimates  <- c(estimates, vc_from_sp)
  }

  # Add residual variance
  components <- c(components, "residual")
  estimates  <- c(estimates, sig2)

  data.table::data.table(
    component = components,
    estimate  = estimates,
    se        = NA_real_,
    z_ratio   = NA_real_,
    bound     = ""
  )
}


#' Extract variance components from a gamm model (nlme lme component)
#'
#' @param gamm_model A gamm result list with `$lme` and `$gam` components.
#' @return data.table with columns: component, estimate, se, z_ratio, bound.
#' @noRd
.mgcv_extract_vc_gamm <- function(gamm_model) {
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("nlme is required for gamm variance component extraction.",
         call. = FALSE)
  }

  lme_model <- gamm_model[["lme"]]

  # VarCorr returns a character matrix with variance/SD info
  vc_raw <- tryCatch(
    nlme::VarCorr(lme_model),
    error = function(e) {
      warning("Failed to extract VarCorr from gamm lme component: ",
              conditionMessage(e), call. = FALSE)
      NULL
    }
  )

  if (is.null(vc_raw)) {
    return(data.table::data.table(
      component = character(0L),
      estimate  = numeric(0L),
      se        = numeric(0L),
      z_ratio   = numeric(0L),
      bound     = character(0L)
    ))
  }

  # VarCorr is a character matrix; parse Variance column
  comp_names <- rownames(vc_raw)
  var_col    <- suppressWarnings(as.numeric(vc_raw[, "Variance"]))

  # Remove NA rows (headers / labels)
  valid <- !is.na(var_col)
  comp_names <- comp_names[valid]
  var_vals   <- var_col[valid]

  data.table::data.table(
    component = comp_names,
    estimate  = var_vals,
    se        = NA_real_,
    z_ratio   = NA_real_,
    bound     = ""
  )
}


# ---- BLUP extraction ---------------------------------------------------------

#' Extract random effects (BLUPs) from an mgcv model
#'
#' For gam/bam: uses `predict(model, type = "terms")` to extract smooth
#' term contributions, or `mgcv::predict.gam()` with `type = "iterms"`.
#' For gamm: uses `nlme::ranef()` from the lme component.
#'
#' @param mgcv_model A fitted mgcv model object (gam, bam, or gamm list).
#' @param terms Character vector of random terms to extract. These should
#'   match the smooth term labels in the model (e.g., `"s(block)"`,
#'   `"s(variety_f,by=Zrange_1)"`).
#'
#' @return A data.table with columns: `term`, `level`, `blup`, `se`.
#' @noRd
.mgcv_extract_blups <- function(mgcv_model, terms) {
  if (!.has_mgcv()) {
    stop("mgcv is required for BLUP extraction.", call. = FALSE)
  }

  # --- gamm: extract from lme component ---
  if (is.list(mgcv_model) && !is.null(mgcv_model[["lme"]])) {
    return(.mgcv_extract_blups_gamm(mgcv_model, terms))
  }

  # --- gam / bam: extract smooth contributions ---
  results <- vector("list", length(terms))

  # Get all smooth term labels from the model
  smooth_labels <- vapply(mgcv_model$smooth, function(s) s$label,
                          character(1L))

  for (i in seq_along(terms)) {
    tm <- terms[i]

    # Find matching smooth(s) — allow partial matching
    matched_idx <- which(
      grepl(tm, smooth_labels, fixed = TRUE) |
      smooth_labels == tm
    )

    if (length(matched_idx) == 0L) {
      # Try matching by component variable names
      matched_idx <- which(vapply(mgcv_model$smooth, function(s) {
        any(grepl(tm, s$term, fixed = TRUE))
      }, logical(1L)))
    }

    if (length(matched_idx) == 0L) next

    for (mi in matched_idx) {
      sm <- mgcv_model$smooth[[mi]]
      label <- sm$label

      # Extract coefficients for this smooth
      first_par <- sm$first.para
      last_par  <- sm$last.para
      coef_idx  <- first_par:last_par

      coefs <- stats::coef(mgcv_model)[coef_idx]
      coef_names <- names(coefs)
      if (is.null(coef_names)) {
        coef_names <- paste0(label, "_", seq_along(coefs))
      }

      # Standard errors from the Bayesian posterior covariance
      se_vals <- tryCatch({
        Vp <- mgcv_model$Vp  # Bayesian posterior covariance of coefficients
        if (!is.null(Vp)) {
          sqrt(diag(Vp)[coef_idx])
        } else {
          rep(NA_real_, length(coef_idx))
        }
      }, error = function(e) rep(NA_real_, length(coef_idx)))

      results[[length(results) + 1L]] <- data.table::data.table(
        term  = label,
        level = coef_names,
        blup  = unname(coefs),
        se    = se_vals
      )
    }
  }

  # Remove NULL entries
  results <- Filter(Negate(is.null), results)

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


#' Extract random effects from a gamm model (lme component)
#'
#' @param gamm_model A gamm result list.
#' @param terms Character vector of terms to extract.
#' @return data.table with columns: term, level, blup, se.
#' @noRd
.mgcv_extract_blups_gamm <- function(gamm_model, terms) {
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("nlme is required for gamm BLUP extraction.", call. = FALSE)
  }

  lme_model <- gamm_model[["lme"]]
  re <- tryCatch(
    nlme::ranef(lme_model),
    error = function(e) {
      warning("Failed to extract random effects from gamm: ",
              conditionMessage(e), call. = FALSE)
      NULL
    }
  )

  if (is.null(re)) {
    return(data.table::data.table(
      term  = character(0L),
      level = character(0L),
      blup  = numeric(0L),
      se    = numeric(0L)
    ))
  }

  results <- vector("list", length(terms))

  # ranef() returns a list (one element per grouping factor) or a data.frame
  if (is.data.frame(re)) {
    # Single grouping factor
    re <- list(re)
    names(re) <- "default"
  }

  for (i in seq_along(terms)) {
    tm <- terms[i]

    for (grp_name in names(re)) {
      if (!grepl(tm, grp_name, fixed = TRUE) && tm != grp_name) next

      grp_re <- re[[grp_name]]
      if (is.data.frame(grp_re)) {
        for (col in names(grp_re)) {
          results[[length(results) + 1L]] <- data.table::data.table(
            term  = paste0(grp_name, ":", col),
            level = rownames(grp_re),
            blup  = grp_re[[col]],
            se    = NA_real_
          )
        }
      }
    }
  }

  results <- Filter(Negate(is.null), results)

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


# ---- Convergence checking ----------------------------------------------------

#' Check convergence of an mgcv model
#'
#' For gam/bam: checks `model$converged` or `model$outer.info$conv`.
#' For gamm: checks convergence of the underlying lme component.
#'
#' @param mgcv_model A fitted mgcv model (gam, bam, or gamm list).
#' @return Logical scalar: TRUE if converged, FALSE otherwise.
#' @noRd
.mgcv_check_convergence <- function(mgcv_model) {
  # --- gamm: check lme convergence ---
  if (is.list(mgcv_model) && !is.null(mgcv_model[["lme"]])) {
    lme_model <- mgcv_model[["lme"]]
    # lme stores convergence info differently; check for error status
    converged <- tryCatch({
      # If lme converged, apVar is not a character error message
      apvar <- lme_model$apVar
      !is.character(apvar)
    }, error = function(e) FALSE)
    return(isTRUE(converged))
  }

  # --- gam / bam ---
  # Method 1: direct converged flag
  if (!is.null(mgcv_model$converged)) {
    return(isTRUE(mgcv_model$converged))
  }

  # Method 2: outer iteration convergence info
  if (!is.null(mgcv_model$outer.info)) {
    conv <- mgcv_model$outer.info$conv
    if (!is.null(conv)) {
      # conv is a character string; "full convergence" or similar
      if (is.character(conv)) {
        return(grepl("converge", conv, ignore.case = TRUE))
      }
      return(isTRUE(conv))
    }
  }

  # Method 3: check mgcv.conv (some versions)
  if (!is.null(mgcv_model$mgcv.conv)) {
    fully_converged <- mgcv_model$mgcv.conv$fully.converged
    if (!is.null(fully_converged)) {
      return(isTRUE(fully_converged))
    }
  }

  # If no convergence info found, assume converged if model exists
  # (mgcv throws errors on non-convergence for most failures)
  TRUE
}


# ---- AR1 detection helpers ---------------------------------------------------

#' Detect whether the model specification requests AR1 residual correlation
#'
#' Checks the `rcov` element of model_spec for AR1-related keywords.
#'
#' @param model_spec Named list with model specification.
#' @return Logical scalar.
#' @noRd
.mgcv_detect_ar1 <- function(model_spec) {
  rcov <- model_spec[["rcov"]]
  if (is.null(rcov)) return(FALSE)

  rcov_str <- if (is.character(rcov)) {
    rcov
  } else if (inherits(rcov, "formula")) {
    deparse(rcov, width.cutoff = 500L)
  } else {
    return(FALSE)
  }

  grepl("ar1|AR1|corAR1|ar\\.1", rcov_str, ignore.case = FALSE)
}


#' Parse AR1 correlation form from rcov specification
#'
#' Attempts to extract the grouping formula for nlme::corAR1() from
#' funcrop's rcov specification. Expected formats:
#' \itemize{
#'   \item `~ ar1(row):id(col)` -> form = `~ row | col`
#'   \item `~ ar1(time, col)` -> form = `~ time | col`
#' }
#'
#' @param rcov Formula or character string with rcov specification.
#' @return A one-sided formula for corAR1(form = ...) or NULL.
#' @noRd
.mgcv_parse_ar1_form <- function(rcov) {
  if (is.null(rcov)) return(NULL)

  rcov_str <- if (inherits(rcov, "formula")) {
    deparse(rcov, width.cutoff = 500L)
  } else if (is.character(rcov)) {
    rcov
  } else {
    return(NULL)
  }

  # Remove the leading "~"
  rcov_str <- sub("^\\s*~\\s*", "", rcov_str)

  # Pattern: ar1(var1):id(var2) -> ~ var1 | var2
  match <- regmatches(
    rcov_str,
    regexec("ar1\\(([^)]+)\\)\\s*:\\s*id\\(([^)]+)\\)", rcov_str,
            ignore.case = TRUE)
  )[[1L]]

  if (length(match) == 3L) {
    return(stats::as.formula(paste("~", match[2L], "|", match[3L])))
  }

  # Pattern: ar1(var1, var2) -> ~ var1 | var2
  match2 <- regmatches(
    rcov_str,
    regexec("ar1\\(\\s*([^,]+)\\s*,\\s*([^)]+)\\s*\\)", rcov_str,
            ignore.case = TRUE)
  )[[1L]]

  if (length(match2) == 3L) {
    return(stats::as.formula(paste("~", trimws(match2[2L]), "|",
                                   trimws(match2[3L]))))
  }

  # Could not parse — return NULL and let the user provide ar1_form in control
  warning(
    "Could not parse AR1 grouping structure from rcov specification.\n",
    "Provide `control = list(ar1_form = ~ time | group)` explicitly.",
    call. = FALSE
  )
  NULL
}


# ---- Capabilities ------------------------------------------------------------

#' Report mgcv engine capabilities
#'
#' Returns a named list describing which funcrop model features are supported,
#' partially supported, or unsupported by the mgcv backend. Used by the
#' engine dispatcher for capability-based routing.
#'
#' @return Named list with elements:
#'   \describe{
#'     \item{`engine`}{Character: `"mgcv"`.}
#'     \item{`supported`}{Character vector of fully supported features.}
#'     \item{`partial`}{Named character vector of partially supported features
#'       with notes.}
#'     \item{`unsupported`}{Character vector of unsupported features.}
#'   }
#' @noRd
.mgcv_capabilities <- function() {
  list(
    engine = "mgcv",
    supported = c(
      "bspline_random",
      "tensor",
      "spatial_spline",
      "glmm",
      "bam"
    ),
    partial = c(
      ar1     = "AR1 residuals via gamm() only (requires nlme)",
      genomic = "Not natively supported; user must pre-compute and supply columns"
    ),
    unsupported = c(
      "vm",
      "fa",
      "us",
      "bayesian"
    )
  )
}


# ---- Engine registration -----------------------------------------------------

#' Register the mgcv engine with funcrop's engine dispatcher
#'
#' Called during package initialisation (or on-demand) to register mgcv as
#' an available estimation engine. Populates the internal engine registry
#' with fitting, extraction, and capability functions.
#'
#' @return Invisible NULL. Side effect: registers mgcv in the engine registry.
#' @noRd
.register_mgcv_engine <- function() {
  # Only register if mgcv is available

  if (!.has_mgcv()) {
    return(invisible(NULL))
  }

  # Build registration entry
  engine_entry <- list(
    name         = "mgcv",
    fit          = .mgcv_fit,
    extract_vc   = .mgcv_extract_vc,
    extract_blups = .mgcv_extract_blups,
    check_convergence = .mgcv_check_convergence,
    build_formula = .mgcv_build_formula,
    capabilities  = .mgcv_capabilities
  )

  # Register in the package-level engine registry (stored in an environment)
  # The registry is expected to be a named list in .funcrop_env$engines
  if (exists(".funcrop_env", envir = asNamespace("funcrop"),
             inherits = FALSE)) {
    env <- get(".funcrop_env", envir = asNamespace("funcrop"))
    if (is.null(env[["engines"]])) {
      env[["engines"]] <- list()
    }
    env[["engines"]][["mgcv"]] <- engine_entry
  }

  .msg("mgcv engine registered.")
  invisible(NULL)
}
