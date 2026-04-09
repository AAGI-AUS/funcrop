# funcrop: Engine registration infrastructure (v0.2.0)
#
# Provides the registration mechanism for estimation engines. Each engine
# registers itself with a standard set of methods. This file does NOT
# provide dispatch or resolution — those remain in backend.R for backward
# compatibility. This file provides:
#   - .funcrop_engine_registry (environment for registered engines)
#   - .register_engine() (registration function)
#   - .known_engine_features (feature vocabulary)
#   - .engines_supporting() (capability query)

# ---- Engine registry (package-level environment) ----------------------------

# Created at source time. Each entry is a named list of engine methods.
.funcrop_engine_registry <- new.env(parent = emptyenv())

# Known features that engines can declare support for.
.known_engine_features <- c(
  "bspline_random", "vm", "at", "fa", "us",
  "ar1", "ar1ar1", "tensor", "genomic", "bayesian",
  "glmm", "bam", "spatial_spline"
)


# ---- Engine registration ----------------------------------------------------

#' Register an estimation engine
#'
#' Adds an engine to the funcrop registry with a standard set of methods.
#' Called internally by backend initialisers or by users supplying custom
#' engines.
#'
#' @param name Character string. Unique engine identifier.
#' @param fit Function: `function(model_spec, data, control, ...)`.
#' @param extract_vc Function: `function(model)` -> data.table of VCs.
#' @param extract_blups Function: `function(model, terms)` -> data.table.
#' @param extract_fitted Function: `function(model, newdata)` -> fitted.
#' @param predict_fn Function: `function(model, newdata)` -> predictions.
#' @param convergence Function: `function(model)` -> list(converged, ...).
#' @param loglik Function: `function(model)` -> numeric.
#' @param capabilities Function: `function()` -> list(name, supported).
#'
#' @return Invisible TRUE on success.
#' @noRd
.register_engine <- function(name, fit, extract_vc, extract_blups,
                             extract_fitted, predict_fn, convergence,
                             loglik, capabilities) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop("`name` must be a non-empty single character string.", call. = FALSE)
  }
  name <- tolower(name)

  required_methods <- list(
    fit            = fit,
    extract_vc     = extract_vc,
    extract_blups  = extract_blups,
    extract_fitted = extract_fitted,
    predict_fn     = predict_fn,
    convergence    = convergence,
    loglik         = loglik,
    capabilities   = capabilities
  )

  for (method_name in names(required_methods)) {
    if (!is.function(required_methods[[method_name]])) {
      stop(sprintf("Engine '%s': `%s` must be a function.", name, method_name),
           call. = FALSE)
    }
  }

  # Validate capabilities output structure
  caps <- capabilities()
  if (!is.list(caps) || !all(c("name", "supported") %in% names(caps))) {
    stop(sprintf(
      "Engine '%s': capabilities() must return list with 'name' and 'supported'.",
      name), call. = FALSE)
  }

  # Warn about unknown features (typo protection)
  unknown <- setdiff(caps[["supported"]], .known_engine_features)
  if (length(unknown) > 0L) {
    warning(sprintf("Engine '%s' declares unknown features: %s.",
                    name, paste(unknown, collapse = ", ")),
            call. = FALSE)
  }

  assign(name, required_methods, envir = .funcrop_engine_registry)
  invisible(TRUE)
}


# ---- Capability query -------------------------------------------------------

#' Find engines supporting a given set of features
#'
#' Returns names of all registered engines whose declared capabilities
#' include all requested features.
#'
#' @param features Character vector of feature strings.
#' @return Character vector of engine names (possibly empty).
#' @noRd
.engines_supporting <- function(features) {
  if (!is.character(features) || length(features) == 0L) {
    return(character(0L))
  }
  registered <- ls(.funcrop_engine_registry)
  matching <- character(0L)
  for (eng_name in registered) {
    eng <- get(eng_name, envir = .funcrop_engine_registry, inherits = FALSE)
    caps <- eng[["capabilities"]]()
    if (all(features %in% caps[["supported"]])) {
      matching <- c(matching, eng_name)
    }
  }
  matching
}
