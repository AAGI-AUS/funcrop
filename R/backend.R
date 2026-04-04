# funcrop: Engine dispatcher -- dual-backend architecture
#
# Provides the user-facing API for selecting and querying estimation engines,
# plus the internal dispatch mechanism that routes model fitting to the
# appropriate backend (ASReml-R v4.2 or bayesreml).
#
# Engine auto-detection order: ASReml first (faster REML), bayesreml fallback.

# ---- Exported functions ------------------------------------------------------

#' List Available Estimation Engines
#'
#' Returns a character vector of estimation engines currently installed and
#' available for use with funcrop. At least one engine must be installed for
#' model fitting.
#'
#' @return Character vector of available engine names. Possible values are
#'   `"asreml"` (REML via ASReml-R v4.2+) and `"bayesreml"` (Bayesian MCMC
#'   via greta). Returns `character(0)` if no engines are installed.
#'
#' @details
#' funcrop supports two estimation backends:
#' \describe{
#'   \item{`asreml`}{ASReml-R v4.2+ -- restricted maximum likelihood (REML)
#'     estimation. Commercial licence required. Faster for large models.}
#'   \item{`bayesreml`}{bayesreml -- Bayesian MCMC estimation via greta/TensorFlow.
#'     Open-source. Provides full posterior distributions and credible intervals.}
#' }
#'
#' Both packages are listed in Suggests; neither is required at install time.
#' However, at least one must be installed before fitting any model.
#'
#' @examples
#' funcrop_engines()
#'
#' @export
funcrop_engines <- function() {
  engines <- character(0L)
  if (.has_asreml()) {
    engines <- c(engines, "asreml")
  }
  if (.has_bayesreml()) {
    engines <- c(engines, "bayesreml")
  }
  engines
}

#' Get or Set the Default Estimation Engine
#'
#' Query or set the default engine used by funcrop fitting functions. When
#' called without arguments, returns the current default. When called with an
#' engine name, sets it as the session default.
#'
#' @param engine Character string specifying the engine to use as default:
#'   `"asreml"` or `"bayesreml"`. If `NULL` (the default), returns the
#'   current default engine without changing it.
#'
#' @return Character string -- the current (or newly set) default engine name.
#'   Returned invisibly when setting.
#'
#' @details
#' The default engine is stored in `options("funcrop.engine")`. If the option
#' is not set, auto-detection is used: ASReml-R if installed, otherwise
#' bayesreml. If neither is installed, an error is raised.
#'
#' Setting an engine that is not installed raises an error immediately, so
#' downstream fitting functions can rely on the default being valid.
#'
#' @examples
#' # Query current default
#' \dontrun{
#' funcrop_default_engine()
#'
#' # Set bayesreml as default
#' funcrop_default_engine("bayesreml")
#' }
#'
#' @export
funcrop_default_engine <- function(engine = NULL) {
  if (is.null(engine)) {
    # --- Query mode ---
    current <- getOption("funcrop.engine", default = NULL)
    if (!is.null(current)) {
      # Validate it is still installed (user may have removed the package)
      if (!.engine_available(current)) {
        warning(
          sprintf(
            "Previously set engine '%s' is no longer available. Auto-detecting.",
            current
          ),
          call. = FALSE
        )
        options(funcrop.engine = NULL)
        current <- NULL
      }
    }
    if (is.null(current)) {
      current <- .auto_detect_engine()
    }
    return(current)
  }


  # --- Set mode ---
  .validate_engine(engine)
  options(funcrop.engine = engine)
  .msg(sprintf("funcrop default engine set to '%s'.", engine))
  invisible(engine)
}


# ---- Internal dispatch functions ---------------------------------------------

#' Resolve the engine to use for a given fitting call
#'
#' Handles the "auto" sentinel by checking the user option and falling back to
#' auto-detection. Validates the resolved engine is actually installed.
#'
#' @param engine Character string: `"auto"`, `"asreml"`, or `"bayesreml"`.
#' @return A single character string: `"asreml"` or `"bayesreml"`.
#' @noRd
.resolve_engine <- function(engine = c("auto", "asreml", "bayesreml")) {
  engine <- match.arg(engine)

  if (engine == "auto") {
    # Check user option first
    opt <- getOption("funcrop.engine", default = NULL)
    if (!is.null(opt) && .engine_available(opt)) {
      return(opt)
    }
    # Auto-detect
    return(.auto_detect_engine())
  }

  # Explicit engine requested -- validate it is installed
  .validate_engine(engine)
  engine
}


#' Central dispatch function for model fitting
#'
#' Routes the model fitting request to the appropriate backend based on the
#' resolved engine. Both backends receive a standardised `model_spec` list and
#' return a standardised raw result.
#'
#' @param engine Character string: `"asreml"` or `"bayesreml"` (already
#'   resolved, not `"auto"`).
#' @param model_spec Named list describing the model. Expected elements:
#'   \describe{
#'     \item{`fixed`}{Formula for fixed effects (e.g., `yield ~ lin(time)`).}
#'     \item{`random`}{Formula for random effects (e.g., `~ variety + variety:spline(time)`).}
#'     \item{`rcov`}{Residual covariance formula (e.g., `~ units`).}
#'     \item{`known_matrices`}{Named list of relationship/design matrices to be
#'       passed via `vm()` or equivalent (e.g., genomic relationship matrix G,
#'       B-spline design matrix Z_spline).}
#'     \item{`start_values`}{Optional named vector of starting values for
#'       variance parameters (ASReml only).}
#'     \item{`maxiter`}{Maximum iterations for REML (default 50).}
#'     \item{`workspace`}{ASReml workspace in bytes (default 128e6).}
#'   }
#' @param data A data.frame or data.table containing the response and predictor
#'   variables.
#' @param ... Additional arguments passed to the backend-specific fitting
#'   function (e.g., `mcmc_control` for bayesreml).
#'
#' @return A named list (raw result from the backend) containing at minimum:
#'   \describe{
#'     \item{`model`}{The fitted model object (asreml or bayesreml class).}
#'     \item{`converged`}{Logical: did the optimisation converge?}
#'     \item{`engine`}{Character string identifying which engine was used.}
#'   }
#'   Additional elements are engine-specific (e.g., `log_lik` for ASReml,
#'   `draws` and `summary` for bayesreml).
#'
#' @noRd
.dispatch_fit <- function(engine, model_spec, data, ...) {
  # Defensive checks
  stopifnot(
    is.character(engine), length(engine) == 1L,
    engine %in% c("asreml", "bayesreml")
  )
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

  .msg(sprintf("Fitting model using '%s' engine...", engine))

  raw_result <- switch(
    engine,
    asreml = .asreml_fit(model_spec = model_spec, data = data, ...),
    bayesreml = .bayesreml_fit(model_spec = model_spec, data = data, ...),
    stop(sprintf("Unknown engine '%s'.", engine), call. = FALSE)
  )

  # Tag the result with the engine used

raw_result[["engine"]] <- engine
  raw_result
}


# ---- Auto-detection helper ---------------------------------------------------

#' Auto-detect the best available engine
#'
#' Preference order: ASReml (faster REML) > bayesreml (Bayesian, open-source).
#' Raises an informative error if no engine is available.
#'
#' @return Character string: `"asreml"` or `"bayesreml"`.
#' @noRd
.auto_detect_engine <- function() {
  if (.has_asreml()) {
    return("asreml")
  }
  if (.has_bayesreml()) {
    return("bayesreml")
  }
  stop(
    "No estimation engine available.\n",
    "Install at least one of:\n",
    "  - asreml (>= 4.2.0): https://vsni.co.uk/software/asreml-r\n",
    "  - bayesreml (>= 0.1.0): install.packages('bayesreml')\n",
    "See ?funcrop_engines for details.",
    call. = FALSE
  )
}


#' Check if a specific engine is available
#'
#' @param engine Character string: `"asreml"` or `"bayesreml"`.
#' @return Logical scalar.
#' @noRd
.engine_available <- function(engine) {
  switch(
    engine,
    asreml    = .has_asreml(),
    bayesreml = .has_bayesreml(),
    FALSE
  )
}
