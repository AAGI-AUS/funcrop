# funcrop: S3 class definitions
# Defines fda_data, fda_model, and fda_comparison classes.

# ==============================================================================
# fda_data -- functional observation data container
# ==============================================================================

#' Create a functional data object
#'
#' Constructs a validated `fda_data` object wrapping functional observations
#' (e.g., repeated measurements of a secondary trait over time) into a
#' structured [data.table::data.table] with metadata. Designed for crop variety
#' trial data where secondary functional traits (NDVI, canopy temperature,
#' grain-fill rate) are related to a primary trait (e.g., yield).
#'
#' @param time Numeric vector of time points (or other continuous dimension
#'   such as thermal time, depth, wavelength).
#' @param value Numeric vector of measured values (the secondary / functional
#'   trait, e.g., NDVI).
#' @param id Character or factor identifying the observational unit (e.g.,
#'   plot, plant).
#' @param group Optional character or factor for grouping (e.g., variety,
#'   genotype). Default `NULL`.
#' @param spatial_row Optional integer vector of spatial row positions. Default
#'   `NULL`.
#' @param spatial_col Optional integer vector of spatial column positions.
#'   Default `NULL`.
#' @param primary_trait Optional numeric vector of primary trait values (e.g.,
#'   yield). Must contain exactly one value per unique `id`. Recycled
#'   internally to match the observation-level data. Default `NULL`.
#' @param primary_trait_name Optional character string naming the primary trait
#'   (e.g., `"yield_t_ha"`). Default `NULL`.
#' @param ... Additional named vectors of the same length as `time`, stored as
#'   extra columns in the underlying data.table.
#'
#' @return An object of class `fda_data`, inheriting from `data.table` and
#'   `data.frame`. Contains an attribute `"fda_meta"` with metadata (trait
#'   names, id column, etc.).
#'
#' @section Validation:
#' - `time`, `value`, and `id` must be non-NA and of equal length.
#' - `time` and `value` must be numeric and finite.
#' - If supplied, `primary_trait` must have exactly one value per unique `id`.
#' - Spatial columns must be integer-valued.
#'
#' @examples
#' # Simulate simple functional data
#' set.seed(42)
#' ids <- rep(paste0("plot_", 1:5), each = 10)
#' times <- rep(seq(100, 190, by = 10), times = 5)
#' vals <- sin(times / 50) + rnorm(50, sd = 0.1)
#' grps <- rep(c("var_A", "var_B", "var_A", "var_B", "var_A"), each = 10)
#' yields <- rep(c(3.2, 4.1, 3.8, 4.5, 3.0), each = 10)
#'
#' fd <- fda_data(
#'   time = times,
#'   value = vals,
#'   id = ids,
#'   group = grps,
#'   primary_trait = yields,
#'   primary_trait_name = "yield_t_ha"
#' )
#' print(fd)
#' summary(fd)
#'
#' @export
fda_data <- function(time, value, id, group = NULL, spatial_row = NULL,
                     spatial_col = NULL, primary_trait = NULL,
                     primary_trait_name = NULL, ...) {

  # --- Validate core vectors --------------------------------------------------
  n <- length(time)
  .validate_numeric(time, "time")
  .validate_numeric(value, "value")
  id <- .validate_factor(id, "id")

  if (length(value) != n) {
    stop(sprintf("`value` length (%d) must match `time` length (%d).",
                 length(value), n), call. = FALSE)
  }
  if (length(id) != n) {
    stop(sprintf("`id` length (%d) must match `time` length (%d).",
                 length(id), n), call. = FALSE)
  }

  # Build data.table
  dt <- data.table::data.table(
    id    = id,
    time  = time,
    value = value
  )

  # --- Optional: group --------------------------------------------------------
  if (!is.null(group)) {
    group <- .validate_factor(group, "group")
    if (length(group) != n) {
      stop(sprintf("`group` length (%d) must match `time` length (%d).",
                   length(group), n), call. = FALSE)
    }
    dt[, group := group]
  }

  # --- Optional: spatial coordinates ------------------------------------------
  if (!is.null(spatial_row)) {
    if (length(spatial_row) != n) {
      stop(sprintf("`spatial_row` length (%d) must match `time` length (%d).",
                   length(spatial_row), n), call. = FALSE)
    }
    .validate_numeric(spatial_row, "spatial_row")
    if (!all(spatial_row == as.integer(spatial_row))) {
      stop("`spatial_row` must contain integer values.", call. = FALSE)
    }
    dt[, spatial_row := as.integer(spatial_row)]
  }

  if (!is.null(spatial_col)) {
    if (length(spatial_col) != n) {
      stop(sprintf("`spatial_col` length (%d) must match `time` length (%d).",
                   length(spatial_col), n), call. = FALSE)
    }
    .validate_numeric(spatial_col, "spatial_col")
    if (!all(spatial_col == as.integer(spatial_col))) {
      stop("`spatial_col` must contain integer values.", call. = FALSE)
    }
    dt[, spatial_col := as.integer(spatial_col)]
  }

  # --- Optional: primary trait ------------------------------------------------
  if (!is.null(primary_trait)) {
    .validate_numeric(primary_trait, "primary_trait")
    uid <- unique(id)
    n_ids <- length(uid)

    if (length(primary_trait) == n) {
      # Check: one unique value per id
      pt_dt <- data.table::data.table(id = id, pt = primary_trait)
      n_unique_per_id <- pt_dt[, data.table::uniqueN(pt), by = id]
      bad <- n_unique_per_id[V1 > 1L]
      if (nrow(bad) > 0L) {
        stop(sprintf(
          "`primary_trait` must have a single value per `id`. Offending ids: %s.",
          paste(head(bad$id, 5L), collapse = ", ")), call. = FALSE)
      }
      dt[, primary_trait := primary_trait]
    } else if (length(primary_trait) == n_ids) {
      # Merge by id -- assume same order as unique(id)
      pt_map <- data.table::data.table(id = uid, primary_trait = primary_trait)
      dt <- merge(dt, pt_map, by = "id", sort = FALSE)
    } else {
      stop(sprintf(
        paste0("`primary_trait` length (%d) must equal either the number of ",
               "observations (%d) or the number of unique ids (%d)."),
        length(primary_trait), n, n_ids), call. = FALSE)
    }
  }

  if (!is.null(primary_trait_name)) {
    if (!is.character(primary_trait_name) || length(primary_trait_name) != 1L) {
      stop("`primary_trait_name` must be a single character string.",
           call. = FALSE)
    }
  }

  # --- Extra columns via ... --------------------------------------------------
  dots <- list(...)
  if (length(dots) > 0L) {
    for (nm in names(dots)) {
      if (is.null(nm) || nm == "") {
        stop("All arguments in `...` must be named.", call. = FALSE)
      }
      vec <- dots[[nm]]
      if (length(vec) != n) {
        stop(sprintf("Extra column `%s` length (%d) must match `time` (%d).",
                     nm, length(vec), n), call. = FALSE)
      }
      data.table::set(dt, j = nm, value = vec)
    }
  }

  # --- Set key for efficient operations ---------------------------------------
  data.table::setkeyv(dt, c("id", "time"))

  # --- Attach metadata --------------------------------------------------------
  meta <- list(
    primary_trait_name = primary_trait_name,
    has_primary_trait  = !is.null(primary_trait),
    has_group          = !is.null(group),
    has_spatial        = !is.null(spatial_row) && !is.null(spatial_col),
    n_ids              = data.table::uniqueN(dt$id),
    n_timepoints       = data.table::uniqueN(dt$time),
    time_range         = range(dt$time),
    value_range        = range(dt$value),
    extra_cols         = names(dots)
  )

  data.table::setattr(dt, "fda_meta", meta)
  data.table::setattr(dt, "class", c("fda_data", "data.table", "data.frame"))
  dt
}

# ---- fda_data methods --------------------------------------------------------

#' Print method for fda_data
#'
#' @param x An `fda_data` object.
#' @param ... Additional arguments (ignored).
#' @return Invisible `x`.
#' @export
print.fda_data <- function(x, ...) {
  meta <- attr(x, "fda_meta")
  cat("-- fda_data ---------------------------------\n")
  cat(sprintf("  Observations : %d rows, %d unique ids\n",
              nrow(x), meta$n_ids))
  cat(sprintf("  Time range   : [%.2f, %.2f] (%d unique points)\n",
              meta$time_range[1L], meta$time_range[2L], meta$n_timepoints))
  cat(sprintf("  Value range  : [%.4f, %.4f]\n",
              meta$value_range[1L], meta$value_range[2L]))
  if (meta$has_group) {
    grp_levels <- levels(x$group)
    n_grp <- length(grp_levels)
    shown <- paste(head(grp_levels, 5L), collapse = ", ")
    if (n_grp > 5L) shown <- paste0(shown, ", ...")
    cat(sprintf("  Groups       : %d (%s)\n", n_grp, shown))
  }
  if (meta$has_primary_trait) {
    pt_nm <- meta$primary_trait_name %||% "primary_trait"
    cat(sprintf("  Primary trait: %s\n", pt_nm))
  }
  if (meta$has_spatial) {
    cat("  Spatial      : row + col coordinates present\n")
  }
  cat("---------------------------------------------\n")
  # Print first few rows via data.table method
  NextMethod()
  invisible(x)
}

#' Summary method for fda_data
#'
#' @param object An `fda_data` object.
#' @param ... Additional arguments (ignored).
#' @return Invisible list of summary statistics.
#' @export
summary.fda_data <- function(object, ...) {
  meta <- attr(object, "fda_meta")
  # Work on a plain data.table to avoid [.fda_data dispatch issues
  dt <- as.data.table.fda_data(object)

  cat("== fda_data summary ==========================\n\n")

  cat(sprintf("Observations : %d\n", nrow(dt)))
  cat(sprintf("Unique ids   : %d\n", meta$n_ids))
  cat(sprintf("Time points  : %d unique in [%.2f, %.2f]\n",
              meta$n_timepoints, meta$time_range[1L], meta$time_range[2L]))

  # Per-id observation counts
  obs_per_id <- dt[, .N, by = id]$N
  cat(sprintf("Obs per id   : median = %d, range = [%d, %d]\n",
              as.integer(stats::median(obs_per_id)),
              min(obs_per_id), max(obs_per_id)))

  cat(sprintf("\nValue (secondary trait):\n"))
  cat(sprintf("  Mean = %.4f, SD = %.4f\n",
              mean(dt$value), stats::sd(dt$value)))
  cat(sprintf("  Range = [%.4f, %.4f]\n",
              meta$value_range[1L], meta$value_range[2L]))

  if (meta$has_group) {
    cat(sprintf("\nGroups: %d levels\n",
                length(levels(dt$group))))
    grp_summary <- dt[, .(n_ids = data.table::uniqueN(id),
                           mean_value = mean(value),
                           sd_value = stats::sd(value)),
                       by = group]
    cat("  Per-group summary:\n")
    print(grp_summary, row.names = FALSE)
  }

  if (meta$has_primary_trait) {
    pt_nm <- meta$primary_trait_name %||% "primary_trait"
    pt_vals <- unique(dt[, .(id, primary_trait)])$primary_trait
    cat(sprintf("\nPrimary trait (%s):\n", pt_nm))
    cat(sprintf("  Mean = %.4f, SD = %.4f, Range = [%.4f, %.4f]\n",
                mean(pt_vals), stats::sd(pt_vals),
                min(pt_vals), max(pt_vals)))
  }

  cat("\n==============================================\n")
  invisible(meta)
}

#' Plot method for fda_data (spaghetti plot)
#'
#' Produces a spaghetti plot of functional curves. If `ggplot2` is available,
#' uses it with a colourblind-safe viridis palette; otherwise falls back to
#' base graphics.
#'
#' @param x An `fda_data` object.
#' @param max_curves Maximum number of individual curves to draw (default 50).
#'   Set to `Inf` to plot all.
#' @param alpha Line transparency (default 0.4).
#' @param ... Additional arguments passed to the plotting function.
#' @return Invisibly returns the ggplot object (if ggplot2 is available) or
#'   `NULL`.
#' @export
plot.fda_data <- function(x, max_curves = 50L, alpha = 0.4, ...) {
  meta <- attr(x, "fda_meta")
  use_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

  # Optionally subsample curves for readability
  all_ids <- unique(x$id)
  if (length(all_ids) > max_curves && is.finite(max_curves)) {
    set.seed(1L)
    sampled_ids <- sample(all_ids, max_curves)
    plot_dt <- x[id %in% sampled_ids]
    subtitle <- sprintf("Showing %d of %d curves (random sample)",
                        max_curves, length(all_ids))
  } else {
    plot_dt <- data.table::copy(x)
    subtitle <- NULL
  }

  if (use_ggplot) {
    p <- ggplot2::ggplot(
      plot_dt,
      ggplot2::aes(x = time, y = value, group = id)
    )

    if (meta$has_group) {
      p <- p + ggplot2::geom_line(
        ggplot2::aes(colour = group), alpha = alpha, ...
      )
      if (requireNamespace("viridis", quietly = TRUE)) {
        p <- p + viridis::scale_colour_viridis(discrete = TRUE)
      }
    } else {
      p <- p + ggplot2::geom_line(alpha = alpha, colour = "steelblue", ...)
    }

    p <- p +
      ggplot2::labs(
        x = "Time",
        y = "Value",
        title = "Functional data: spaghetti plot",
        subtitle = subtitle,
        colour = "Group"
      ) +
      ggplot2::theme_minimal(base_size = 12)

    print(p)
    return(invisible(p))
  }

  # --- Base graphics fallback -------------------------------------------------
  if (meta$has_group) {
    grp_levels <- levels(plot_dt$group)
    n_grp <- length(grp_levels)
    pal <- grDevices::hcl.colors(n_grp, palette = "viridis")
    col_map <- stats::setNames(pal, grp_levels)
  }

  plot(range(plot_dt$time), range(plot_dt$value),
       type = "n", xlab = "Time", ylab = "Value",
       main = "Functional data: spaghetti plot")

  plot_ids <- unique(plot_dt$id)
  for (pid in plot_ids) {
    sub <- plot_dt[id == pid]
    col_val <- if (meta$has_group) {
      grDevices::adjustcolor(col_map[as.character(sub$group[1L])], alpha.f = alpha)
    } else {
      grDevices::adjustcolor("steelblue", alpha.f = alpha)
    }
    lines(sub$time, sub$value, col = col_val)
  }

  if (meta$has_group) {
    legend("topright", legend = grp_levels, col = pal, lty = 1L,
           bty = "n", cex = 0.8)
  }

  invisible(NULL)
}

#' Subset method for fda_data
#'
#' Subsets an `fda_data` object while preserving the class and updating
#' metadata. Accepts the same arguments as [data.table::data.table]
#' subsetting (row expressions, column selections, etc.).
#'
#' @param x An `fda_data` object.
#' @param ... Arguments passed to `[.data.table` (e.g., row filter
#'   expressions, column selections, `by` grouping).
#' @return An `fda_data` object (if essential columns remain) or a
#'   `data.table`.
#' @export
`[.fda_data` <- function(x, ...) {

  # Temporarily strip fda_data class so data.table's [ handles NSE correctly
  old_class <- class(x)
  old_meta  <- attr(x, "fda_meta")
  data.table::setattr(x, "class", c("data.table", "data.frame"))

  # Reconstruct the [.data.table call and evaluate in the parent frame so

  # that data.table NSE can resolve variables from the caller's environment.
  .fda_dt_x <- x
  sc <- sys.call()
  sc[[1L]] <- quote(`[`)
  sc[[2L]] <- quote(.fda_dt_x)
  pf <- parent.frame()
  pf$.fda_dt_x <- x
  result <- eval(sc, envir = pf)
  rm(".fda_dt_x", envir = pf)

  # Restore original object's class (setattr modifies in place)
  data.table::setattr(x, "class", old_class)
  data.table::setattr(x, "fda_meta", old_meta)

  # If the result is still a data.table with the essential columns, re-class it
  if (is.data.frame(result) &&
      all(c("id", "time", "value") %in% names(result))) {
    # Rebuild metadata
    meta <- list(
      primary_trait_name = old_meta$primary_trait_name,
      has_primary_trait  = "primary_trait" %in% names(result),
      has_group          = "group" %in% names(result),
      has_spatial        = all(c("spatial_row", "spatial_col") %in%
                                 names(result)),
      n_ids              = data.table::uniqueN(result$id),
      n_timepoints       = data.table::uniqueN(result$time),
      time_range         = range(result$time),
      value_range        = range(result$value),
      extra_cols         = old_meta$extra_cols
    )
    data.table::setattr(result, "fda_meta", meta)
    data.table::setattr(result, "class",
                        c("fda_data", "data.table", "data.frame"))
  }
  result
}

#' Coerce fda_data to data.table
#'
#' @param x An `fda_data` object.
#' @param keep.rownames Ignored (for compatibility).
#' @param ... Additional arguments (ignored).
#' @return A plain [data.table::data.table] without `fda_data` class or
#'   metadata.
#' @export
as.data.table.fda_data <- function(x, keep.rownames = FALSE, ...) {
  dt <- data.table::copy(x)
  data.table::setattr(dt, "fda_meta", NULL)
  data.table::setattr(dt, "class", c("data.table", "data.frame"))
  dt
}


# ==============================================================================
# fda_model -- fitted FDA model container
# ==============================================================================

#' Create a new fda_model object (internal constructor)
#'
#' Backend-agnostic container for fitted functional data analysis model results.
#' Intended to be called by fitting functions (e.g., `fit_fda()`) rather than
#' directly by users.
#'
#' @param fitted_curves A [data.table::data.table] with columns: `id`, `group`,
#'   `time`, `fitted`, and optionally `se`, `ci_lower`, `ci_upper`.
#' @param coefficient_function A list with components: `time`, `beta`, `se`,
#'   `ci_lower`, `ci_upper` -- the estimated coefficient function beta(t).
#' @param variance_components A [data.table::data.table] of estimated variance
#'   components.
#' @param predictions A [data.table::data.table] of predictions for the primary
#'   trait.
#' @param residuals A numeric vector or data.table of model residuals.
#' @param basis A list describing the basis (type, degree, knots, etc.).
#' @param data The original `fda_data` object used for fitting.
#' @param engine Character: `"asreml"` or `"bayesreml"`.
#' @param call The matched call from the fitting function.
#' @param extras A list for backend-specific extras (e.g., posterior draws,
#'   raw model object).
#'
#' @return An object of class `fda_model`.
#' @noRd
new_fda_model <- function(fitted_curves, coefficient_function,
                          variance_components, predictions,
                          residuals, basis, data, engine, call,
                          extras = list()) {

  # --- Minimal validation -----------------------------------------------------
  stopifnot(
    is.data.table(fitted_curves),
    is.list(coefficient_function),
    is.data.table(variance_components),
    is.data.table(predictions),
    is.character(engine) && length(engine) == 1L
  )

  structure(
    list(
      fitted_curves          = fitted_curves,
      coefficient_function   = coefficient_function,
      variance_components    = variance_components,
      predictions            = predictions,
      residuals              = residuals,
      basis                  = basis,
      data                   = data,
      engine                 = engine,
      call                   = call,
      extras                 = extras
    ),
    class = "fda_model"
  )
}

# ---- fda_model methods -------------------------------------------------------

#' Print method for fda_model
#'
#' @param x An `fda_model` object.
#' @param ... Additional arguments (ignored).
#' @return Invisible `x`.
#' @export
print.fda_model <- function(x, ...) {
  cat("-- fda_model --------------------------------\n")
  cat(sprintf("  Engine           : %s\n", x$engine))
  if (!is.null(x$basis)) {
    btype <- x$basis$type %||% "unknown"
    cat(sprintf("  Basis            : %s\n", btype))
  }
  cat(sprintf("  Fitted curves    : %d rows (%d ids)\n",
              nrow(x$fitted_curves),
              data.table::uniqueN(x$fitted_curves$id)))
  cat(sprintf("  Variance comps   : %d\n", nrow(x$variance_components)))
  cat(sprintf("  Predictions      : %d rows\n", nrow(x$predictions)))
  if (length(x$extras) > 0L) {
    cat(sprintf("  Extras           : %s\n",
                paste(names(x$extras), collapse = ", ")))
  }
  cat("---------------------------------------------\n")
  invisible(x)
}

#' Summary method for fda_model
#'
#' @param object An `fda_model` object.
#' @param ... Additional arguments (ignored).
#' @return Invisible list with summary components.
#' @export
summary.fda_model <- function(object, ...) {
  cat("== fda_model summary =========================\n\n")
  cat("Call:\n")
  print(object$call)
  cat(sprintf("\nEngine: %s\n", object$engine))

  # Coefficient function summary
  cf <- object$coefficient_function
  if (!is.null(cf) && length(cf$beta) > 0L) {
    cat("\nCoefficient function beta(t):\n")
    cat(sprintf("  Time range  : [%.2f, %.2f]\n",
                min(cf$time), max(cf$time)))
    cat(sprintf("  beta range  : [%.4f, %.4f]\n",
                min(cf$beta), max(cf$beta)))
    cat(sprintf("  Evaluations : %d time points\n", length(cf$time)))
  }

  # Variance components
  cat("\nVariance components:\n")
  print(object$variance_components)

  # Prediction summary
  if (nrow(object$predictions) > 0L) {
    cat(sprintf("\nPredictions: %d entries\n", nrow(object$predictions)))
  }

  cat("\n==============================================\n")
  invisible(list(
    engine     = object$engine,
    coef_fn    = object$coefficient_function,
    var_comps  = object$variance_components
  ))
}

#' Extract coefficient function beta(t) from fda_model
#'
#' @param object An `fda_model` object.
#' @param ... Additional arguments (ignored).
#' @return A [data.table::data.table] with columns `time`, `beta`, `se`,
#'   `ci_lower`, `ci_upper`.
#' @export
coef.fda_model <- function(object, ...) {
  cf <- object$coefficient_function
  data.table::data.table(
    time     = cf$time,
    beta     = cf$beta,
    se       = cf$se,
    ci_lower = cf$ci_lower,
    ci_upper = cf$ci_upper
  )
}

#' Extract fitted curves from fda_model
#'
#' @param object An `fda_model` object.
#' @param ... Additional arguments (ignored).
#' @return A [data.table::data.table] with fitted curve data.
#' @export
fitted.fda_model <- function(object, ...) {
  object$fitted_curves
}

#' Extract residuals from fda_model
#'
#' @param object An `fda_model` object.
#' @param ... Additional arguments (ignored).
#' @return Residuals (numeric vector or data.table).
#' @export
residuals.fda_model <- function(object, ...) {
  object$residuals
}

#' Predict method for fda_model
#'
#' @param object An `fda_model` object.
#' @param newdata Optional new data for prediction.
#' @param ... Additional arguments (passed to backend).
#' @return A [data.table::data.table] of predictions.
#' @export
predict.fda_model <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    return(object$predictions)
  }
  # Delegate to backend-specific prediction (to be implemented per engine)
  predict_fn <- object$extras$predict_function
  if (is.null(predict_fn)) {
    stop("Prediction with new data is not yet supported for this model.",
         call. = FALSE)
  }
  predict_fn(object, newdata, ...)
}

#' Extract variance-covariance information from fda_model
#'
#' @param object An `fda_model` object.
#' @param ... Additional arguments (ignored).
#' @return A [data.table::data.table] of variance components, or a matrix if
#'   available in `extras$vcov_matrix`.
#' @export
vcov.fda_model <- function(object, ...) {
  if (!is.null(object$extras$vcov_matrix)) {
    return(object$extras$vcov_matrix)
  }
  object$variance_components
}

#' Plot method for fda_model
#'
#' Produces diagnostic and result plots. If `ggplot2` is available, generates
#' a multi-panel figure showing: (1) the coefficient function beta(t) with
#' confidence band, and (2) fitted vs observed overlay. Falls back to base
#' graphics otherwise.
#'
#' @param x An `fda_model` object.
#' @param which Character: `"coef"` for coefficient function, `"fitted"` for
#'   fitted curves, `"both"` (default).
#' @param ... Additional arguments passed to plotting functions.
#' @return Invisible list of ggplot objects (if ggplot2 available), or NULL.
#' @export
plot.fda_model <- function(x, which = c("both", "coef", "fitted"), ...) {
  which <- match.arg(which)
  use_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
  plots <- list()

  # --- Coefficient function plot ----------------------------------------------
  if (which %in% c("both", "coef")) {
    cf_dt <- coef.fda_model(x)

    if (use_ggplot) {
      p_coef <- ggplot2::ggplot(cf_dt, ggplot2::aes(x = time, y = beta)) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
          fill = "steelblue", alpha = 0.2
        ) +
        ggplot2::geom_line(colour = "steelblue", linewidth = 1) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                            colour = "grey40") +
        ggplot2::labs(x = "Time", y = expression(beta(t)),
                      title = "Coefficient function") +
        ggplot2::theme_minimal(base_size = 12)
      plots$coef <- p_coef
      if (which == "coef") print(p_coef)
    } else {
      plot(cf_dt$time, cf_dt$beta, type = "l", col = "steelblue", lwd = 2,
           xlab = "Time", ylab = expression(beta(t)),
           main = "Coefficient function",
           ylim = range(c(cf_dt$ci_lower, cf_dt$ci_upper), na.rm = TRUE))
      if (!is.null(cf_dt$ci_lower)) {
        graphics::polygon(
          c(cf_dt$time, rev(cf_dt$time)),
          c(cf_dt$ci_lower, rev(cf_dt$ci_upper)),
          col = grDevices::adjustcolor("steelblue", alpha.f = 0.2),
          border = NA
        )
        lines(cf_dt$time, cf_dt$beta, col = "steelblue", lwd = 2)
      }
      abline(h = 0, lty = 2, col = "grey40")
    }
  }

  # --- Fitted curves plot -----------------------------------------------------
  if (which %in% c("both", "fitted")) {
    fc <- x$fitted_curves
    if (use_ggplot) {
      p_fit <- ggplot2::ggplot(fc, ggplot2::aes(x = time, group = id)) +
        ggplot2::geom_line(ggplot2::aes(y = fitted, colour = group),
                           alpha = 0.5) +
        ggplot2::labs(x = "Time", y = "Fitted value",
                      title = "Fitted functional curves",
                      colour = "Group") +
        ggplot2::theme_minimal(base_size = 12)
      if (requireNamespace("viridis", quietly = TRUE)) {
        p_fit <- p_fit + viridis::scale_colour_viridis(discrete = TRUE)
      }
      plots$fitted <- p_fit
      if (which == "fitted") print(p_fit)
    } else {
      plot(range(fc$time), range(fc$fitted, na.rm = TRUE),
           type = "n", xlab = "Time", ylab = "Fitted value",
           main = "Fitted functional curves")
      for (uid in unique(fc$id)) {
        sub <- fc[id == uid]
        lines(sub$time, sub$fitted,
              col = grDevices::adjustcolor("steelblue", alpha.f = 0.4))
      }
    }
  }

  # --- Combine if both -------------------------------------------------------
  if (which == "both" && use_ggplot && length(plots) == 2L) {
    print(plots$coef)
    grDevices::devAskNewPage(TRUE)
    print(plots$fitted)
    grDevices::devAskNewPage(FALSE)
  }

  invisible(plots)
}


# ==============================================================================
# fda_comparison -- comparison of multiple models
# ==============================================================================

#' Create a new fda_comparison object (internal constructor)
#'
#' Stores the results of comparing multiple FDA model fits (e.g., different
#' engines, basis specifications, or covariance structures) alongside
#' comparison metrics.
#'
#' @param models A named list of `fda_model` objects.
#' @param metrics A [data.table::data.table] of comparison metrics (e.g.,
#'   AIC, BIC, RMSE, correlation) with one row per model.
#' @param label Character string describing the comparison.
#'
#' @return An object of class `fda_comparison`.
#' @noRd
new_fda_comparison <- function(models, metrics, label) {
  stopifnot(
    is.list(models),
    all(vapply(models, inherits, logical(1L), "fda_model")),
    is.data.table(metrics),
    is.character(label) && length(label) == 1L
  )

  if (is.null(names(models))) {
    names(models) <- paste0("model_", seq_along(models))
  }

  structure(
    list(
      models  = models,
      metrics = metrics,
      label   = label
    ),
    class = "fda_comparison"
  )
}

# ---- fda_comparison methods --------------------------------------------------

#' Print method for fda_comparison
#'
#' @param x An `fda_comparison` object.
#' @param ... Additional arguments (ignored).
#' @return Invisible `x`.
#' @export
print.fda_comparison <- function(x, ...) {
  cat("-- fda_comparison ---------------------------\n")
  cat(sprintf("  Label  : %s\n", x$label))
  cat(sprintf("  Models : %d (%s)\n",
              length(x$models),
              paste(names(x$models), collapse = ", ")))
  cat("\nComparison metrics:\n")
  print(x$metrics)
  cat("---------------------------------------------\n")
  invisible(x)
}

#' Summary method for fda_comparison
#'
#' @param object An `fda_comparison` object.
#' @param ... Additional arguments (ignored).
#' @return Invisible metrics data.table.
#' @export
summary.fda_comparison <- function(object, ...) {
  cat("== fda_comparison summary ====================\n\n")
  cat(sprintf("Comparison: %s\n\n", object$label))

  # Per-model engine info
  engines <- vapply(object$models, function(m) m$engine, character(1L))
  cat("Models:\n")
  for (nm in names(object$models)) {
    cat(sprintf("  %s -- engine: %s\n", nm, engines[nm]))
  }

  cat("\nMetrics:\n")
  print(object$metrics)

  # Highlight best model for each numeric metric
  numeric_cols <- names(object$metrics)[
    vapply(object$metrics, is.numeric, logical(1L))
  ]
  if (length(numeric_cols) > 0L) {
    cat("\nBest model per metric (lowest value):\n")
    for (col in numeric_cols) {
      best_idx <- which.min(object$metrics[[col]])
      if (length(best_idx) == 1L) {
        cat(sprintf("  %-20s : %s (%.4f)\n",
                    col,
                    object$metrics$model[best_idx] %||%
                      names(object$models)[best_idx],
                    object$metrics[[col]][best_idx]))
      }
    }
  }

  cat("\n==============================================\n")
  invisible(object$metrics)
}

#' Plot method for fda_comparison
#'
#' Produces overlay plots of coefficient functions from competing models for
#' visual comparison, plus a metric comparison bar chart.
#'
#' @param x An `fda_comparison` object.
#' @param ... Additional arguments (ignored).
#' @return Invisible list of ggplot objects (if available), or NULL.
#' @export
plot.fda_comparison <- function(x, ...) {
  use_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
  plots <- list()

  # --- Overlay coefficient functions ------------------------------------------
  coef_list <- lapply(names(x$models), function(nm) {
    cf <- coef.fda_model(x$models[[nm]])
    cf[, model := nm]
    cf
  })
  coef_all <- data.table::rbindlist(coef_list)

  if (use_ggplot) {
    p_coef <- ggplot2::ggplot(
      coef_all,
      ggplot2::aes(x = time, y = beta, colour = model)
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper, fill = model),
        alpha = 0.1, colour = NA
      ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          colour = "grey40") +
      ggplot2::labs(
        x = "Time", y = expression(beta(t)),
        title = sprintf("Coefficient functions: %s", x$label),
        colour = "Model", fill = "Model"
      ) +
      ggplot2::theme_minimal(base_size = 12)

    if (requireNamespace("viridis", quietly = TRUE)) {
      p_coef <- p_coef +
        viridis::scale_colour_viridis(discrete = TRUE) +
        viridis::scale_fill_viridis(discrete = TRUE)
    }
    plots$coef <- p_coef
    print(p_coef)
  } else {
    # Base graphics overlay
    plot(range(coef_all$time), range(coef_all$beta, na.rm = TRUE),
         type = "n", xlab = "Time", ylab = expression(beta(t)),
         main = sprintf("Coefficient functions: %s", x$label))
    model_names <- unique(coef_all$model)
    pal <- grDevices::hcl.colors(length(model_names), palette = "viridis")
    for (i in seq_along(model_names)) {
      sub <- coef_all[model == model_names[i]]
      lines(sub$time, sub$beta, col = pal[i], lwd = 2)
    }
    legend("topright", legend = model_names, col = pal, lty = 1L,
           lwd = 2, bty = "n", cex = 0.8)
    abline(h = 0, lty = 2, col = "grey40")
  }

  invisible(plots)
}
