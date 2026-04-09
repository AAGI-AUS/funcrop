# plot.R -- Publication-ready visualisation suite for funcrop
#
# All plotting functions return ggplot objects (invisibly) when ggplot2 is
# available, with a base R fallback that prints a message recommending ggplot2
# for publication-quality output. Uses colourblind-safe viridis palettes
# throughout. Designed to handle both ASReml (REML) and bayesreml (Bayesian)
# model objects transparently.
#
# Author: Max Moldovan
# Licence: GPL (>= 3)


# ---- Internal helpers --------------------------------------------------------

#' Check ggplot2 availability and message if absent
#' @return Logical: TRUE if ggplot2 is available.
#' @noRd
.has_ggplot2 <- function() {
  requireNamespace("ggplot2", quietly = TRUE)
}

#' Check viridis availability
#' @return Logical: TRUE if viridis is available.
#' @noRd
.has_viridis <- function() {
  requireNamespace("viridis", quietly = TRUE)
}

#' Check ggrepel availability
#' @return Logical: TRUE if ggrepel is available.
#' @noRd
.has_ggrepel <- function() {
  requireNamespace("ggrepel", quietly = TRUE)
}

#' Apply consistent funcrop ggplot2 theme
#' @param base_size Numeric; base font size. Default 12.
#' @return A ggplot2 theme object.
#' @noRd
.funcrop_theme <- function(base_size = 12) {
  if (!.has_ggplot2()) return(NULL)
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle    = ggplot2::element_text(colour = "grey40"),
      axis.title       = ggplot2::element_text(face = "bold"),
      legend.position  = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold")
    )
}

#' Apply viridis discrete colour scale if available
#' @param ... Arguments passed to the scale function.
#' @return A ggplot2 scale layer or NULL.
#' @noRd
.scale_colour_viridis_d <- function(...) {
  if (.has_viridis()) {
    viridis::scale_colour_viridis(discrete = TRUE, option = "D", ...)
  } else if (.has_ggplot2()) {
    ggplot2::scale_colour_brewer(palette = "Dark2", ...)
  } else {
    NULL
  }
}

#' Apply viridis fill scale if available (continuous)
#' @param ... Arguments passed to the scale function.
#' @return A ggplot2 scale layer or NULL.
#' @noRd
.scale_fill_viridis_c <- function(...) {
  if (.has_viridis()) {
    viridis::scale_fill_viridis(option = "D", ...)
  } else if (.has_ggplot2()) {
    ggplot2::scale_fill_distiller(palette = "YlGnBu", ...)
  } else {
    NULL
  }
}

#' Apply viridis discrete fill scale if available
#' @param ... Arguments passed to the scale function.
#' @return A ggplot2 scale layer or NULL.
#' @noRd
.scale_fill_viridis_d <- function(...) {
  if (.has_viridis()) {
    viridis::scale_fill_viridis(discrete = TRUE, option = "D", ...)
  } else if (.has_ggplot2()) {
    ggplot2::scale_fill_brewer(palette = "Dark2", ...)
  } else {
    NULL
  }
}


# ==============================================================================
# plot_functional_profiles
# ==============================================================================

#' Plot Fitted Functional Profiles
#'
#' Produces publication-ready plots of variety-specific fitted functional curves
#' from an `fda_model` object. Optionally includes confidence/credible interval
#' ribbons. When ggplot2 is available, uses colourblind-safe viridis palettes
#' and a clean minimal theme; otherwise falls back to base R graphics.
#'
#' @param model An `fda_model` object (e.g., from [fit_functional_profiles()],
#'   [fit_fda_joint()], or [fit_2d_functional()]).
#' @param varieties Character vector of variety/group names to plot, or `NULL`
#'   (default) to plot all varieties. Useful for reducing clutter when many
#'   varieties are present.
#' @param ci Logical; if `TRUE` (default), add confidence/credible interval
#'   ribbons around each fitted curve. Requires `se` or `ci_lower`/`ci_upper`
#'   columns in the fitted curves data.
#' @param n_grid Integer; number of grid points for smooth curve evaluation.
#'   Default 200. Higher values give smoother curves at marginal cost.
#' @param alpha Numeric; ribbon transparency (0--1). Default 0.2.
#' @param line_size Numeric; line width for fitted curves. Default 0.8.
#' @param ... Additional arguments (currently unused).
#'
#' @return A `ggplot` object (invisibly), or `NULL` if ggplot2 is unavailable
#'   (base R plot is drawn directly).
#'
#' @examples
#' \dontrun{
#' model <- fit_functional_profiles(sim_grain_fill, ...)
#' plot_functional_profiles(model)
#' plot_functional_profiles(model, varieties = c("V01", "V05", "V10"))
#' }
#'
#' @export
plot_functional_profiles <- function(model,
                                     varieties = NULL,
                                     ci = TRUE,
                                     n_grid = 200L,
                                     alpha = 0.2,
                                     line_size = 0.8,
                                     ...) {

  if (!inherits(model, "fda_model")) {
    stop("`model` must be an 'fda_model' object.", call. = FALSE)
  }

  fc <- model$fitted_curves
  if (is.null(fc) || nrow(fc) == 0L) {
    stop("No fitted curves available in the model object.", call. = FALSE)
  }

  # Determine group column name
  group_var <- if ("group" %in% names(fc)) "group" else "id"

  # Subset varieties if requested
  if (!is.null(varieties)) {
    fc <- fc[fc[[group_var]] %in% varieties, ]
    if (nrow(fc) == 0L) {
      stop("No matching varieties found in fitted curves.", call. = FALSE)
    }
  }

  # Determine if this is a MET model (has environment facet)
  has_env <- "environment" %in% names(fc)

  # Check for CI columns
  has_ci <- ci && all(c("ci_lower", "ci_upper") %in% names(fc))
  has_se <- ci && !has_ci && "se" %in% names(fc)

  if (ci && has_se && !has_ci) {
    # Compute approximate CIs from SE (Wald-type 95%)
    fc <- data.table::copy(fc)
    fc[, ci_lower := fitted - 1.96 * se]
    fc[, ci_upper := fitted + 1.96 * se]
    has_ci <- TRUE
  }

  # ---- ggplot2 path -----------------------------------------------------------
  if (.has_ggplot2()) {
    p <- ggplot2::ggplot(
      fc,
      ggplot2::aes(
        x      = .data[["time"]],
        y      = .data[["fitted"]],
        colour = .data[[group_var]],
        fill   = .data[[group_var]],
        group  = .data[[group_var]]
      )
    )

    if (has_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data[["ci_lower"]],
          ymax = .data[["ci_upper"]]
        ),
        alpha  = alpha,
        colour = NA
      )
    }

    p <- p +
      ggplot2::geom_line(linewidth = line_size) +
      .scale_colour_viridis_d() +
      .scale_fill_viridis_d() +
      ggplot2::labs(
        x      = "Time",
        y      = "Fitted value",
        colour = "Variety",
        fill   = "Variety",
        title  = "Fitted functional profiles"
      ) +
      .funcrop_theme()

    if (has_env) {
      p <- p + ggplot2::facet_wrap(~ environment)
    }

    print(p)
    return(invisible(p))
  }

  # ---- Base R fallback --------------------------------------------------------
  message("Install ggplot2 for publication-quality plots. Using base R fallback.")

  groups <- unique(fc[[group_var]])
  n_grp  <- length(groups)
  pal    <- grDevices::hcl.colors(n_grp, palette = "viridis")

  plot(range(fc$time, na.rm = TRUE), range(fc$fitted, na.rm = TRUE),
       type = "n", xlab = "Time", ylab = "Fitted value",
       main = "Fitted functional profiles")

  for (i in seq_along(groups)) {
    sub <- fc[fc[[group_var]] == groups[i], ]
    sub <- sub[order(sub$time), ]

    if (has_ci) {
      graphics::polygon(
        c(sub$time, rev(sub$time)),
        c(sub$ci_lower, rev(sub$ci_upper)),
        col = grDevices::adjustcolor(pal[i], alpha.f = alpha),
        border = NA
      )
    }

    lines(sub$time, sub$fitted, col = pal[i], lwd = line_size * 2)
  }

  if (n_grp <= 20L) {
    legend("topright", legend = groups, col = pal, lty = 1L,
           lwd = 2, bty = "n", cex = 0.7)
  }

  invisible(NULL)
}


# ==============================================================================
# plot_coefficient_function
# ==============================================================================

#' Plot the Coefficient Function
#'
#' Plots the estimated coefficient function \eqn{\beta(t)} from a
#' scalar-on-function model (or \eqn{\beta(t, d)} heatmap/contour from a
#' scalar-on-2D-function model). Includes confidence/credible interval ribbons
#' and a horizontal reference line at \eqn{\beta = 0}.
#'
#' @param model An `fda_model` object from [scalar_on_function()] or
#'   [scalar_on_2d_function()].
#' @param ci Logical; if `TRUE` (default), add confidence/credible interval
#'   ribbons (1D) or a separate uncertainty panel (2D).
#' @param n_grid Integer; number of grid points for evaluation. Default 200
#'   (1D) or 50 per dimension (2D).
#' @param ... Additional arguments (currently unused).
#'
#' @return A `ggplot` object (invisibly), or `NULL` for base R fallback.
#'
#' @details
#' For 1D models, the plot shows \eqn{\beta(t)} as a line with a ribbon for
#' the CI/credible interval and a dashed horizontal line at zero. For 2D
#' models (from [scalar_on_2d_function()]), a heatmap of \eqn{\beta(t, d)}
#' is drawn using [ggplot2::geom_tile()] with a viridis continuous fill scale,
#' overlaid with contour lines.
#'
#' @examples
#' \dontrun{
#' sof_model <- scalar_on_function(yields, profiles)
#' plot_coefficient_function(sof_model)
#' }
#'
#' @export
plot_coefficient_function <- function(model,
                                      ci = TRUE,
                                      n_grid = 200L,
                                      ...) {

  if (!inherits(model, "fda_model")) {
    stop("`model` must be an 'fda_model' object.", call. = FALSE)
  }

  cf <- model$coefficient_function
  if (is.null(cf) || length(cf$beta) == 0L) {
    stop("No coefficient function available in the model.", call. = FALSE)
  }

  # Detect if this is a 2D model
  is_2d <- "depth" %in% names(cf)

  if (is_2d) {
    return(.plot_coef_2d(cf, ci = ci, ...))
  }

  # ---- 1D coefficient function -----------------------------------------------
  cf_dt <- data.table::data.table(
    time     = cf$time,
    beta     = cf$beta,
    se       = cf$se %||% rep(NA_real_, length(cf$beta)),
    ci_lower = cf$ci_lower %||% rep(NA_real_, length(cf$beta)),
    ci_upper = cf$ci_upper %||% rep(NA_real_, length(cf$beta))
  )

  # Compute CI from SE if needed
  has_ci <- !all(is.na(cf_dt$ci_lower))
  if (ci && !has_ci && !all(is.na(cf_dt$se))) {
    cf_dt[, ci_lower := beta - 1.96 * se]
    cf_dt[, ci_upper := beta + 1.96 * se]
    has_ci <- TRUE
  }

  if (.has_ggplot2()) {
    p <- ggplot2::ggplot(cf_dt, ggplot2::aes(x = time, y = beta))

    if (ci && has_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        fill = "#31688EFF", alpha = 0.2
      )
    }

    p <- p +
      ggplot2::geom_line(colour = "#31688EFF", linewidth = 1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          colour = "grey40", linewidth = 0.5) +
      ggplot2::labs(
        x     = "Time",
        y     = expression(beta(t)),
        title = "Coefficient function"
      ) +
      .funcrop_theme()

    # Add subtitle indicating engine
    if (!is.null(model$engine)) {
      p <- p + ggplot2::labs(
        subtitle = sprintf("Engine: %s", model$engine)
      )
    }

    print(p)
    return(invisible(p))
  }

  # ---- Base R fallback --------------------------------------------------------
  message("Install ggplot2 for publication-quality plots. Using base R fallback.")

  y_range <- range(c(cf_dt$beta, cf_dt$ci_lower, cf_dt$ci_upper),
                   na.rm = TRUE)

  plot(cf_dt$time, cf_dt$beta, type = "l", col = "#31688E", lwd = 2,
       xlab = "Time", ylab = expression(beta(t)),
       main = "Coefficient function", ylim = y_range)

  if (ci && has_ci) {
    graphics::polygon(
      c(cf_dt$time, rev(cf_dt$time)),
      c(cf_dt$ci_lower, rev(cf_dt$ci_upper)),
      col = grDevices::adjustcolor("#31688E", alpha.f = 0.2),
      border = NA
    )
    lines(cf_dt$time, cf_dt$beta, col = "#31688E", lwd = 2)
  }

  abline(h = 0, lty = 2, col = "grey40")
  invisible(NULL)
}


#' Internal: plot 2D coefficient surface
#' @param cf List with time, depth, beta, se, ci_lower, ci_upper.
#' @param ci Logical; show uncertainty.
#' @param ... Ignored.
#' @return ggplot object (invisible) or NULL.
#' @noRd
.plot_coef_2d <- function(cf, ci = TRUE, ...) {

  cf_dt <- data.table::data.table(
    time  = cf$time,
    depth = cf$depth,
    beta  = cf$beta,
    se    = cf$se %||% rep(NA_real_, length(cf$beta))
  )

  if (.has_ggplot2()) {
    p <- ggplot2::ggplot(cf_dt, ggplot2::aes(x = time, y = depth)) +
      ggplot2::geom_tile(ggplot2::aes(fill = beta)) +
      ggplot2::geom_contour(ggplot2::aes(z = beta),
                            colour = "white", alpha = 0.5,
                            linewidth = 0.3) +
      .scale_fill_viridis_c() +
      ggplot2::labs(
        x     = "Time",
        y     = "Depth",
        fill  = expression(beta(t, d)),
        title = "2D coefficient surface"
      ) +
      .funcrop_theme() +
      ggplot2::coord_fixed(ratio = diff(range(cf_dt$time)) /
                                   diff(range(cf_dt$depth)))

    print(p)
    return(invisible(p))
  }

  # Base R fallback: simple image plot
  message("Install ggplot2 for publication-quality 2D plots. Using base R fallback.")

  t_vals <- sort(unique(cf_dt$time))
  d_vals <- sort(unique(cf_dt$depth))
  z_mat  <- matrix(cf_dt$beta, nrow = length(t_vals), ncol = length(d_vals),
                   byrow = FALSE)

  pal <- grDevices::hcl.colors(256, palette = "viridis")
  graphics::image(t_vals, d_vals, z_mat, col = pal,
                  xlab = "Time", ylab = "Depth",
                  main = "2D coefficient surface")

  invisible(NULL)
}


# ==============================================================================
# plot_gxe_heatmap
# ==============================================================================

#' Heatmap of Genotype-by-Environment Interaction
#'
#' Visualises the variety x environment interaction for functional traits using
#' a heatmap. Can display variance component estimates, FA loadings, or
#' predicted functional summary statistics across environments.
#'
#' @param model An `fda_model` object from a multi-environment trial (MET)
#'   analysis.
#' @param value_col Character; which quantity to display. One of `"predicted"`
#'   (default, variety x environment predicted values), `"variance"` (variance
#'   components), or `"correlation"` (pairwise environment correlations).
#' @param ... Additional arguments (currently unused).
#'
#' @return A `ggplot` object (invisibly), or `NULL` for base R fallback.
#'
#' @examples
#' \dontrun{
#' met_model <- fit_functional_profiles(sim_met_fda, ...)
#' plot_gxe_heatmap(met_model)
#' }
#'
#' @export
plot_gxe_heatmap <- function(model,
                              value_col = c("predicted", "variance",
                                            "correlation"),
                              ...) {

  if (!inherits(model, "fda_model")) {
    stop("`model` must be an 'fda_model' object.", call. = FALSE)
  }

  value_col <- match.arg(value_col)

  # Attempt to extract variety x environment data
  # Look in predictions, fitted_curves, or extras
  hm_dt <- NULL

  if (value_col == "predicted") {
    preds <- model$predictions
    if (!is.null(preds) && nrow(preds) > 0L) {
      # Check for environment column
      if ("environment" %in% names(preds)) {
        hm_dt <- data.table::copy(preds)
        if (!"value" %in% names(hm_dt)) {
          hm_dt[, value := hm_dt[["predicted"]]]
        }
      } else if ("group" %in% names(preds)) {
        # Try fitted_curves for environment info
        fc <- model$fitted_curves
        if (!is.null(fc) && "environment" %in% names(fc)) {
          # Summarise fitted curves: mean fitted value per group x environment
          hm_dt <- fc[, .(value = mean(fitted, na.rm = TRUE)),
                       by = c("id", "environment")]
          data.table::setnames(hm_dt, "id", "group")
        }
      }
    }

    if (is.null(hm_dt) || nrow(hm_dt) == 0L) {
      # Fallback: try to compute from fitted_curves
      fc <- model$fitted_curves
      if (!is.null(fc) && "environment" %in% names(fc)) {
        group_var <- if ("group" %in% names(fc)) "group" else "id"
        hm_dt <- fc[, .(value = mean(fitted, na.rm = TRUE)),
                     by = c(group_var, "environment")]
        if (group_var != "group") {
          data.table::setnames(hm_dt, group_var, "group")
        }
      }
    }
  } else if (value_col == "variance") {
    vc <- model$variance_components
    if (!is.null(vc) && nrow(vc) > 0L) {
      hm_dt <- data.table::copy(vc)
      hm_dt[, value := estimate]
    }
  } else {
    # correlation: extract or compute environment correlation matrix
    if (!is.null(model$extras$env_correlation)) {
      cor_mat <- model$extras$env_correlation
      envs <- rownames(cor_mat)
      hm_dt <- data.table::CJ(env1 = envs, env2 = envs)
      hm_dt[, value := as.numeric(cor_mat)]
    }
  }

  if (is.null(hm_dt) || nrow(hm_dt) == 0L) {
    stop(
      "Cannot extract variety x environment data for heatmap. ",
      "Ensure the model is from a multi-environment analysis.",
      call. = FALSE
    )
  }

  # ---- ggplot2 path -----------------------------------------------------------
  if (.has_ggplot2()) {

    # Determine axes based on available columns
    if (all(c("group", "environment") %in% names(hm_dt))) {
      p <- ggplot2::ggplot(
        hm_dt,
        ggplot2::aes(x = environment, y = group, fill = value)
      ) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
        .scale_fill_viridis_c() +
        ggplot2::labs(
          x     = "Environment",
          y     = "Variety",
          fill  = "Value",
          title = "Genotype x Environment interaction"
        ) +
        .funcrop_theme() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    } else if (all(c("env1", "env2") %in% names(hm_dt))) {
      # Correlation heatmap
      p <- ggplot2::ggplot(
        hm_dt,
        ggplot2::aes(x = env1, y = env2, fill = value)
      ) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)),
                           size = 3, colour = "white") +
        .scale_fill_viridis_c() +
        ggplot2::labs(
          x     = "Environment",
          y     = "Environment",
          fill  = "Correlation",
          title = "Environment correlation matrix"
        ) +
        .funcrop_theme() +
        ggplot2::coord_fixed()
    } else {
      stop(
        "Unexpected data structure for GxE heatmap. ",
        "Expected 'group'+'environment' or 'env1'+'env2' columns.",
        call. = FALSE
      )
    }

    print(p)
    return(invisible(p))
  }

  # ---- Base R fallback --------------------------------------------------------
  message("Install ggplot2 for publication-quality heatmaps. Using base R fallback.")

  if (all(c("group", "environment") %in% names(hm_dt))) {
    groups <- unique(hm_dt$group)
    envs   <- unique(hm_dt$environment)
    z_mat  <- matrix(NA_real_, nrow = length(groups), ncol = length(envs))
    for (i in seq_along(groups)) {
      for (j in seq_along(envs)) {
        val <- hm_dt[group == groups[i] & environment == envs[j], value]
        if (length(val) > 0L) z_mat[i, j] <- val[1L]
      }
    }
    pal <- grDevices::hcl.colors(256, palette = "viridis")
    graphics::image(seq_along(envs), seq_along(groups), t(z_mat),
                    col = pal, axes = FALSE,
                    xlab = "Environment", ylab = "Variety",
                    main = "GxE interaction heatmap")
    graphics::axis(1, at = seq_along(envs), labels = envs, las = 2, cex.axis = 0.7)
    graphics::axis(2, at = seq_along(groups), labels = groups, las = 1, cex.axis = 0.5)
  }

  invisible(NULL)
}


# ==============================================================================
# plot_fa_biplot
# ==============================================================================

#' Factor-Analytic Biplot for GxE Interaction
#'
#' Produces a biplot from the factor-analytic (FA) decomposition of
#' genotype-by-environment interaction effects estimated from B-spline
#' coefficients. Environment loadings and variety scores are overlaid in a
#' two-dimensional space showing the dominant patterns of GxE interaction.
#'
#' @param model An `fda_model` object containing FA structure in extras
#'   (e.g., `model$extras$fa_loadings` and `model$extras$fa_scores`), or
#'   with variance components that can be decomposed.
#' @param factors Integer vector of length 2 specifying which FA factors to
#'   plot. Default `c(1, 2)`.
#' @param point_size Numeric; size of points. Default 3.
#' @param text_size Numeric; size of labels. Default 3.5.
#' @param ... Additional arguments (currently unused).
#'
#' @return A `ggplot` object (invisibly), or `NULL` for base R fallback.
#'
#' @details
#' The biplot displays:
#' - **Environment loadings** as points/arrows, showing how environments
#'   differentiate varieties.
#' - **Variety scores** as labelled points, showing relative variety performance
#'   across the FA factors.
#'
#' If the `ggrepel` package is available, labels are automatically positioned
#' to avoid overlapping. Otherwise, [ggplot2::geom_text()] is used.
#'
#' @examples
#' \dontrun{
#' met_model <- fit_functional_profiles(sim_met_fda, ...)
#' plot_fa_biplot(met_model)
#' }
#'
#' @export
plot_fa_biplot <- function(model,
                           factors = c(1L, 2L),
                           point_size = 3,
                           text_size = 3.5,
                           ...) {

  if (!inherits(model, "fda_model")) {
    stop("`model` must be an 'fda_model' object.", call. = FALSE)
  }

  # ---- Extract FA structure ---------------------------------------------------

  fa_loadings <- model$extras$fa_loadings
  fa_scores   <- model$extras$fa_scores

  if (is.null(fa_loadings) || is.null(fa_scores)) {
    # Attempt to extract from variance components or raw model
    fa_result <- .extract_fa_structure(model)
    if (is.null(fa_result)) {
      stop(
        "No factor-analytic structure found in the model. ",
        "Ensure the model was fitted with an FA covariance structure for GxE.",
        call. = FALSE
      )
    }
    fa_loadings <- fa_result$loadings
    fa_scores   <- fa_result$scores
  }

  # Validate factor selection
  n_factors <- ncol(fa_loadings)
  if (length(factors) != 2L || any(factors > n_factors)) {
    stop(sprintf(
      "Model has %d FA factors. `factors` must be a length-2 integer vector ",
      "with values in 1:%d.", n_factors, n_factors
    ), call. = FALSE)
  }

  # Prepare data for plotting
  f1 <- factors[1L]
  f2 <- factors[2L]

  env_names <- rownames(fa_loadings)
  if (is.null(env_names)) env_names <- paste0("Env_", seq_len(nrow(fa_loadings)))

  var_names <- rownames(fa_scores)
  if (is.null(var_names)) var_names <- paste0("Var_", seq_len(nrow(fa_scores)))

  load_dt <- data.table::data.table(
    label = env_names,
    x     = fa_loadings[, f1],
    y     = fa_loadings[, f2],
    type  = "Environment"
  )

  score_dt <- data.table::data.table(
    label = var_names,
    x     = fa_scores[, f1],
    y     = fa_scores[, f2],
    type  = "Variety"
  )

  plot_dt <- data.table::rbindlist(list(load_dt, score_dt))

  # ---- ggplot2 path -----------------------------------------------------------
  if (.has_ggplot2()) {

    p <- ggplot2::ggplot(
      plot_dt,
      ggplot2::aes(x = x, y = y, colour = type, shape = type)
    ) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          colour = "grey60", linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                          colour = "grey60", linewidth = 0.3)

    # Add labels with ggrepel if available, else geom_text
    if (.has_ggrepel()) {
      p <- p + ggrepel::geom_text_repel(
        ggplot2::aes(label = label),
        size          = text_size,
        max.overlaps  = 20,
        show.legend   = FALSE
      )
    } else {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = label),
        size    = text_size,
        vjust   = -0.5,
        show.legend = FALSE
      )
    }

    # Add arrows from origin to environment loadings
    p <- p + ggplot2::geom_segment(
      data = load_dt,
      ggplot2::aes(x = 0, y = 0, xend = x, yend = y),
      arrow     = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
      linewidth = 0.4,
      alpha     = 0.5,
      colour    = "#440154FF",
      inherit.aes = FALSE
    )

    p <- p +
      ggplot2::scale_shape_manual(values = c("Environment" = 17,
                                              "Variety" = 16)) +
      .scale_colour_viridis_d() +
      ggplot2::labs(
        x      = sprintf("Factor %d", f1),
        y      = sprintf("Factor %d", f2),
        colour = "",
        shape  = "",
        title  = "Factor-analytic biplot",
        subtitle = sprintf("FA factors %d vs %d", f1, f2)
      ) +
      .funcrop_theme() +
      ggplot2::coord_fixed()

    print(p)
    return(invisible(p))
  }

  # ---- Base R fallback --------------------------------------------------------
  message("Install ggplot2 for publication-quality biplots. Using base R fallback.")

  x_range <- range(plot_dt$x, 0)
  y_range <- range(plot_dt$y, 0)

  plot(x_range, y_range, type = "n",
       xlab = sprintf("Factor %d", f1),
       ylab = sprintf("Factor %d", f2),
       main = "FA biplot", asp = 1)
  abline(h = 0, v = 0, lty = 2, col = "grey60")

  # Environments: triangles + arrows
  points(load_dt$x, load_dt$y, pch = 17, col = "#440154", cex = 1.5)
  graphics::arrows(0, 0, load_dt$x * 0.95, load_dt$y * 0.95,
                   col = "#440154", length = 0.1)
  text(load_dt$x, load_dt$y, labels = load_dt$label,
       pos = 3, cex = 0.7, col = "#440154")

 # Varieties: circles
  points(score_dt$x, score_dt$y, pch = 16, col = "#31688E", cex = 1.2)
  text(score_dt$x, score_dt$y, labels = score_dt$label,
       pos = 3, cex = 0.6, col = "#31688E")

  legend("topright",
         legend = c("Environment", "Variety"),
         pch    = c(17, 16),
         col    = c("#440154", "#31688E"),
         bty    = "n")

  invisible(NULL)
}


#' Extract FA structure from an fda_model (internal)
#'
#' Attempts to extract factor-analytic loadings and scores from the model's
#' variance components or raw backend model object.
#'
#' @param model An `fda_model` object.
#' @return List with `loadings` and `scores` matrices, or NULL if unavailable.
#' @noRd
.extract_fa_structure <- function(model) {
  raw_model <- model$extras$raw_model
  if (is.null(raw_model)) return(NULL)

  if (model$engine == "asreml" && .has_asreml()) {
    # Try to extract FA parameters from ASReml model
    fa_terms <- tryCatch({
      vc <- summary(raw_model)$varcomp
      fa_rows <- grep("^fa\\(|!fa", rownames(vc))
      if (length(fa_rows) == 0L) return(NULL)
      vc[fa_rows, ]
    }, error = function(e) NULL)

    if (is.null(fa_terms)) return(NULL)

    # Parse FA loadings and specific variances from ASReml output
    # This is backend-specific and depends on the naming convention
    loading_rows <- grep("!fa.*\\.load", rownames(fa_terms))
    if (length(loading_rows) > 0L) {
      loadings <- matrix(fa_terms[loading_rows, "component"],
                        ncol = 1L)
      rownames(loadings) <- sub(".*!fa.*:", "",
                                 rownames(fa_terms)[loading_rows])
      # Scores would need to be extracted from BLUPs
      return(list(loadings = loadings, scores = matrix(0, 1, 1)))
    }
  }

  if (model$engine == "bayesreml" && .has_bayesreml()) {
    # bayesreml v0.1.0: FA structure stored in fit$extras
    tryCatch({
      fa_post <- raw_model$extras$fa_structure
      if (!is.null(fa_post)) {
        return(list(
          loadings = fa_post$loadings,
          scores   = fa_post$scores
        ))
      }
    }, error = function(e) NULL)
  }

  NULL
}


# ==============================================================================
# plot_backend_comparison
# ==============================================================================

#' Compare REML and Bayesian Model Fits
#'
#' Produces a two-panel comparison of results from ASReml-R (REML) and
#' bayesreml (Bayesian) backends. The left panel compares fitted functional
#' curves; the right panel compares coefficient function estimates. REML
#' results are shown as point estimates with standard error bands, while
#' Bayesian results show the posterior mean with credible intervals.
#'
#' @param model_reml An `fda_model` object fitted with the ASReml backend
#'   (`engine = "asreml"`).
#' @param model_bayes An `fda_model` object fitted with the bayesreml backend
#'   (`engine = "bayesreml"`).
#' @param varieties Character vector of variety names to include in the
#'   curves panel, or `NULL` (default) for all. Plots may be cluttered if
#'   many varieties are included; consider limiting to 5--10.
#' @param ... Additional arguments (currently unused).
#'
#' @return A `ggplot` object (invisibly), or `NULL` for base R fallback.
#'   When ggplot2 is available, both panels are arranged side by side.
#'
#' @details
#' Colour coding:
#' - ASReml (REML): solid lines with dashed CI bands
#' - bayesreml (Bayesian): semi-transparent ribbons for credible intervals
#'
#' @examples
#' \dontrun{
#' model_a <- fit_functional_profiles(..., engine = "asreml")
#' model_b <- fit_functional_profiles(..., engine = "bayesreml")
#' plot_backend_comparison(model_a, model_b)
#' }
#'
#' @export
plot_backend_comparison <- function(model_reml,
                                    model_bayes,
                                    varieties = NULL,
                                    ...) {

  if (!inherits(model_reml, "fda_model")) {
    stop("`model_reml` must be an 'fda_model' object.", call. = FALSE)
  }
  if (!inherits(model_bayes, "fda_model")) {
    stop("`model_bayes` must be an 'fda_model' object.", call. = FALSE)
  }

  # Colour palette for the two backends
  col_reml  <- "#440154FF"   # viridis dark purple
  col_bayes <- "#35B779FF"   # viridis green

  # ---- Prepare fitted curves data ---------------------------------------------

  fc_reml <- model_reml$fitted_curves
  fc_bayes <- model_bayes$fitted_curves

  if (!is.null(fc_reml) && nrow(fc_reml) > 0L) {
    fc_reml <- data.table::copy(fc_reml)
    fc_reml[, backend := "ASReml (REML)"]
  }
  if (!is.null(fc_bayes) && nrow(fc_bayes) > 0L) {
    fc_bayes <- data.table::copy(fc_bayes)
    fc_bayes[, backend := "bayesreml (Bayesian)"]
  }

  # Ensure matching columns
  group_var <- if ("group" %in% names(fc_reml)) "group" else "id"

  if (!is.null(varieties)) {
    if (!is.null(fc_reml)) {
      fc_reml <- fc_reml[fc_reml[[group_var]] %in% varieties, ]
    }
    if (!is.null(fc_bayes)) {
      fc_bayes <- fc_bayes[fc_bayes[[group_var]] %in% varieties, ]
    }
  }

  # ---- Prepare coefficient function data --------------------------------------

  cf_reml_dt <- NULL
  cf_bayes_dt <- NULL

  cf_r <- model_reml$coefficient_function
  if (!is.null(cf_r) && length(cf_r$beta) > 0L) {
    cf_reml_dt <- data.table::data.table(
      time     = cf_r$time,
      beta     = cf_r$beta,
      ci_lower = cf_r$ci_lower %||% (cf_r$beta - 1.96 * (cf_r$se %||% 0)),
      ci_upper = cf_r$ci_upper %||% (cf_r$beta + 1.96 * (cf_r$se %||% 0)),
      backend  = "ASReml (REML)"
    )
  }

  cf_b <- model_bayes$coefficient_function
  if (!is.null(cf_b) && length(cf_b$beta) > 0L) {
    cf_bayes_dt <- data.table::data.table(
      time     = cf_b$time,
      beta     = cf_b$beta,
      ci_lower = cf_b$ci_lower %||% (cf_b$beta - 1.96 * (cf_b$se %||% 0)),
      ci_upper = cf_b$ci_upper %||% (cf_b$beta + 1.96 * (cf_b$se %||% 0)),
      backend  = "bayesreml (Bayesian)"
    )
  }

  # ---- ggplot2 path -----------------------------------------------------------

  if (.has_ggplot2()) {
    plots <- list()

    # Panel 1: Fitted curves comparison
    if (!is.null(fc_reml) && !is.null(fc_bayes) &&
        nrow(fc_reml) > 0L && nrow(fc_bayes) > 0L) {

      # Ensure consistent columns
      keep_cols <- intersect(names(fc_reml), names(fc_bayes))
      keep_cols <- intersect(keep_cols,
                              c(group_var, "time", "fitted", "ci_lower",
                                "ci_upper", "se", "backend"))

      # Add CI columns if missing
      for (dt_ref in list(fc_reml, fc_bayes)) {
        if (!"ci_lower" %in% names(dt_ref) && "se" %in% names(dt_ref)) {
          dt_ref[, ci_lower := fitted - 1.96 * se]
          dt_ref[, ci_upper := fitted + 1.96 * se]
        }
      }

      fc_combined <- data.table::rbindlist(
        list(fc_reml, fc_bayes),
        use.names = TRUE, fill = TRUE
      )

      p1 <- ggplot2::ggplot(
        fc_combined,
        ggplot2::aes(
          x         = time,
          y         = fitted,
          colour    = backend,
          fill      = backend,
          group     = interaction(.data[[group_var]], backend),
          linetype  = backend
        )
      )

      if (all(c("ci_lower", "ci_upper") %in% names(fc_combined))) {
        p1 <- p1 + ggplot2::geom_ribbon(
          ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
          alpha = 0.1, colour = NA
        )
      }

      p1 <- p1 +
        ggplot2::geom_line(linewidth = 0.6, alpha = 0.7) +
        ggplot2::scale_colour_manual(
          values = c("ASReml (REML)" = col_reml,
                     "bayesreml (Bayesian)" = col_bayes)
        ) +
        ggplot2::scale_fill_manual(
          values = c("ASReml (REML)" = col_reml,
                     "bayesreml (Bayesian)" = col_bayes)
        ) +
        ggplot2::scale_linetype_manual(
          values = c("ASReml (REML)" = "solid",
                     "bayesreml (Bayesian)" = "longdash")
        ) +
        ggplot2::labs(
          x        = "Time",
          y        = "Fitted value",
          colour   = "Backend",
          fill     = "Backend",
          linetype = "Backend",
          title    = "Fitted curves"
        ) +
        .funcrop_theme()

      plots$curves <- p1
    }

    # Panel 2: Coefficient function comparison
    if (!is.null(cf_reml_dt) && !is.null(cf_bayes_dt)) {
      cf_combined <- data.table::rbindlist(
        list(cf_reml_dt, cf_bayes_dt),
        use.names = TRUE
      )

      p2 <- ggplot2::ggplot(
        cf_combined,
        ggplot2::aes(x = time, y = beta, colour = backend, fill = backend)
      ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
          alpha = 0.15, colour = NA
        ) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                            colour = "grey40", linewidth = 0.4) +
        ggplot2::scale_colour_manual(
          values = c("ASReml (REML)" = col_reml,
                     "bayesreml (Bayesian)" = col_bayes)
        ) +
        ggplot2::scale_fill_manual(
          values = c("ASReml (REML)" = col_reml,
                     "bayesreml (Bayesian)" = col_bayes)
        ) +
        ggplot2::labs(
          x      = "Time",
          y      = expression(beta(t)),
          colour = "Backend",
          fill   = "Backend",
          title  = "Coefficient function"
        ) +
        .funcrop_theme()

      plots$coef <- p2
    }

    # Display panels
    if (length(plots) == 2L) {
      # Print side by side if possible (suggest patchwork, but print separately)
      if (requireNamespace("patchwork", quietly = TRUE)) {
        combined <- plots$curves + plots$coef +
          patchwork::plot_layout(guides = "collect") &
          .funcrop_theme()
        print(combined)
        return(invisible(combined))
      } else {
        # Print sequentially
        print(plots$curves)
        grDevices::devAskNewPage(TRUE)
        print(plots$coef)
        grDevices::devAskNewPage(FALSE)
        return(invisible(plots))
      }
    } else if (length(plots) == 1L) {
      print(plots[[1L]])
      return(invisible(plots[[1L]]))
    }

    return(invisible(plots))
  }

  # ---- Base R fallback --------------------------------------------------------
  message("Install ggplot2 for publication-quality comparison plots. Using base R fallback.")

  old_par <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(old_par), add = TRUE)

  # Panel 1: Fitted curves
  if (!is.null(fc_reml) && !is.null(fc_bayes) &&
      nrow(fc_reml) > 0L && nrow(fc_bayes) > 0L) {
    y_range <- range(c(fc_reml$fitted, fc_bayes$fitted), na.rm = TRUE)
    x_range <- range(c(fc_reml$time, fc_bayes$time), na.rm = TRUE)

    plot(x_range, y_range, type = "n",
         xlab = "Time", ylab = "Fitted value",
         main = "Fitted curves comparison")

    for (gid in unique(fc_reml[[group_var]])) {
      sub_r <- fc_reml[fc_reml[[group_var]] == gid, ]
      sub_r <- sub_r[order(sub_r$time), ]
      lines(sub_r$time, sub_r$fitted,
            col = grDevices::adjustcolor(col_reml, 0.5), lwd = 1.5)
    }
    for (gid in unique(fc_bayes[[group_var]])) {
      sub_b <- fc_bayes[fc_bayes[[group_var]] == gid, ]
      sub_b <- sub_b[order(sub_b$time), ]
      lines(sub_b$time, sub_b$fitted,
            col = grDevices::adjustcolor(col_bayes, 0.5), lwd = 1.5, lty = 2)
    }
    legend("topright",
           legend = c("ASReml", "bayesreml"),
           col = c(col_reml, col_bayes), lty = c(1, 2), lwd = 2,
           bty = "n", cex = 0.8)
  }

  # Panel 2: Coefficient function
  if (!is.null(cf_reml_dt) && !is.null(cf_bayes_dt)) {
    y_range <- range(c(cf_reml_dt$ci_lower, cf_reml_dt$ci_upper,
                       cf_bayes_dt$ci_lower, cf_bayes_dt$ci_upper),
                     na.rm = TRUE)

    plot(range(cf_reml_dt$time), y_range, type = "n",
         xlab = "Time", ylab = expression(beta(t)),
         main = "Coefficient function comparison")
    abline(h = 0, lty = 2, col = "grey40")

    # REML
    graphics::polygon(
      c(cf_reml_dt$time, rev(cf_reml_dt$time)),
      c(cf_reml_dt$ci_lower, rev(cf_reml_dt$ci_upper)),
      col = grDevices::adjustcolor(col_reml, 0.15), border = NA
    )
    lines(cf_reml_dt$time, cf_reml_dt$beta, col = col_reml, lwd = 2)

    # Bayesian
    graphics::polygon(
      c(cf_bayes_dt$time, rev(cf_bayes_dt$time)),
      c(cf_bayes_dt$ci_lower, rev(cf_bayes_dt$ci_upper)),
      col = grDevices::adjustcolor(col_bayes, 0.15), border = NA
    )
    lines(cf_bayes_dt$time, cf_bayes_dt$beta, col = col_bayes, lwd = 2, lty = 2)

    legend("topright",
           legend = c("ASReml", "bayesreml"),
           col = c(col_reml, col_bayes), lty = c(1, 2), lwd = 2,
           bty = "n", cex = 0.8)
  }

  invisible(NULL)
}
