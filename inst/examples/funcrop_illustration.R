###############################################################################
#
#  funcrop: Full Illustration -- ASReml (REML) vs bayesreml (Bayesian)
#
#  Progresses through 7 models (simple -> complex) showing classical and
#  Bayesian approaches side-by-side. Each model block documents:
#    - The statistical model (fixed, random, covariance structures)
#    - What FDA adds beyond standard approaches
#    - Biological/agronomic interpretation in the grain-fill context
#
#  Author: Maksym Gaidashenko
#  Date:   2026-04-04
#  Licence: GPL (>= 3)
#
###############################################################################

# ---- 0. Setup ----------------------------------------------------------------

library(funcrop)
library(data.table)
library(Matrix)
library(ggplot2)
library(viridis)

# Conditional loads -- script runs even if one backend is missing
HAS_ASREML   <- requireNamespace("asreml", quietly = TRUE)
HAS_BAYESREML <- requireNamespace("bayesreml", quietly = TRUE)
if (HAS_ASREML)   library(asreml)
if (HAS_BAYESREML) library(bayesreml)

# Simulated grain-fill trial shipped with funcrop
data(sim_grain_fill)
dt <- copy(sim_grain_fill)  # 480 rows: 20 varieties x 3 blocks x 8 times

# Ensure factors
dt[, variety := factor(variety)]
dt[, block   := factor(block)]
dt[, row     := as.integer(row)]
dt[, col     := as.integer(col)]
dt[, plot_id := factor(plot_id)]

cat("\n===== Dataset overview =====\n")
cat("Rows:", nrow(dt), " | Varieties:", uniqueN(dt$variety),
    " | Blocks:", uniqueN(dt$block), " | Time points:", uniqueN(dt$time), "\n")
cat("Time points:", sort(unique(dt$time)), "\n")
cat("Grid: rows 1-", max(dt$row), " x cols 1-", max(dt$col), "\n\n")

# Yield in sim_grain_fill is per-variety (identical for all plots of a variety).
# Add plot-level noise so RCBD models are estimable.
set.seed(42)
plot_ids <- unique(dt$plot_id)
yield_noise <- data.table(
  plot_id = plot_ids,
  yield_noise = rnorm(length(plot_ids), 0, 0.8)
)
dt <- merge(dt, yield_noise, by = "plot_id", sort = FALSE)
dt[, yield_plot := yield + yield_noise]  # plot-level yield with noise

# Variety-mean yields (for scalar-on-function)
yield_dt <- dt[, .(yield = mean(yield_plot)), by = variety]
setkey(yield_dt, variety)

# Storage for cross-model comparison
results <- list()


###############################################################################
#  MODEL 0 -- Baseline: RCBD Spatial Model for Yield (No FDA)
###############################################################################
#
#  STATISTICAL MODEL:
#    yield_ij = mu + variety_i + block_j + epsilon_ij
#
#    - Fixed:  mu (intercept), block
#    - Random: variety ~ N(0, sigma_v^2)
#    - Residual: epsilon ~ N(0, sigma_e^2)
#
#  NO FUNCTIONAL DATA here -- this is the standard approach breeders use.
#  Yield is analysed as a single scalar observation per plot.
#
#  BIOLOGICAL INTERPRETATION:
#    Estimates variety means for yield. Block accounts for field gradient.
#    This is the benchmark: can FDA of the secondary trait (grain-fill)
#    improve upon this simple analysis?
#
###############################################################################

cat(strrep("=", 70), "\nMODEL 0: Baseline RCBD for yield (no FDA)\n", strrep("=", 70), "\n")

# One yield per plot (with plot-level noise for model estimability)
yield_plot_dt <- unique(dt[, .(plot_id, variety, block, row, col, yield_plot)])

# ---- 0a. ASReml ----
if (HAS_ASREML) {
  cat("\n--- ASReml ---\n")
  m0_asr <- asreml(
    fixed    = yield_plot ~ block,
    random   = ~ variety,
    data     = yield_plot_dt,
    trace    = FALSE
  )
  cat("Converged:", m0_asr$converge, "\n")

  # Variance components
  vc0_asr <- summary(m0_asr)$varcomp
  cat("\nVariance components:\n")
  print(vc0_asr)

  # Variety BLUPs
  blup0_asr <- predict(m0_asr, classify = "variety")$pvals
  setDT(blup0_asr)
  cat("\nTop 5 varieties (BLUPs):\n")
  print(head(blup0_asr[order(-predicted.value)], 5))

  results$m0_asr <- list(
    vc     = vc0_asr,
    blups  = blup0_asr,
    loglik = m0_asr$loglik
  )
}

# ---- 0b. bayesreml ----
if (HAS_BAYESREML) {
  cat("\n--- bayesreml ---\n")
  m0_bay <- bayesreml(
    fixed     = yield_plot ~ block,
    random    = ~ variety,
    data      = yield_plot_dt,
    n_samples = 2000,
    warmup    = 1000,
    chains    = 4,
    verbose   = FALSE,
    mcmc_verbose = FALSE
  )

  cat("\nVariance components (posterior):\n")
  print(m0_bay$extras$variance_comps)

  # IMPORTANT: bayesreml uses non-centred parameterisation internally:
  #   u_variety = variety_raw * sigma_variety
  # Only hyperparameters (sigma_variety, sigma_e) and fixed effects are in
  # the MCMC draws. Variety random effects (BLUPs) must be computed post-hoc
  # via greta::calculate() on the u_variety greta array.
  blup0_bay <- tryCatch({
    # Attempt to use greta::calculate for posterior variety effects
    if (!is.null(m0_bay$greta$greta_arrays$u_variety) &&
        !is.null(m0_bay$greta$draws)) {
      u_post <- greta::calculate(m0_bay$greta$greta_arrays$u_variety,
                                  values = m0_bay$greta$draws)
      u_means <- colMeans(do.call(rbind, u_post))
      data.table(variety = levels(yield_plot_dt$variety),
                 predicted.value = as.numeric(u_means))
    } else {
      cat("  (Variety BLUPs not available in draws -- need greta::calculate)\n")
      NULL
    }
  }, error = function(e) {
    cat("  (Could not extract BLUPs:", e$message, ")\n")
    NULL
  })

  if (!is.null(blup0_bay) && nrow(blup0_bay) > 0) {
    cat("\nTop 5 varieties (posterior means):\n")
    print(head(blup0_bay[order(-predicted.value)], 5))
  } else {
    cat("\n  Variety BLUPs: not directly in MCMC draws (non-centred param).\n")
    cat("  Use greta::calculate(model$greta$greta_arrays$u_variety, values = draws)\n")
  }

  results$m0_bay <- list(
    vc      = m0_bay$extras$variance_comps,
    blups   = blup0_bay,
    summary = m0_bay$extras$summary
  )
}


###############################################################################
#  MODEL 1 -- B-Spline Smoothing of Grain-Fill Curves (Per-Plot)
###############################################################################
#
#  STATISTICAL MODEL:
#    grain_weight_ij(t) = f_plot_ij(t) + epsilon_ij(t)
#
#    where f_plot(t) = sum_k alpha_k * B_k(t)  [B-spline representation]
#
#    - Fixed:  B-spline coefficients for each plot (OLS fit)
#    - This is NOT yet a mixed model -- just curve fitting.
#
#  WHAT FDA ADDS:
#    Instead of fitting a parametric logistic/Gompertz to each plot and
#    extracting rate/max parameters (the standard approach), we represent
#    each curve nonparametrically as a linear combination of B-spline
#    basis functions. This is more flexible:
#      - No assumption about curve shape (handles non-monotone, plateaus)
#      - Smoothness controlled by penalty, not by parametric form
#      - Basis coefficients become the "functional data" for downstream models
#
#  BIOLOGICAL INTERPRETATION:
#    Each plot's grain-fill trajectory is a smooth curve. We represent it
#    with ~13 B-spline basis functions (10 knots, degree 3). The penalty
#    prevents overfitting to the 8 time points.
#
###############################################################################

cat("\n", strrep("=", 70), "\nMODEL 1: B-spline smoothing of grain-fill curves\n",
    strrep("=", 70), "\n")

# ---- 1a. Build B-spline basis ----
times <- sort(unique(dt$time))
basis <- bspline_basis(times, n_knots = 4, degree = 3, penalty_order = 2)

cat("Basis functions:", basis$n_basis,
    " | Knots:", length(basis$knots),
    " | Degree:", basis$degree, "\n")

# Mixed model reparameterisation: B = [X_null, Z_range]
# Null space (polynomial) -> fixed effects
# Range space (wiggly) -> random effects with identity covariance
# This is the key link: penalised spline <-> mixed model
mm <- make_Zspline(basis, constraint = "decompose")
cat("Fixed (null space) cols:", ncol(mm$X),
    " | Random (range space) cols:", ncol(mm$Z), "\n")

# ---- 1b. Fit per-plot OLS curves (no mixed model yet) ----
# Create basis matrix for all observations
B_all <- bspline_basis(dt$time, n_knots = 4, degree = 3,
                        boundary = basis$boundary)$B

# OLS fit per plot: alpha_plot = (B'B)^{-1} B' y_plot
plots <- unique(dt$plot_id)
alpha_hat <- matrix(NA, nrow = length(plots), ncol = ncol(B_all))
rownames(alpha_hat) <- as.character(plots)

for (p in seq_along(plots)) {
  idx <- which(dt$plot_id == plots[p])
  y_p <- dt$grain_weight[idx]
  B_p <- B_all[idx, ]
  alpha_hat[p, ] <- solve(crossprod(B_p), crossprod(B_p, y_p))
}

# Evaluate on a fine grid for plotting
t_fine <- seq(min(times), max(times), length.out = 200)
B_fine <- bspline_basis(t_fine, n_knots = 4, degree = 3,
                         boundary = basis$boundary)$B
fitted_curves <- B_fine %*% t(alpha_hat)  # 200 x n_plots

cat("Fitted", ncol(fitted_curves), "plot-level curves\n")

# ---- 1c. Visualise a subset ----
# Pick 4 varieties, 1 block each
show_vars <- levels(dt$variety)[c(1, 5, 10, 15)]
show_plots <- dt[variety %in% show_vars & block == "B1",
                  unique(as.character(plot_id))]

p1 <- ggplot() +
  theme_minimal(base_size = 12) +
  labs(title = "Model 1: Per-Plot B-Spline Smoothing of Grain-Fill",
       x = "Days after anthesis", y = "Grain weight (g)",
       colour = "Plot (variety)") +
  scale_colour_viridis_d()

for (i in seq_along(show_plots)) {
  pid <- show_plots[i]
  vid <- dt[plot_id == pid, unique(as.character(variety))]
  curve_dt <- data.table(time = t_fine,
                          fitted = fitted_curves[, pid],
                          plot_id = pid, variety = vid)
  obs_dt <- dt[plot_id == pid, .(time, grain_weight, variety)]
  p1 <- p1 +
    geom_line(data = curve_dt, aes(x = time, y = fitted, colour = variety)) +
    geom_point(data = obs_dt, aes(x = time, y = grain_weight, colour = variety),
               size = 2)
}
print(p1)
cat("Plot 1 rendered.\n")


###############################################################################
#  MODEL 2 -- Variety-Specific Functional Profiles (Mixed Model)
###############################################################################
#
#  STATISTICAL MODEL:
#    grain_weight_ijk(t) = mu(t) + f_variety_i(t) + block_j + epsilon_ijk(t)
#
#    where:
#      mu(t) = sum_k beta_k * B_k(t)        [population mean curve, FIXED]
#      f_variety_i(t) = sum_k u_ik * B_k(t)  [variety deviation, RANDOM]
#      u_i ~ N(0, sigma_v^2 * I)             [iid across varieties]
#      epsilon ~ N(0, sigma_e^2)              [iid residual]
#
#  WHAT FDA ADDS BEYOND MODEL 0:
#    Model 0 analyses yield (scalar). Model 2 analyses the full grain-fill
#    TRAJECTORY. Each variety gets a smooth curve, not just a mean. We can
#    see WHERE in time varieties differ -- early fill, late fill, peak rate.
#    The P-spline mixed model formulation estimates smoothness via REML
#    (ASReml) or posterior (bayesreml).
#
#  BIOLOGICAL INTERPRETATION:
#    Varieties differ not just in final grain weight but in their filling
#    dynamics. Some fill fast early then plateau; others fill steadily.
#    The variety-specific curves f_variety_i(t) capture these dynamics.
#    Block effects remove spatial field trends.
#
###############################################################################

cat("\n", strrep("=", 70), "\nMODEL 2: Variety-specific functional profiles (mixed model)\n",
    strrep("=", 70), "\n")

# ---- 2a. Construct design matrices ----

# B-spline basis evaluated at observation time points
B_obs <- bspline_basis(dt$time, n_knots = 4, degree = 3,
                        boundary = basis$boundary)$B
n_basis <- ncol(B_obs)

# Mixed model reparameterisation
mm <- make_Zspline(basis, constraint = "decompose")
n_fixed_spline <- ncol(mm$X)
n_rand_spline  <- ncol(mm$Z)

# Evaluate X and Z at observation times
B_obs_full <- bspline_basis(dt$time, n_knots = 4, degree = 3,
                             boundary = basis$boundary)$B
# The reparameterisation transforms: B = [X_null | Z_range] %*% T
# We need X_null and Z_range evaluated at observation times
# For simplicity, use the full basis B and let ASReml handle the penalty

# Create variety:basis interaction columns
# Each variety gets its own set of B-spline coefficients
n_var <- nlevels(dt$variety)
cat("Constructing Z matrix:", n_var, "varieties x", n_basis, "basis =",
    n_var * n_basis, "random effects\n")

# Variety-specific spline: interaction of variety factor with each B-spline column
# Add basis columns to data
for (k in seq_len(n_basis)) {
  set(dt, j = paste0("B", k), value = B_obs[, k])
}

# ---- 2b. ASReml ----
if (HAS_ASREML) {
  cat("\n--- ASReml ---\n")

  # Population mean curve (fixed B-spline) + block
  # Variety-specific deviations (random B-spline coefficients)
  # Using idv() to keep coefficients iid within variety
  #
  # Fixed:  grain_weight ~ block + B1 + B2 + ... + B8
  #         (B-spline basis columns model the population mean curve)
  # Random: ~ variety:B1 + variety:B2 + ... + variety:B8
  #         (variety-specific deviations from the mean curve)

  fixed_form <- as.formula(
    paste("grain_weight ~ block +",
          paste0("B", 1:n_basis, collapse = " + "))
  )
  random_form <- as.formula(
    paste("~",
          paste0("variety:B", 1:n_basis, collapse = " + "))
  )

  m2_asr <- asreml(
    fixed   = fixed_form,
    random  = random_form,
    data    = dt,
    trace   = FALSE
  )
  cat("Converged:", m2_asr$converge, "\n")

  vc2_asr <- summary(m2_asr)$varcomp
  cat("\nVariance components (first 5):\n")
  print(head(vc2_asr, 5))
  cat("... (", nrow(vc2_asr), "total components)\n")

  # Extract variety-specific B-spline BLUPs
  # Each variety:Bk term gives a BLUP for that coefficient
  coef2_asr <- coef(m2_asr)$random
  cat("\nRandom coefficients extracted:", length(coef2_asr), "\n")

  results$m2_asr <- list(
    vc      = vc2_asr,
    loglik  = m2_asr$loglik,
    coef    = coef2_asr
  )
}

# ---- 2c. bayesreml ----
if (HAS_BAYESREML) {
  cat("\n--- bayesreml ---\n")

  # bayesreml uses the same formula syntax.
  # The B-spline columns are treated as covariates.
  m2_bay <- bayesreml(
    fixed     = as.formula(
      paste("grain_weight ~ block +",
            paste0("B", 1:n_basis, collapse = " + "))
    ),
    random    = as.formula(
      paste("~",
            paste0("variety:B", 1:n_basis, collapse = " + "))
    ),
    data      = dt,
    n_samples = 2000,
    warmup    = 1000,
    chains    = 4,
    verbose   = FALSE,
    mcmc_verbose = FALSE
  )

  cat("\nVariance components (posterior):\n")
  print(m2_bay$extras$variance_comps)

  results$m2_bay <- list(
    vc      = m2_bay$extras$variance_comps,
    summary = m2_bay$extras$summary
  )
}

# ---- 2d. Reconstruct and visualise variety curves ----
if (HAS_ASREML) {
  cat("\nReconstructing variety curves from ASReml BLUPs...\n")

  # Population mean curve (fixed effect B-spline coefficients)
  # ASReml v4.2: coef()$fixed is a matrix; rownames are coefficient names
  fixed_mat <- coef(m2_asr)$fixed
  b_rows <- grep("^B[0-9]", rownames(fixed_mat))
  beta_spline <- fixed_mat[b_rows, 1]

  # Mean curve on fine grid
  B_fine <- bspline_basis(t_fine, n_knots = 4, degree = 3,
                           boundary = basis$boundary)$B
  mean_curve <- as.numeric(B_fine %*% beta_spline)

  # Variety deviations: parse random coefficients
  # ASReml v4.2: coef()$random is a matrix with rownames "variety_V01:B1", etc.
  rand_mat <- coef(m2_asr)$random
  rand_names <- rownames(rand_mat)
  var_curves <- data.table()

  for (v in levels(dt$variety)) {
    v_idx <- grep(paste0("variety_", v, ":B"), rand_names)
    if (length(v_idx) == n_basis) {
      u_v <- rand_mat[v_idx, 1]
      variety_curve <- mean_curve + as.numeric(B_fine %*% u_v)
      var_curves <- rbind(var_curves,
                          data.table(time = t_fine, fitted = variety_curve,
                                     variety = v))
    }
  }

  if (nrow(var_curves) > 0) {
    # Plot selected varieties
    show_dt <- var_curves[variety %in% show_vars]
    mean_dt <- data.table(time = t_fine, fitted = mean_curve)

    p2 <- ggplot(show_dt, aes(x = time, y = fitted, colour = variety)) +
      geom_line(linewidth = 0.8) +
      geom_line(data = mean_dt, aes(x = time, y = fitted),
                colour = "black", linewidth = 1.2, linetype = "dashed") +
      # Overlay raw data
      geom_point(data = dt[variety %in% show_vars],
                 aes(x = time, y = grain_weight, colour = variety),
                 alpha = 0.4, size = 1) +
      scale_colour_viridis_d() +
      theme_minimal(base_size = 12) +
      labs(title = "Model 2: Variety-Specific Grain-Fill Curves (ASReml)",
           subtitle = "Dashed = population mean; solid = variety curve",
           x = "Days after anthesis", y = "Grain weight (g)")
    print(p2)
    cat("Plot 2 rendered.\n")
  }
}


###############################################################################
#  MODEL 3 -- Scalar-on-Function: Yield ~ Grain-Fill Trajectory
###############################################################################
#
#  STATISTICAL MODEL (Two-Stage):
#    Stage 1: Estimate variety curves f_hat_v(t) from Model 2
#    Stage 2: yield_v = alpha + integral[ beta(t) * f_hat_v(t) ] dt + epsilon_v
#
#    where:
#      beta(t) = sum_k b_k * B_k(t)   [coefficient function, penalised]
#      The integral is approximated as:
#        integral[ beta(t) * f_hat(t) ] dt ~= alpha_v' * J * b
#        where J_{jk} = integral[ B_j(t) * B_k(t) ] dt
#              alpha_v = B-spline coefficients of variety v's curve
#              b = B-spline coefficients of beta(t)
#
#    - Fixed:  intercept + functional covariate (C = Alpha %*% J, penalised)
#    - Random: none (or genomic/pedigree if available)
#    - Residual: epsilon ~ N(0, sigma_e^2)
#
#  WHAT FDA ADDS:
#    Standard approach: fit logistic to grain-fill, extract rate/Wmax,
#    run yield ~ rate + Wmax. This assumes yield depends on those specific
#    parameters and the relationship is linear.
#
#    FDA approach: yield depends on the ENTIRE curve shape. beta(t) tells us
#    WHEN during grain-fill the trajectory most influences yield. If beta(t)
#    is large at t=25 (mid-fill) and near zero at t=10 (early), then mid-fill
#    dynamics drive yield. No parametric assumption about the curve or the
#    yield-curve relationship.
#
#  BIOLOGICAL INTERPRETATION:
#    beta(t) is the "importance function" -- it reveals which phase of grain
#    filling matters most for yield. Breeders can target selection pressure
#    on the critical growth window identified by beta(t).
#
###############################################################################

cat("\n", strrep("=", 70), "\nMODEL 3: Scalar-on-function -- yield ~ grain-fill curve\n",
    strrep("=", 70), "\n")

# ---- 3a. Compute functional covariate matrix ----

# Use per-variety mean curves (averaged over blocks) from OLS (Model 1)
# alpha_hat is plots x basis. Average over blocks to get variety means.
plot_var_map <- unique(dt[, .(plot_id, variety)])
alpha_variety <- matrix(NA, nrow = n_var, ncol = n_basis)
rownames(alpha_variety) <- levels(dt$variety)

for (v in levels(dt$variety)) {
  v_plots <- plot_var_map[variety == v, as.character(plot_id)]
  alpha_variety[v, ] <- colMeans(alpha_hat[v_plots, , drop = FALSE])
}

# Inner product matrix J: J_{jk} = integral B_j(t) * B_k(t) dt
# Numerical integration via Simpson's rule on a fine grid
t_quad <- seq(min(times), max(times), length.out = 500)
B_quad <- bspline_basis(t_quad, n_knots = 4, degree = 3,
                         boundary = basis$boundary)$B
dt_quad <- diff(t_quad[1:2])  # uniform spacing
J <- crossprod(B_quad) * dt_quad  # (n_basis x n_basis) Gram matrix

cat("Inner product matrix J:", nrow(J), "x", ncol(J), "\n")
cat("Condition number:", kappa(J), "\n")

# Functional covariate: C = Alpha %*% J (n_var x n_basis)
# C_vk = sum_j alpha_vj * J_jk
# Then: yield_v ~ intercept + C_v %*% b
C_mat <- alpha_variety %*% as.matrix(J)
cat("Functional covariate C:", nrow(C_mat), "x", ncol(C_mat), "\n")

# Build regression data.frame
reg_dt <- merge(
  data.table(variety = levels(dt$variety)),
  yield_dt, by = "variety", sort = FALSE
)
# Add C columns
for (k in seq_len(n_basis)) {
  set(reg_dt, j = paste0("C", k), value = C_mat[, k])
}

# ---- 3b. ASReml (scalar-on-function) ----
if (HAS_ASREML) {
  cat("\n--- ASReml ---\n")

  # Fixed effects: intercept + C1..C8 (functional covariate columns)
  # Ideally we'd penalise beta(t) via a random effect formulation,
  # but for simplicity we fit as fixed first, then add penalty.
  fixed_sof <- as.formula(
    paste("yield ~ 1 +", paste0("C", 1:n_basis, collapse = " + "))
  )

  m3_asr <- asreml(
    fixed = fixed_sof,
    data  = reg_dt,
    trace = FALSE
  )
  cat("Converged:", m3_asr$converge, "\n")

  # Extract beta(t) coefficient function
  # ASReml v4.2: coef()$fixed is a matrix with rownames
  m3_fixed <- coef(m3_asr)$fixed
  b_hat_asr <- m3_fixed[paste0("C", 1:n_basis), 1]
  beta_t_asr <- as.numeric(B_fine %*% b_hat_asr)

  cat("\nFixed effects (B-spline coefficients for beta(t)):\n")
  print(round(b_hat_asr, 4))

  # Model fit
  cat("Log-likelihood:", m3_asr$loglik, "\n")

  # Predictions
  pred3_asr <- predict(m3_asr, classify = "variety")
  if (!is.null(pred3_asr$pvals)) {
    pred3_dt <- as.data.table(pred3_asr$pvals)
    cat("\nPredicted yields (top 5):\n")
    print(head(pred3_dt[order(-predicted.value)], 5))
  }

  results$m3_asr <- list(
    beta_coef = b_hat_asr,
    beta_t    = data.table(time = t_fine, beta = beta_t_asr),
    loglik    = m3_asr$loglik
  )
}

# ---- 3c. bayesreml (scalar-on-function) ----
if (HAS_BAYESREML) {
  cat("\n--- bayesreml ---\n")

  m3_bay <- bayesreml(
    fixed     = as.formula(
      paste("yield ~ 1 +", paste0("C", 1:n_basis, collapse = " + "))
    ),
    data      = reg_dt,
    n_samples = 3000,
    warmup    = 1500,
    chains    = 4,
    verbose   = FALSE,
    mcmc_verbose = FALSE
  )

  # Extract beta(t)
  bay_summary <- m3_bay$extras$summary
  b_hat_bay <- bay_summary[grep("^beta_C", bay_summary$param), "mean"]
  if (length(b_hat_bay) == n_basis) {
    beta_t_bay <- as.numeric(B_fine %*% b_hat_bay)
  } else {
    # Try extracting from fixed coefficients
    b_hat_bay <- coef(m3_bay)[paste0("C", 1:n_basis)]
    beta_t_bay <- as.numeric(B_fine %*% b_hat_bay)
  }

  cat("\nBeta(t) coefficients (posterior means):\n")
  print(round(b_hat_bay, 4))

  results$m3_bay <- list(
    beta_coef = b_hat_bay,
    beta_t    = data.table(time = t_fine, beta = beta_t_bay),
    summary   = bay_summary
  )
}

# ---- 3d. Visualise beta(t) ----
if (exists("beta_t_asr") || exists("beta_t_bay")) {
  p3_dt <- data.table()
  if (exists("beta_t_asr")) {
    p3_dt <- rbind(p3_dt,
                   data.table(time = t_fine, beta = beta_t_asr,
                              method = "ASReml (REML)"))
  }
  if (exists("beta_t_bay")) {
    p3_dt <- rbind(p3_dt,
                   data.table(time = t_fine, beta = beta_t_bay,
                              method = "bayesreml (Bayesian)"))
  }

  p3 <- ggplot(p3_dt, aes(x = time, y = beta, colour = method)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("ASReml (REML)" = "#440154",
                                   "bayesreml (Bayesian)" = "#21918c")) +
    theme_minimal(base_size = 12) +
    labs(title = "Model 3: Coefficient Function beta(t)",
         subtitle = "How each phase of grain-fill influences yield",
         x = "Days after anthesis",
         y = expression(beta(t)),
         colour = "Method")
  print(p3)
  cat("Plot 3 rendered.\n")
}


###############################################################################
#  MODEL 4 -- Penalised Scalar-on-Function (Mixed Model Formulation)
###############################################################################
#
#  STATISTICAL MODEL:
#    Same as Model 3 but with PENALISED beta(t):
#
#    yield_v = alpha + C_v' * b + epsilon_v
#    b = [b_fixed, b_random]'
#
#    The P-spline penalty on beta(t) is enforced by decomposing b into:
#      - b_fixed: null space (polynomial trend in beta(t)) -- FIXED
#      - b_random: range space (wiggly deviations) -- RANDOM, var = sigma_b^2
#
#    So: yield ~ 1 + C_null * b_fixed + C_range * b_random
#    where C_null = C %*% T_X,  C_range = C %*% T_Z
#    and T_X, T_Z come from make_Zspline() decomposition.
#
#  WHAT FDA ADDS:
#    The penalty prevents overfitting beta(t) to the 20 variety yield
#    observations (n=20 with p=8 basis columns is ill-conditioned).
#    Smoothness of beta(t) is estimated from data via REML or posterior.
#
#  BIOLOGICAL INTERPRETATION:
#    Same as Model 3 but the smoothed beta(t) gives a more reliable picture
#    of which growth phases matter. Prevents spurious wiggles from noise.
#
###############################################################################

cat("\n", strrep("=", 70), "\nMODEL 4: Penalised scalar-on-function (P-spline on beta(t))\n",
    strrep("=", 70), "\n")

# ---- 4a. Decompose functional covariate ----
# C_null: functional covariate projected onto null space of penalty
# C_range: projected onto range space
C_null  <- C_mat %*% mm$X  # n_var x penalty_order
C_range <- C_mat %*% mm$Z  # n_var x (n_basis - penalty_order)

cat("C_null:", nrow(C_null), "x", ncol(C_null), "\n")
cat("C_range:", nrow(C_range), "x", ncol(C_range), "\n")

# Add to regression data
for (k in seq_len(ncol(C_null))) {
  set(reg_dt, j = paste0("Cnull", k), value = C_null[, k])
}
for (k in seq_len(ncol(C_range))) {
  set(reg_dt, j = paste0("Crange", k), value = C_range[, k])
}

# ---- 4b. ASReml (penalised) ----
if (HAS_ASREML) {
  cat("\n--- ASReml ---\n")

  # Fixed: intercept + null-space columns (unpenalised polynomial)
  # Random: range-space columns (penalised, variance = sigma_b^2)
  fixed_pen  <- as.formula(
    paste("yield ~ 1 +",
          paste0("Cnull", 1:ncol(C_null), collapse = " + "))
  )
  random_pen <- as.formula(
    paste("~",
          paste0("Crange", 1:ncol(C_range), collapse = " + "))
  )

  m4_asr <- asreml(
    fixed   = fixed_pen,
    random  = random_pen,
    data    = reg_dt,
    trace   = FALSE
  )
  cat("Converged:", m4_asr$converge, "\n")

  vc4_asr <- summary(m4_asr)$varcomp
  cat("\nVariance components:\n")
  print(vc4_asr)

  # Reconstruct beta(t) from fixed + random parts
  m4_fixed <- coef(m4_asr)$fixed
  m4_random <- coef(m4_asr)$random
  b_null_asr  <- m4_fixed[paste0("Cnull", 1:ncol(C_null)), 1]
  b_range_asr <- m4_random[paste0("Crange", 1:ncol(C_range)), 1]

  if (length(b_null_asr) > 0 && length(b_range_asr) > 0) {
    # Transform back to original basis: b = T_X %*% b_null + T_Z %*% b_range
    # But we need [X_fine, Z_fine] evaluated at fine grid
    B_fine_null  <- B_fine %*% mm$X
    B_fine_range <- B_fine %*% mm$Z
    beta_t_pen_asr <- as.numeric(
      B_fine_null %*% b_null_asr + B_fine_range %*% b_range_asr
    )
    cat("Penalised beta(t) reconstructed.\n")
    results$m4_asr <- list(
      beta_t = data.table(time = t_fine, beta = beta_t_pen_asr),
      vc     = vc4_asr
    )
  }
}

# ---- 4c. bayesreml (penalised) ----
if (HAS_BAYESREML) {
  cat("\n--- bayesreml ---\n")

  m4_bay <- bayesreml(
    fixed     = as.formula(
      paste("yield ~ 1 +",
            paste0("Cnull", 1:ncol(C_null), collapse = " + "))
    ),
    random    = as.formula(
      paste("~",
            paste0("Crange", 1:ncol(C_range), collapse = " + "))
    ),
    data      = reg_dt,
    n_samples = 3000,
    warmup    = 1500,
    chains    = 4,
    verbose   = FALSE,
    mcmc_verbose = FALSE
  )

  cat("\nVariance components:\n")
  print(m4_bay$extras$variance_comps)

  # Reconstruct beta(t)
  b_null_bay  <- coef(m4_bay)[paste0("Cnull", 1:ncol(C_null))]
  b_range_bay <- m4_bay$extras$blups
  # Extract range coefficients
  if (is.list(b_range_bay)) {
    b_r <- unlist(b_range_bay[grep("Crange", names(b_range_bay))])
    if (length(b_r) == ncol(C_range)) {
      B_fine_null  <- B_fine %*% mm$X
      B_fine_range <- B_fine %*% mm$Z
      beta_t_pen_bay <- as.numeric(
        B_fine_null %*% b_null_bay + B_fine_range %*% b_r
      )
      results$m4_bay <- list(
        beta_t = data.table(time = t_fine, beta = beta_t_pen_bay),
        vc     = m4_bay$extras$variance_comps
      )
    }
  }
}

# ---- 4d. Compare penalised vs unpenalised beta(t) ----
p4_dt <- data.table()
if (!is.null(results$m3_asr)) {
  p4_dt <- rbind(p4_dt, cbind(results$m3_asr$beta_t, method = "Unpenalised (ASReml)"))
}
if (!is.null(results$m4_asr)) {
  p4_dt <- rbind(p4_dt, cbind(results$m4_asr$beta_t, method = "Penalised (ASReml)"))
}
if (!is.null(results$m4_bay)) {
  p4_dt <- rbind(p4_dt, cbind(results$m4_bay$beta_t, method = "Penalised (bayesreml)"))
}

if (nrow(p4_dt) > 0) {
  p4 <- ggplot(p4_dt, aes(x = time, y = beta, colour = method, linetype = method)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_colour_viridis_d(option = "C") +
    theme_minimal(base_size = 12) +
    labs(title = "Model 4: Penalised vs Unpenalised beta(t)",
         subtitle = "Penalty smooths out spurious wiggles",
         x = "Days after anthesis", y = expression(beta(t)))
  print(p4)
  cat("Plot 4 rendered.\n")
}


###############################################################################
#  MODEL 5 -- funcrop High-Level API: Two-Stage FDA
###############################################################################
#
#  STATISTICAL MODEL:
#    Stage 1: fit_functional_profiles() -- Model 2 via funcrop API
#    Stage 2: scalar_on_function()      -- Model 4 via funcrop API
#
#  This demonstrates the streamlined funcrop workflow where the user
#  specifies the model in domain terms, not matrix algebra.
#
#  WHAT funcrop ADDS:
#    All the basis construction, Z-matrix Kronecker products, penalty
#    decomposition, BLUP extraction, beta(t) reconstruction, and
#    visualisation are automated. The user just says:
#      "Fit variety-specific grain-fill curves, then regress yield on them."
#
###############################################################################

cat("\n", strrep("=", 70), "\nMODEL 5: funcrop high-level API (two-stage FDA)\n",
    strrep("=", 70), "\n")

# Clean the extra B-spline columns before passing to funcrop
dt_clean <- dt[, .(plot_id, variety, block, row, col, time,
                    grain_weight, yield)]

cat("Using funcrop::fit_functional_profiles() + scalar_on_function()\n")
cat("Engine: ASReml if available, else bayesreml\n\n")

# The funcrop API handles everything:
for (eng in c("asreml", "bayesreml")) {
  if (eng == "asreml" && !HAS_ASREML) next
  if (eng == "bayesreml" && !HAS_BAYESREML) next

  cat(sprintf("\n--- funcrop with engine = '%s' ---\n", eng))

  tryCatch({
    # Stage 1: Variety-specific functional profiles
    fit_s1 <- fit_functional_profiles(
      data       = dt_clean,
      time_col   = "time",
      value_col  = "grain_weight",
      id_col     = "plot_id",
      group_col  = "variety",
      n_knots    = 4,
      degree     = 3,
      spatial    = "none",
      engine     = eng
    )

    cat("Stage 1 complete. Engine:", fit_s1$engine, "\n")

    # Stage 2: Relate yield to functional profiles
    fit_s2 <- scalar_on_function(
      primary_trait       = yield_dt[["yield"]],
      functional_profiles = fit_s1,
      engine              = eng
    )

    cat("Stage 2 complete.\n")

    # Store
    results[[paste0("m5_", eng)]] <- list(
      stage1 = fit_s1,
      stage2 = fit_s2
    )

    # Plot using funcrop's built-in visualisation
    tryCatch({
      plot_functional_profiles(fit_s1,
                                varieties = show_vars,
                                ci = TRUE)
    }, error = function(e) cat("Plot error:", e$message, "\n"))

    tryCatch({
      plot_coefficient_function(fit_s2, ci = TRUE)
    }, error = function(e) cat("Plot error:", e$message, "\n"))

  }, error = function(e) {
    cat("Error with engine '", eng, "':", e$message, "\n")
  })
}


###############################################################################
#  MODEL 6 -- MET-FDA with GxE (Multi-Environment Extension)
###############################################################################
#
#  STATISTICAL MODEL:
#    For each environment e:
#      y_ijk(t) = mu_e(t) + f_{v,e}(t) + block_j + spatial_e(row,col) + eps
#
#    Across environments:
#      vec(alpha_v) ~ N(0, Sigma_GxE (x) I_spline)
#      Sigma_GxE modelled as FA1: Lambda * Lambda' + Psi (diagonal)
#
#    - Fixed:  environment-specific intercepts, environment:basis interactions
#    - Random: variety:environment B-spline coefficients with FA structure
#    - GxE:   captured by FA loadings on variety spline coefficients
#
#  WHAT FDA ADDS BEYOND STANDARD MET:
#    Standard MET: yield_ve = mu_e + g_v + (ge)_ve with FA on (ge).
#    FDA-MET: the GxE interaction operates on the functional SHAPE, not just
#    the scalar yield. We can see that variety V01 fills grain faster in
#    environment E1 but slower in E3 -- the interaction is on the trajectory.
#
#  BIOLOGICAL INTERPRETATION:
#    Some varieties' grain-fill response is stable across environments
#    (low GxE on the curve). Others are sensitive -- their filling rate or
#    duration shifts across environments. The FA loadings reveal which
#    environments are similar in how they affect grain-fill dynamics.
#
###############################################################################

cat("\n", strrep("=", 70), "\nMODEL 6: MET-FDA with GxE (using sim_met_fda)\n",
    strrep("=", 70), "\n")

data(sim_met_fda)
met <- copy(sim_met_fda)
met[, variety     := factor(variety)]
met[, environment := factor(environment)]
met[, rep         := factor(rep)]
met[, iblock      := factor(iblock)]

cat("MET data:", nrow(met), "rows |", uniqueN(met$variety), "varieties |",
    uniqueN(met$environment), "environments\n\n")

# For this demonstration, we show the manual ASReml approach for a single
# environment, then the funcrop API for the full MET.

# ---- 6a. Single-environment functional model (ASReml) ----
if (HAS_ASREML) {
  cat("--- Single environment (E1) via ASReml ---\n")
  met_e1 <- met[environment == "E1"]

  # B-spline basis for NDVI time series
  met_times <- sort(unique(met$time))
  met_basis <- bspline_basis(met_times, n_knots = 3, degree = 3)
  B_met <- bspline_basis(met_e1$time, n_knots = 3, degree = 3,
                          boundary = met_basis$boundary)$B

  for (k in seq_len(ncol(B_met))) {
    set(met_e1, j = paste0("MB", k), value = B_met[, k])
  }

  n_mbasis <- ncol(B_met)
  m6_e1 <- asreml(
    fixed  = as.formula(
      paste("ndvi ~ rep +", paste0("MB", 1:n_mbasis, collapse = " + "))
    ),
    random = as.formula(
      paste("~", paste0("variety:MB", 1:n_mbasis, collapse = " + "))
    ),
    data   = met_e1,
    trace  = FALSE
  )
  cat("E1 model converged:", m6_e1$converge, "\n")
  cat("Variance components:\n")
  print(head(summary(m6_e1)$varcomp, 5))

  results$m6_e1_asr <- list(
    vc     = summary(m6_e1)$varcomp,
    loglik = m6_e1$loglik
  )
}

# ---- 6b. Full MET-FDA via funcrop API ----
cat("\n--- Full MET-FDA via funcrop ---\n")
for (eng in c("asreml", "bayesreml")) {
  if (eng == "asreml" && !HAS_ASREML) next
  if (eng == "bayesreml" && !HAS_BAYESREML) next

  cat(sprintf("\nfuncrop MET-FDA with engine = '%s'\n", eng))
  tryCatch({
    fit_met <- fit_fda_met(
      data            = met,
      environment_col = "environment",
      time_col        = "time",
      value_col       = "ndvi",
      id_col          = "plot_id",
      group_col       = "variety",
      primary_col     = "yield",
      n_knots         = 3,
      degree          = 3,
      gxe_structure   = "diag",  # Start simple; FA1 for larger datasets
      spatial         = "none",
      two_stage       = TRUE,
      engine          = eng
    )

    cat("MET-FDA complete. Engine:", fit_met$engine, "\n")
    results[[paste0("m6_met_", eng)]] <- fit_met

  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
}


###############################################################################
#  COMPARISON AND RECONCILIATION
###############################################################################

cat("\n\n", strrep("=", 70), "\n")
cat("COMPARISON: ASReml (REML) vs bayesreml (Bayesian)\n")
cat(strrep("=", 70), "\n\n")

# ---- Model 0: Baseline Yield ----
if (!is.null(results$m0_asr) && !is.null(results$m0_bay)) {
  cat("--- Model 0: Baseline RCBD for Yield ---\n")

  # Variance components
  asr_sigv <- results$m0_asr$vc["variety", "component"]
  asr_sige <- results$m0_asr$vc["units!R", "component"]
  bay_vc   <- results$m0_bay$vc

  cat(sprintf("  sigma_variety^2:  ASReml = %.4f  |  bayesreml = %s\n",
              asr_sigv,
              ifelse(is.data.frame(bay_vc),
                     sprintf("%.4f", bay_vc$estimate[grep("variety", bay_vc$component)]),
                     "N/A")))
  cat(sprintf("  sigma_resid^2:    ASReml = %.4f  |  bayesreml = %s\n",
              asr_sige,
              ifelse(is.data.frame(bay_vc),
                     sprintf("%.4f", bay_vc$estimate[grep("units|resid", bay_vc$component)]),
                     "N/A")))

  # BLUP comparison
  asr_blups <- results$m0_asr$blups
  bay_blups <- results$m0_bay$blups
  if (!is.null(asr_blups) && !is.null(bay_blups)) {
    merged <- merge(asr_blups, bay_blups,
                    by = "variety", suffixes = c("_asr", "_bay"))
    if ("predicted.value_asr" %in% names(merged) &&
        "predicted.value_bay" %in% names(merged)) {
      blup_cor <- cor(merged$predicted.value_asr, merged$predicted.value_bay)
      cat(sprintf("  BLUP correlation: r = %.4f\n", blup_cor))
    }
  }
}

# ---- Model 3: Beta(t) comparison ----
if (!is.null(results$m3_asr) && !is.null(results$m3_bay)) {
  cat("\n--- Model 3: Unpenalised Scalar-on-Function ---\n")
  b_asr <- results$m3_asr$beta_coef
  b_bay <- results$m3_bay$beta_coef
  if (length(b_asr) == length(b_bay)) {
    cat("  B-spline coefficient comparison:\n")
    comp_dt <- data.table(
      basis   = paste0("b", 1:length(b_asr)),
      ASReml  = round(as.numeric(b_asr), 5),
      bayesreml = round(as.numeric(b_bay), 5),
      diff    = round(as.numeric(b_asr) - as.numeric(b_bay), 5)
    )
    print(comp_dt)
    cat(sprintf("  Max absolute difference: %.5f\n", max(abs(comp_dt$diff))))
    cat(sprintf("  Correlation: %.4f\n", cor(comp_dt$ASReml, comp_dt$bayesreml)))
  }
}

# ---- Summary table ----
cat("\n\n")
cat(strrep("=", 70), "\n")
cat("SUMMARY TABLE: Key Results Across All Models\n")
cat(strrep("=", 70), "\n\n")

summary_rows <- list()

# Populate from available results
if (!is.null(results$m0_asr)) {
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    Model = "M0: Baseline RCBD", Engine = "ASReml",
    LogLik = round(results$m0_asr$loglik, 2),
    Note = "Scalar yield, no FDA"
  )
}
if (!is.null(results$m0_bay)) {
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    Model = "M0: Baseline RCBD", Engine = "bayesreml",
    LogLik = NA_real_,
    Note = "Bayesian -- no log-lik, use WAIC"
  )
}
if (!is.null(results$m2_asr)) {
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    Model = "M2: Variety curves", Engine = "ASReml",
    LogLik = round(results$m2_asr$loglik, 2),
    Note = paste(nrow(results$m2_asr$vc), "variance components")
  )
}
if (!is.null(results$m3_asr)) {
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    Model = "M3: Yield ~ f(t)", Engine = "ASReml",
    LogLik = round(results$m3_asr$loglik, 2),
    Note = "Unpenalised beta(t)"
  )
}
if (!is.null(results$m4_asr)) {
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    Model = "M4: Yield ~ f(t) penalised", Engine = "ASReml",
    LogLik = NA_real_,
    Note = paste("sigma_beta^2 =",
                 round(results$m4_asr$vc[1, "component"], 4))
  )
}
if (!is.null(results$m6_e1_asr)) {
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    Model = "M6: MET single env (E1)", Engine = "ASReml",
    LogLik = round(results$m6_e1_asr$loglik, 2),
    Note = "NDVI functional model"
  )
}

if (length(summary_rows) > 0) {
  summary_table <- rbindlist(summary_rows, fill = TRUE)
  print(summary_table)
} else {
  cat("No results to summarise (check that ASReml/bayesreml are available).\n")
}


###############################################################################
#  EXPECTED INCONSISTENCIES AND RECONCILIATION
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("RECONCILIATION: Expected Differences Between REML and Bayesian\n")
cat(strrep("=", 70), "\n\n")

cat("
1. VARIANCE COMPONENTS
   - ASReml (REML): point estimates via restricted maximum likelihood.
     Can hit boundary (zero) for small variance components.
   - bayesreml (Bayesian): posterior means with lognormal priors.
     Priors shrink toward prior mean; never exactly zero.
   - RECONCILIATION: For well-identified components, REML and posterior
     means should be close. Divergence signals either weak data or
     strong prior influence. Increase n_samples and check Rhat < 1.05.

2. FIXED EFFECT ESTIMATES (including beta(t) coefficients)
   - ASReml: generalised least squares conditional on REML variance estimates.
   - bayesreml: posterior means marginalised over variance uncertainty.
   - RECONCILIATION: Typically very close (< 5%% difference) for models
     with moderate data. Bayesian estimates have slightly wider CIs
     because they account for variance component uncertainty.

3. BLUPs vs POSTERIOR MEANS
   - ASReml BLUPs: conditional modes E[u|y, theta_hat].
   - bayesreml: posterior means E[u|y] (integrated over theta).
   - RECONCILIATION: Nearly identical for large datasets. For small
     datasets, Bayesian estimates are slightly more shrunk toward zero
     due to prior regularisation.

4. MODEL FIT STATISTICS
   - ASReml: log-likelihood, AIC, BIC (frequentist).
   - bayesreml: WAIC, LOO-IC (Bayesian). Not directly comparable to AIC.
   - RECONCILIATION: Both penalise complexity but with different
     philosophies. Use each within its own framework for model selection.

5. CONVERGENCE
   - ASReml: REML iterations; may not converge for complex models.
   - bayesreml: MCMC; convergence assessed via Rhat (< 1.05) and
     ESS (> 400). Slow for high-dimensional random effects.
   - RECONCILIATION: If ASReml converges but bayesreml has high Rhat,
     increase warmup/chains. If ASReml doesn't converge, bayesreml
     with informative priors may succeed.

6. PENALTY / SMOOTHING PARAMETER
   - ASReml: sigma_b^2 estimated by REML (ratio sigma_b^2 / sigma_e^2
     controls smoothness). Can be zero (unpenalised) or very large
     (fully penalised).
   - bayesreml: sigma_b has a lognormal prior. Prior choice matters
     when n_basis >> n_observations.
   - RECONCILIATION: For Model 4 (20 obs, 6 random B-spline effects),
     the prior matters. Sensitivity analysis recommended.

7. COMPUTATIONAL COST
   - ASReml: fast (seconds for simple models, minutes for MET-FA).
   - bayesreml: slow (minutes for simple, hours for MET-FA with many
     varieties and time points). Use ASReml for exploration,
     bayesreml for final inference with full uncertainty.
")

cat("\n===== Script complete =====\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Session info:\n")
print(sessionInfo())
