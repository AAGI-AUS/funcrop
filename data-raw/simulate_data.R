# simulate_data.R — Generate simulated datasets for funcrop testing & vignettes
#
# Fully reproducible simulation of two trial datasets:
#   1. sim_grain_fill  — single-trial grain-fill experiment (RCBD)
#   2. sim_met_fda     — multi-environment trial with functional NDVI trait
#
# Author: Max Moldovan
# Date:   2026-04-04
# Licence: GPL (>= 3)

library(data.table)

# ==============================================================================
# Helper: Generate AR1 correlation matrix
# ==============================================================================

.ar1_cor <- function(n, phi) {

  d <- abs(outer(seq_len(n), seq_len(n), "-"))
  phi^d
}

# ==============================================================================
# Dataset 1: sim_grain_fill
# ==============================================================================
# Simulates a grain-fill trial:
#   - 20 varieties, RCBD with 3 blocks, 60 plots
#   - 12 rows x 5 columns spatial layout
#   - Grain weight measured at 8 time points (logistic growth)
#   - AR1 spatial correlation in rows (phi = 0.3)
#   - Primary trait: yield (related to grain-fill parameters)

set.seed(20250415)

n_var <- 20L
n_block <- 3L
n_plot <- n_var * n_block
time_pts <- c(10, 15, 20, 25, 30, 35, 40, 45)
n_time <- length(time_pts)

n_row <- 12L
n_col <- 5L

# --- Variety-level true parameters ---
wmax <- rnorm(n_var, mean = 40, sd = 5)
rate <- rnorm(n_var, mean = 0.15, sd = 0.02)
tmid <- rnorm(n_var, mean = 25, sd = 2)

# --- Block effects ---
block_eff <- rnorm(n_block, mean = 0, sd = 2)

# --- Allocate plots: RCBD (randomise variety order within each block) ---
plot_dt <- rbindlist(lapply(seq_len(n_block), function(b) {
  data.table(
    variety = paste0("V", sprintf("%02d", sample(seq_len(n_var)))),
    block   = paste0("B", b)
  )
}))
plot_dt[, plot_id := paste0("P", sprintf("%03d", .I))]

# --- Assign spatial coordinates (row, col) ---
# 60 plots laid out in a 12 x 5 grid, filled column-major
plot_dt[, row := rep(seq_len(n_row), length.out = n_plot)]
plot_dt[, col := rep(seq_len(n_col), each = n_row, length.out = n_plot)]

# --- Generate AR1 correlated spatial noise ---
# AR1 in rows only (phi = 0.3), independent across columns
phi_row <- 0.3
sigma_e <- 1.5
R_row <- .ar1_cor(n_row, phi_row)
L_row <- chol(R_row)  # upper Cholesky factor

# Generate spatially correlated residuals for each (column, time) combination
# n_plot x n_time residual matrix
spatial_noise <- matrix(0, nrow = n_plot, ncol = n_time)
for (tt in seq_len(n_time)) {
  # For each column, draw AR1-correlated residuals across its rows
  col_noise <- matrix(0, nrow = n_row, ncol = n_col)
  for (cc in seq_len(n_col)) {
    z <- rnorm(n_row)
    col_noise[, cc] <- as.numeric(crossprod(L_row, z))  # t(L) %*% z
  }
  # Map to plots (plots are ordered by row within column)
  spatial_noise[, tt] <- sigma_e * as.vector(col_noise)
}

# --- Expand to long format (plot x time) ---
# Create variety index lookup
var_idx <- as.integer(gsub("V", "", plot_dt$variety))
block_idx <- as.integer(gsub("B", "", plot_dt$block))

sim_grain_fill <- rbindlist(lapply(seq_len(n_plot), function(i) {
  vi <- var_idx[i]
  bi <- block_idx[i]

  # True logistic grain-fill curve for this variety
  true_curve <- wmax[vi] / (1 + exp(-rate[vi] * (time_pts - tmid[vi])))

  # Observed = true + block effect + spatially correlated noise
  observed <- true_curve + block_eff[bi] + spatial_noise[i, ]

  data.table(
    plot_id      = plot_dt$plot_id[i],
    variety      = plot_dt$variety[i],
    block        = plot_dt$block[i],
    row          = plot_dt$row[i],
    col          = plot_dt$col[i],
    time         = time_pts,
    grain_weight = round(observed, 2)
  )
}))

# --- Primary trait: yield ---
yield_dt <- data.table(
  plot_id = plot_dt$plot_id,
  yield   = round(
    5 + 0.1 * wmax[var_idx] + 0.5 * rate[var_idx] * 100 + rnorm(n_plot, 0, 1),
    2
  )
)

sim_grain_fill <- merge(sim_grain_fill, yield_dt, by = "plot_id")

# Ensure data.table class and clean column order
setcolorder(sim_grain_fill,
            c("plot_id", "variety", "block", "row", "col",
              "time", "grain_weight", "yield"))
setkey(sim_grain_fill, plot_id, time)


# ==============================================================================
# Dataset 2: sim_met_fda
# ==============================================================================
# Simulates a multi-environment trial with functional NDVI trait:
#   - 30 varieties across 4 environments
#   - Alpha-lattice: 3 reps, 5 incomplete blocks per rep, 6 plots per block
#   - NDVI measured at 6 time points (Gaussian decay stay-green model)
#   - GxE via FA1 structure (1 factor)
#   - AR1 x AR1 spatial correlation per environment
#   - Primary trait: yield per environment

set.seed(20250416)

n_var_met   <- 30L
n_env       <- 4L
env_names   <- paste0("E", seq_len(n_env))
n_rep       <- 3L
n_iblock    <- 5L   # incomplete blocks per rep
n_per_block <- 6L   # plots per incomplete block
n_plot_env  <- n_var_met * n_rep  # 90 plots per environment
time_ndvi   <- c(70, 80, 90, 100, 110, 120)
n_time_ndvi <- length(time_ndvi)

# --- Variety true parameters (main effects) ---
ndvi_max_main  <- runif(n_var_met, 0.7, 0.95)
decay_main     <- runif(n_var_met, 0.3, 0.8)
t_onset_main   <- runif(n_var_met, 80, 100)

# --- GxE via FA1 structure ---
env_loadings <- c(0.8, 0.5, -0.3, -0.7)
var_scores   <- rnorm(n_var_met, 0, 1)

# GxE effects matrix: n_var x n_env
gxe_matrix <- outer(var_scores, env_loadings)

# --- Environment-specific residual SDs ---
sigma_env <- c(0.03, 0.04, 0.035, 0.05)

# --- AR1 x AR1 spatial parameters ---
phi_r <- 0.4
phi_c <- 0.3

# --- Layout: determine row/col for n_plot_env = 90 plots ---
# Use 15 rows x 6 columns = 90 plots per environment
n_row_met <- 15L
n_col_met <- 6L

# --- Generate data per environment ---
sim_met_fda <- rbindlist(lapply(seq_len(n_env), function(e) {

  # Alpha-lattice design: randomise variety allocation

  # 3 reps x 5 iblocks x 6 plots = 90 plots
  design_dt <- rbindlist(lapply(seq_len(n_rep), function(r) {
    # Randomly assign all 30 varieties to 5 incomplete blocks of 6
    perm <- sample(seq_len(n_var_met))
    iblk <- rep(seq_len(n_iblock), each = n_per_block)
    data.table(
      variety = paste0("G", sprintf("%02d", perm)),
      rep     = paste0("R", r),
      iblock  = paste0("IB", iblk)
    )
  }))

  design_dt[, plot_id := paste0(env_names[e], "_P", sprintf("%03d", .I))]
  design_dt[, row := rep(seq_len(n_row_met), length.out = n_plot_env)]
  design_dt[, col := rep(seq_len(n_col_met), each = n_row_met,
                         length.out = n_plot_env)]

  # --- AR1 x AR1 spatial noise ---
  R_r <- .ar1_cor(n_row_met, phi_r)
  R_c <- .ar1_cor(n_col_met, phi_c)
  # Kronecker: R_row (x) R_col for the full n_row x n_col grid
  # Since plots map to (row, col), generate noise on the grid
  L_r <- chol(R_r)
  L_c <- chol(R_c)

  spatial_noise_env <- matrix(0, nrow = n_plot_env, ncol = n_time_ndvi)
  for (tt in seq_len(n_time_ndvi)) {
    # Generate matrix normal: L_r' %*% Z %*% L_c where Z is n_row x n_col iid
    Z <- matrix(rnorm(n_row_met * n_col_met), nrow = n_row_met)
    corr_mat <- crossprod(L_r, Z) %*% L_c
    # Flatten grid to plot order (column-major: row varies fastest)
    spatial_noise_env[, tt] <- sigma_env[e] * as.vector(corr_mat)
  }

  # --- Variety index lookup ---
  vi <- as.integer(gsub("G", "", design_dt$variety))

  # --- Expand to long format ---
  rbindlist(lapply(seq_len(n_plot_env), function(i) {
    v <- vi[i]

    # Variety NDVI parameters modified by GxE
    ndvi_max_i <- ndvi_max_main[v] + 0.05 * gxe_matrix[v, e]
    decay_i    <- decay_main[v] + 0.1 * gxe_matrix[v, e]
    t_onset_i  <- t_onset_main[v]

    # True NDVI curve: Gaussian decay after t_onset
    ndvi_true <- ifelse(
      time_ndvi <= t_onset_i,
      ndvi_max_i,
      ndvi_max_i * exp(-decay_i * (time_ndvi - t_onset_i)^2 / 1000)
    )

    # Observed NDVI = true + spatial noise
    ndvi_obs <- ndvi_true + spatial_noise_env[i, ]
    # Clamp to [0, 1] for realism
    ndvi_obs <- pmin(pmax(ndvi_obs, 0), 1)

    data.table(
      plot_id     = design_dt$plot_id[i],
      variety     = design_dt$variety[i],
      environment = env_names[e],
      rep         = design_dt$rep[i],
      iblock      = design_dt$iblock[i],
      row         = design_dt$row[i],
      col         = design_dt$col[i],
      time        = time_ndvi,
      ndvi        = round(ndvi_obs, 4)
    )
  }))
}))

# --- Primary trait: yield per environment ---
# Yield = grand mean + variety effect + GxE + environment noise
yield_met <- rbindlist(lapply(seq_len(n_env), function(e) {
  vi <- as.integer(gsub("G", "", sim_met_fda[environment == env_names[e] &
                                               time == time_ndvi[1],
                                             variety]))
  plot_ids <- sim_met_fda[environment == env_names[e] &
                            time == time_ndvi[1], plot_id]
  n_plt <- length(plot_ids)

  # Yield model: intercept + variety main + GxE + noise
  yield_vals <- 3.5 +
    0.8 * (ndvi_max_main[vi] - mean(ndvi_max_main)) +
    0.5 * gxe_matrix[cbind(vi, e)] +
    rnorm(n_plt, 0, 0.3)

  data.table(
    plot_id = plot_ids,
    yield   = round(yield_vals, 2)
  )
}))

sim_met_fda <- merge(sim_met_fda, yield_met, by = "plot_id")

setcolorder(sim_met_fda,
            c("plot_id", "variety", "environment", "rep", "iblock",
              "row", "col", "time", "ndvi", "yield"))
setkey(sim_met_fda, environment, plot_id, time)


# ==============================================================================
# Save datasets
# ==============================================================================

usethis::use_data(sim_grain_fill, sim_met_fda, overwrite = TRUE)

cat("Datasets saved:\n")
cat("  sim_grain_fill:", nrow(sim_grain_fill), "rows x",
    ncol(sim_grain_fill), "cols\n")
cat("  sim_met_fda:   ", nrow(sim_met_fda), "rows x",
    ncol(sim_met_fda), "cols\n")
