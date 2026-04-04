# data.R -- Documentation for simulated datasets in funcrop
#
# Roxygen2 documentation for package datasets.
# See data-raw/simulate_data.R for generation code.

#' Simulated Grain-Fill Trial Data
#'
#' A simulated randomised complete block design (RCBD) trial with 20 varieties
#' across 3 blocks (60 plots), measuring grain weight over 8 time points during
#' grain filling. Includes spatially correlated residuals (AR1 in rows) and a
#' primary yield trait related to grain-fill parameters.
#'
#' The true grain-fill curves follow a logistic growth model:
#' \deqn{w(t) = W_{\max} / (1 + \exp(-\text{rate} \times (t - t_{\text{mid}})))}
#'
#' @format A \code{data.table} with 480 rows (60 plots x 8 time points) and
#'   8 columns:
#' \describe{
#'   \item{plot_id}{Character. Unique plot identifier (e.g., \code{"P001"}).}
#'   \item{variety}{Character. Variety identifier (\code{"V01"} to
#'     \code{"V20"}).}
#'   \item{block}{Character. Block identifier (\code{"B1"} to \code{"B3"}).}
#'   \item{row}{Integer. Row position in the spatial layout (1--12).}
#'   \item{col}{Integer. Column position in the spatial layout (1--5).}
#'   \item{time}{Numeric. Days after anthesis (10, 15, 20, 25, 30, 35, 40,
#'     45).}
#'   \item{grain_weight}{Numeric. Observed grain weight (mg) at each time
#'     point, including block effects and spatially correlated noise.}
#'   \item{yield}{Numeric. Plot-level yield (t/ha), related to underlying
#'     grain-fill parameters via
#'     \eqn{5 + 0.1 W_{\max} + 0.5 \times \text{rate} \times 100 + \epsilon}.}
#' }
#'
#' @details
#' The spatial layout is a 12 x 5 grid. Residual noise has AR1 correlation
#' along rows (\eqn{\phi = 0.3}) with marginal standard deviation 1.5.
#' Block effects are drawn from \eqn{N(0, 2^2)}.
#'
#' @source Simulated data. See \code{data-raw/simulate_data.R} for the full
#'   generation code. Seed: \code{set.seed(20250415)}.
#'
#' @examples
#' data(sim_grain_fill)
#' str(sim_grain_fill)
#'
#' # Number of plots and time points
#' sim_grain_fill[, .N, by = plot_id]
#' sim_grain_fill[, unique(time)]
"sim_grain_fill"


#' Simulated Multi-Environment Trial with Functional NDVI Trait
#'
#' A simulated multi-environment trial (MET) with 30 varieties across 4
#' environments, using an alpha-lattice design (3 replicates, 5 incomplete
#' blocks per replicate). NDVI (normalised difference vegetation index) is
#' measured at 6 time points as a functional stay-green trait. Includes
#' genotype-by-environment (GxE) interaction via a factor-analytic (FA1)
#' structure, AR1 x AR1 spatial correlation, and environment-specific residual
#' variances.
#'
#' The true NDVI curves follow a Gaussian decay model:
#' \deqn{\text{NDVI}(t) = \text{NDVI}_{\max} \exp\left(-\frac{\text{decay}
#'   \times (t - t_{\text{onset}})^2}{1000}\right)}
#' for \eqn{t > t_{\text{onset}}}, and \eqn{\text{NDVI}_{\max}} otherwise.
#'
#' @format A \code{data.table} with 2160 rows (30 varieties x 4 environments x
#'   3 reps x 6 time points) and 10 columns:
#' \describe{
#'   \item{plot_id}{Character. Unique plot identifier within each environment
#'     (e.g., \code{"E1_P001"}).}
#'   \item{variety}{Character. Variety identifier (\code{"G01"} to
#'     \code{"G30"}).}
#'   \item{environment}{Character. Environment identifier (\code{"E1"} to
#'     \code{"E4"}).}
#'   \item{rep}{Character. Replicate identifier (\code{"R1"} to \code{"R3"}).}
#'   \item{iblock}{Character. Incomplete block identifier (\code{"IB1"} to
#'     \code{"IB5"}).}
#'   \item{row}{Integer. Row position in the spatial layout (1--15).}
#'   \item{col}{Integer. Column position in the spatial layout (1--6).}
#'   \item{time}{Numeric. Days after sowing (70, 80, 90, 100, 110, 120).}
#'   \item{ndvi}{Numeric. Observed NDVI value at each time point, clamped to
#'     the range 0 to 1.}
#'   \item{yield}{Numeric. Plot-level yield (t/ha), incorporating variety main
#'     effects, GxE interaction, and environment-specific noise.}
#' }
#'
#' @details
#' The GxE structure uses a factor-analytic model with 1 factor:
#' environment loadings \eqn{(0.8, 0.5, -0.3, -0.7)} and variety scores
#' \eqn{\sim N(0, 1)}. Spatial correlation is AR1 x AR1 with
#' \eqn{\phi_{\text{row}} = 0.4}, \eqn{\phi_{\text{col}} = 0.3} per
#' environment. Residual standard deviations vary by environment:
#' \eqn{(0.03, 0.04, 0.035, 0.05)}.
#'
#' The spatial layout per environment is a 15 x 6 grid (90 plots per
#' environment, 360 plots total).
#'
#' @source Simulated data. See \code{data-raw/simulate_data.R} for the full
#'   generation code. Seed: \code{set.seed(20250416)}.
#'
#' @examples
#' data(sim_met_fda)
#' str(sim_met_fda)
#'
#' # Plots per environment
#' sim_met_fda[, data.table::uniqueN(plot_id), by = environment]
#'
#' # Time points
#' sim_met_fda[, unique(time)]
"sim_met_fda"
