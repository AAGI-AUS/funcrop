#' @keywords internal
"_PACKAGE"

#' @importFrom data.table data.table setDT setkey setkeyv copy := .N .SD
#'   .GRP .NGRP fifelse fcase rbindlist setnames setcolorder
#' @importFrom Matrix Matrix sparseMatrix crossprod tcrossprod Diagonal
#'   bdiag t solve
#' @importFrom splines bs splineDesign
#' @importFrom stats model.matrix model.frame terms formula fitted
#'   residuals predict coef vcov confint logLik nobs update anova
#'   quantile median var cov sd dnorm pnorm qnorm setNames
#'   as.formula na.omit complete.cases
#' @importFrom methods is new as
#' @importFrom utils head tail
#' @importFrom graphics par plot lines abline legend points text
#' @importFrom data.table is.data.table uniqueN setattr
#' @importFrom grDevices colorRampPalette
NULL

# Required for data.table := and [.data.table to work inside package code
.datatable.aware <- TRUE

# Suppress R CMD check NOTEs for data.table NSE column references
utils::globalVariables(c(
  "component", "ci_lower", "ci_upper", "id", "level", "coef_idx",
  "variety", "blup", "coef_f", "time", "depth", "pt", "V1",
  ".", "..cols", "environment", "fitted", "se", "group",
  "predicted", "beta", "loading", "score", "env_label",
  "var_label", "value", "metric", "model_label", "trait_type",
  "secondary", "primary", "first_id", "plot_id", "ndvi",
  "grain_weight", "yield", "block", "rep", "iblock",
  "spatial_row", "spatial_col", "primary_trait", "n_ids",
  "mean_value", "sd_value", "model", "engine", "aic", "bic",
  "waic", "rmse", "r_squared", "coverage_95",
  "response", "row_id", ".I", "sire", "dam",
  "backend", ".data", "x", "y",
  "type", "label", "estimate", "env1", "env2", "blup_value"
))
