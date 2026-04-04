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
#' @importFrom methods is new
#' @importFrom graphics par plot lines abline legend
#' @importFrom data.table is.data.table uniqueN setattr
#' @importFrom grDevices colorRampPalette
NULL

# Required for data.table := and [.data.table to work inside package code
.datatable.aware <- TRUE
