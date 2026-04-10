# Predict Variety Performance in New Environments

Uses the fitted FA loadings and variety scores from a MET-FDA model to
predict functional profiles and primary trait performance in new or
unobserved environments.

## Usage

``` r
predict_new_env(
  model,
  new_environment_loadings,
  varieties = NULL,
  n_grid = 200L,
  ci = TRUE,
  ...
)
```

## Arguments

- model:

  An `fda_model` object from
  [`fit_fda_met()`](https://biometryhub.github.io/funcrop/reference/fit_fda_met.md).
  Must contain `environment_loadings` and `variety_scores` in its
  extras.

- new_environment_loadings:

  Numeric vector of FA loadings for the new environment. Length must
  equal the number of FA factors (columns in the loadings matrix). These
  can be estimated externally (e.g., from environmental covariates) or
  interpolated from observed environments.

- varieties:

  Character vector of variety IDs to predict for. Default `NULL`
  predicts all varieties in the fitted model.

- n_grid:

  Integer; number of time points for predicted curves. Default 200.

- ci:

  Logical; if `TRUE` (default), compute approximate confidence intervals
  for predicted curves. Requires the specific variance (Psi) from the FA
  model.

- ...:

  Additional arguments (currently unused).

## Value

A data.table with columns:

- variety:

  Variety identifier.

- time:

  Evaluation time point.

- predicted:

  Predicted functional value at this time point.

- se:

  Approximate standard error (if `ci = TRUE`).

- ci_lower:

  Lower 95 percent confidence bound (if `ci = TRUE`).

- ci_upper:

  Upper 95 percent confidence bound (if `ci = TRUE`).

## Details

The prediction mechanism leverages the factor-analytic decomposition.
For a new environment with loadings \\\lambda\_{new}\\, the predicted
B-spline coefficients for variety \\v\\ are:

\$\$\hat{\alpha}\_{v,new} = \mu + \lambda\_{new} \cdot s_v\$\$

where \\s_v\\ is the variety's estimated FA score vector and \\\mu\\ is
the overall mean coefficient vector. The predicted curve is then
obtained by evaluating the B-spline basis at the predicted coefficients.

Uncertainty is approximated from the specific variance \\\Psi\_{new}\\
(extrapolated) and the prediction error variance of the scores.

## References

Kelly, A.M., Smith, A.B., Eccleston, J.A. and Cullis, B.R. (2007). The
accuracy of varietal selection using factor analytic models for
multi-environment plant breeding trials. *Crop Science*, 47(3),
1063–1070.

## See also

[`fit_fda_met()`](https://biometryhub.github.io/funcrop/reference/fit_fda_met.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After fitting a MET model
met_model <- fit_fda_met(...)

# Predict into a new environment with estimated loadings
new_preds <- predict_new_env(
  model = met_model,
  new_environment_loadings = c(0.85),  # FA1 loading for new site
  varieties = c("var_A", "var_B", "var_C")
)
} # }
```
