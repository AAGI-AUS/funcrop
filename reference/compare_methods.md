# Compare Multiple FDA Model Fits

Takes multiple `fda_model` objects and computes comparison metrics
including information criteria (AIC, BIC, WAIC where available),
prediction accuracy (RMSE, R-squared), and coverage of
confidence/credible intervals.

## Usage

``` r
compare_methods(..., labels = NULL)
```

## Arguments

- ...:

  Two or more `fda_model` objects to compare.

- labels:

  Optional character vector of model labels. If `NULL`, labels are
  generated from the call or assigned sequentially.

## Value

An `fda_comparison` object (see `new_fda_comparison()`) containing:

- models:

  Named list of the input `fda_model` objects.

- metrics:

  data.table of comparison metrics with one row per model and columns:
  `model`, `engine`, `aic`, `bic`, `waic`, `rmse`, `r_squared`,
  `coverage_95`.

- label:

  Character description of the comparison.

## Details

Metrics computed:

- AIC, BIC:

  Extracted from the fitted model where available (ASReml). `NA` for
  Bayesian models (use WAIC instead).

- WAIC:

  Widely applicable information criterion. Available for bayesreml
  models only.

- RMSE:

  Root mean squared error of predictions vs observed primary trait
  values.

- R-squared:

  Coefficient of determination for the primary trait predictions.

- Coverage:

  Proportion of observed values falling within the 95\\
  confidence/credible intervals of the fitted curves (where available).

## See also

[`fit_functional_profiles()`](https://biometryhub.github.io/funcrop/reference/fit_functional_profiles.md),
[`scalar_on_function()`](https://biometryhub.github.io/funcrop/reference/scalar_on_function.md),
`new_fda_comparison()`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare two models with different basis sizes
comp <- compare_methods(model_10k, model_20k,
                        labels = c("10 knots", "20 knots"))
print(comp)
summary(comp)
plot(comp)
} # }
```
