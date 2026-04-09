# Compare REML and Bayesian Model Fits

Produces a two-panel comparison of results from ASReml-R (REML) and
bayesreml (Bayesian) backends. The left panel compares fitted functional
curves; the right panel compares coefficient function estimates. REML
results are shown as point estimates with standard error bands, while
Bayesian results show the posterior mean with credible intervals.

## Usage

``` r
plot_backend_comparison(model_reml, model_bayes, varieties = NULL, ...)
```

## Arguments

- model_reml:

  An `fda_model` object fitted with the ASReml backend
  (`engine = "asreml"`).

- model_bayes:

  An `fda_model` object fitted with the bayesreml backend
  (`engine = "bayesreml"`).

- varieties:

  Character vector of variety names to include in the curves panel, or
  `NULL` (default) for all. Plots may be cluttered if many varieties are
  included; consider limiting to 5–10.

- ...:

  Additional arguments (currently unused).

## Value

A `ggplot` object (invisibly), or `NULL` for base R fallback. When
ggplot2 is available, both panels are arranged side by side.

## Details

Colour coding:

- ASReml (REML): solid lines with dashed CI bands

- bayesreml (Bayesian): semi-transparent ribbons for credible intervals

## Examples

``` r
if (FALSE) { # \dontrun{
model_a <- fit_functional_profiles(..., engine = "asreml")
model_b <- fit_functional_profiles(..., engine = "bayesreml")
plot_backend_comparison(model_a, model_b)
} # }
```
