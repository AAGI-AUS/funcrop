# Factor-Analytic Biplot for GxE Interaction

Produces a biplot from the factor-analytic (FA) decomposition of
genotype-by-environment interaction effects estimated from B-spline
coefficients. Environment loadings and variety scores are overlaid in a
two-dimensional space showing the dominant patterns of GxE interaction.

## Usage

``` r
plot_fa_biplot(
  model,
  factors = c(1L, 2L),
  point_size = 3,
  text_size = 3.5,
  ...
)
```

## Arguments

- model:

  An `fda_model` object containing FA structure in extras (e.g.,
  `model$extras$fa_loadings` and `model$extras$fa_scores`), or with
  variance components that can be decomposed.

- factors:

  Integer vector of length 2 specifying which FA factors to plot.
  Default `c(1, 2)`.

- point_size:

  Numeric; size of points. Default 3.

- text_size:

  Numeric; size of labels. Default 3.5.

- ...:

  Additional arguments (currently unused).

## Value

A `ggplot` object (invisibly), or `NULL` for base R fallback.

## Details

The biplot displays:

- **Environment loadings** as points/arrows, showing how environments
  differentiate varieties.

- **Variety scores** as labelled points, showing relative variety
  performance across the FA factors.

If the `ggrepel` package is available, labels are automatically
positioned to avoid overlapping. Otherwise,
[`ggplot2::geom_text()`](https://ggplot2.tidyverse.org/reference/geom_text.html)
is used.

## Examples

``` r
if (FALSE) { # \dontrun{
met_model <- fit_functional_profiles(sim_met_fda, ...)
plot_fa_biplot(met_model)
} # }
```
