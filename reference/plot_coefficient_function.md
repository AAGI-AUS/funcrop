# Plot the Coefficient Function

Plots the estimated coefficient function \\\beta(t)\\ from a
scalar-on-function model (or \\\beta(t, d)\\ heatmap/contour from a
scalar-on-2D-function model). Includes confidence/credible interval
ribbons and a horizontal reference line at \\\beta = 0\\.

## Usage

``` r
plot_coefficient_function(model, ci = TRUE, n_grid = 200L, ...)
```

## Arguments

- model:

  An `fda_model` object from
  [`scalar_on_function()`](https://AAGI-AUS.github.io/funcrop/reference/scalar_on_function.md)
  or
  [`scalar_on_2d_function()`](https://AAGI-AUS.github.io/funcrop/reference/scalar_on_2d_function.md).

- ci:

  Logical; if `TRUE` (default), add confidence/credible interval ribbons
  (1D) or a separate uncertainty panel (2D).

- n_grid:

  Integer; number of grid points for evaluation. Default 200 (1D) or 50
  per dimension (2D).

- ...:

  Additional arguments (currently unused).

## Value

A `ggplot` object (invisibly), or `NULL` for base R fallback.

## Details

For 1D models, the plot shows \\\beta(t)\\ as a line with a ribbon for
the CI/credible interval and a dashed horizontal line at zero. For 2D
models (from
[`scalar_on_2d_function()`](https://AAGI-AUS.github.io/funcrop/reference/scalar_on_2d_function.md)),
a heatmap of \\\beta(t, d)\\ is drawn using
[`ggplot2::geom_tile()`](https://ggplot2.tidyverse.org/reference/geom_tile.html)
with a viridis continuous fill scale, overlaid with contour lines.

## Examples

``` r
if (FALSE) { # \dontrun{
sof_model <- scalar_on_function(yields, profiles)
plot_coefficient_function(sof_model)
} # }
```
