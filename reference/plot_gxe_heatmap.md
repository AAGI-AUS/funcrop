# Heatmap of Genotype-by-Environment Interaction

Visualises the variety x environment interaction for functional traits
using a heatmap. Can display variance component estimates, FA loadings,
or predicted functional summary statistics across environments.

## Usage

``` r
plot_gxe_heatmap(
  model,
  value_col = c("predicted", "variance", "correlation"),
  ...
)
```

## Arguments

- model:

  An `fda_model` object from a multi-environment trial (MET) analysis.

- value_col:

  Character; which quantity to display. One of `"predicted"` (default,
  variety x environment predicted values), `"variance"` (variance
  components), or `"correlation"` (pairwise environment correlations).

- ...:

  Additional arguments (currently unused).

## Value

A `ggplot` object (invisibly), or `NULL` for base R fallback.

## Examples

``` r
if (FALSE) { # \dontrun{
met_model <- fit_functional_profiles(sim_met_fda, ...)
plot_gxe_heatmap(met_model)
} # }
```
