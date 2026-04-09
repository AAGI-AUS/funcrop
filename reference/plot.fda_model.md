# Plot method for fda_model

Produces diagnostic and result plots. If `ggplot2` is available,
generates a multi-panel figure showing: (1) the coefficient function
beta(t) with confidence band, and (2) fitted vs observed overlay. Falls
back to base graphics otherwise.

## Usage

``` r
# S3 method for class 'fda_model'
plot(x, which = c("both", "coef", "fitted"), ...)
```

## Arguments

- x:

  An `fda_model` object.

- which:

  Character: `"coef"` for coefficient function, `"fitted"` for fitted
  curves, `"both"` (default).

- ...:

  Additional arguments passed to plotting functions.

## Value

Invisible list of ggplot objects (if ggplot2 available), or NULL.
