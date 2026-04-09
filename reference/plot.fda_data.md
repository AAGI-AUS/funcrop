# Plot method for fda_data (spaghetti plot)

Produces a spaghetti plot of functional curves. If `ggplot2` is
available, uses it with a colourblind-safe viridis palette; otherwise
falls back to base graphics.

## Usage

``` r
# S3 method for class 'fda_data'
plot(x, max_curves = 50L, alpha = 0.4, ...)
```

## Arguments

- x:

  An `fda_data` object.

- max_curves:

  Maximum number of individual curves to draw (default 50). Set to `Inf`
  to plot all.

- alpha:

  Line transparency (default 0.4).

- ...:

  Additional arguments passed to the plotting function.

## Value

Invisibly returns the ggplot object (if ggplot2 is available) or `NULL`.
