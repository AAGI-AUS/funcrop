# Plot Fitted Functional Profiles

Produces publication-ready plots of variety-specific fitted functional
curves from an `fda_model` object. Optionally includes
confidence/credible interval ribbons. When ggplot2 is available, uses
colourblind-safe viridis palettes and a clean minimal theme; otherwise
falls back to base R graphics.

## Usage

``` r
plot_functional_profiles(
  model,
  varieties = NULL,
  ci = TRUE,
  n_grid = 200L,
  alpha = 0.2,
  line_size = 0.8,
  ...
)
```

## Arguments

- model:

  An `fda_model` object (e.g., from
  [`fit_functional_profiles()`](https://biometryhub.github.io/funcrop/reference/fit_functional_profiles.md),
  [`fit_fda_joint()`](https://biometryhub.github.io/funcrop/reference/fit_fda_joint.md),
  or
  [`fit_2d_functional()`](https://biometryhub.github.io/funcrop/reference/fit_2d_functional.md)).

- varieties:

  Character vector of variety/group names to plot, or `NULL` (default)
  to plot all varieties. Useful for reducing clutter when many varieties
  are present.

- ci:

  Logical; if `TRUE` (default), add confidence/credible interval ribbons
  around each fitted curve. Requires `se` or `ci_lower`/`ci_upper`
  columns in the fitted curves data.

- n_grid:

  Integer; number of grid points for smooth curve evaluation.
  Default 200. Higher values give smoother curves at marginal cost.

- alpha:

  Numeric; ribbon transparency (0–1). Default 0.2.

- line_size:

  Numeric; line width for fitted curves. Default 0.8.

- ...:

  Additional arguments (currently unused).

## Value

A `ggplot` object (invisibly), or `NULL` if ggplot2 is unavailable (base
R plot is drawn directly).

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_functional_profiles(sim_grain_fill, ...)
plot_functional_profiles(model)
plot_functional_profiles(model, varieties = c("V01", "V05", "V10"))
} # }
```
