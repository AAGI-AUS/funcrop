# Fit 2D Functional Profiles via Tensor Product B-Splines

Models a functional trait observed over two continuous indices (e.g.,
time and soil depth, time and wavelength) using tensor product B-splines
with separable penalties. Group-specific (e.g., variety-specific) 2D
surfaces are estimated as random effects.

## Usage

``` r
fit_2d_functional(
  data,
  time_col = "time",
  depth_col = "depth",
  value_col = "value",
  id_col = "id",
  group_col = "variety",
  n_knots1 = 10L,
  n_knots2 = 10L,
  degree1 = 3L,
  degree2 = 3L,
  spatial_row_col = NULL,
  spatial_col_col = NULL,
  engine = "auto",
  ...
)
```

## Arguments

- data:

  A `data.frame` or
  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  in long format with one row per observation (id x dim1 x dim2).

- time_col:

  Character; column name for the first functional dimension (e.g.,
  time). Default `"time"`.

- depth_col:

  Character; column name for the second functional dimension (e.g.,
  depth, wavelength). Default `"depth"`.

- value_col:

  Character; column name for the response. Default `"value"`.

- id_col:

  Character; column name for the observational unit. Default `"id"`.

- group_col:

  Character; column name for the grouping factor (e.g., variety).
  Default `"variety"`.

- n_knots1:

  Integer; number of internal knots for the first dimension. Default 10.

- n_knots2:

  Integer; number of internal knots for the second dimension. Default
  10.

- degree1:

  Integer; B-spline degree for dimension 1. Default 3.

- degree2:

  Integer; B-spline degree for dimension 2. Default 3.

- spatial_row_col:

  Character or `NULL`; spatial row column. Default `NULL`.

- spatial_col_col:

  Character or `NULL`; spatial column column. Default `NULL`.

- engine:

  Character; estimation engine. Default `"auto"`.

- ...:

  Additional arguments passed to the backend.

## Value

An `fda_model` object with `model_type = "2d_functional"` in extras. The
`fitted_curves` component contains fitted 2D surfaces per group, and
`basis` contains the `fda_tensor_basis` object.

## Statistical Model

\$\$y(t, d) = f_v(t, d) + \text{spatial} + \varepsilon\$\$ where
\\f_v(t, d) = \sum\_{j} \sum\_{k} \alpha\_{vjk} B\_{1j}(t) B\_{2k}(d)\\
is a variety-specific tensor product random surface, and the penalty is
the standard Kronecker sum: \\P = P_1 \otimes I_2 + I_1 \otimes P_2\\.

## References

Marx, B.D. and Eilers, P.H.C. (2005). Multidimensional penalized signal
regression. *Technometrics*, 47(1), 13–22.

Wood, S.N. (2006). Low-rank scale-invariant tensor product smooths for
generalized additive mixed models. *Biometrics*, 62(4), 1025–1036.

## See also

[`tensor_bspline_basis()`](https://AAGI-AUS.github.io/funcrop/reference/tensor_bspline_basis.md)
for tensor product basis construction,
[`scalar_on_2d_function()`](https://AAGI-AUS.github.io/funcrop/reference/scalar_on_2d_function.md)
for relating 2D profiles to a scalar response.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulated 2D functional data: time x depth
set.seed(123)
n_obs <- 500
dt2d <- data.table::data.table(
  id      = rep(paste0("P", 1:10), each = 50),
  variety = rep(rep(c("V1", "V2"), each = 5), each = 50),
  time    = runif(n_obs, 0, 100),
  depth   = runif(n_obs, 0, 1),
  value   = rnorm(n_obs)
)

fit2d <- fit_2d_functional(
  data      = dt2d,
  time_col  = "time",
  depth_col = "depth",
  value_col = "value",
  id_col    = "id",
  group_col = "variety",
  n_knots1  = 6, n_knots2 = 5
)
} # }
```
