# Scalar-on-2D-Function Regression

Extends scalar-on-function regression to functional predictors observed
over two continuous dimensions (e.g., time x depth). The model relates a
scalar primary trait to a 2D functional predictor via double
integration:

## Usage

``` r
scalar_on_2d_function(
  primary_trait,
  functional_profiles_2d,
  basis = NULL,
  genomic_matrix = NULL,
  engine = "auto",
  ...
)
```

## Arguments

- primary_trait:

  A named numeric vector of primary trait values (one per group/variety)
  or a
  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  with columns `group` and `value`.

- functional_profiles_2d:

  An `fda_model` object from
  [`fit_2d_functional()`](https://AAGI-AUS.github.io/funcrop/reference/fit_2d_functional.md),
  or a list containing `fitted_surfaces` (a
  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  with columns `group`, `time`, `depth`, `fitted`) and `basis` (an
  `fda_tensor_basis` object).

- basis:

  An `fda_tensor_basis` object. If `NULL`, extracted from
  `functional_profiles_2d`.

- genomic_matrix:

  Optional genomic relationship matrix. Default `NULL`.

- engine:

  Character; estimation engine. Default `"auto"`.

- ...:

  Additional arguments passed to the backend.

## Value

An `fda_model` object where:

- coefficient_function:

  Contains the estimated 2D coefficient surface \\\beta(t, d)\\ as a
  data.table with columns `time`, `depth`, `beta`, `se`, `ci_lower`,
  `ci_upper`.

- predictions:

  Predicted primary trait values for each group.

## Details

\$\$y_i = \alpha + \int_t \int_d \beta(t, d) f_i(t, d)\\ dt\\ dd +
\varepsilon_i\$\$

where \\\beta(t, d)\\ is a 2D coefficient surface estimated via tensor
product B-splines with penalised regression, and \\f_i(t, d)\\ are the
group-specific 2D functional profiles.

## References

Marx, B.D. and Eilers, P.H.C. (2005). Multidimensional penalized signal
regression. *Technometrics*, 47(1), 13–22.

## See also

[`fit_2d_functional()`](https://AAGI-AUS.github.io/funcrop/reference/fit_2d_functional.md)
for 2D profile estimation,
[`scalar_on_function()`](https://AAGI-AUS.github.io/funcrop/reference/scalar_on_function.md)
for the 1D equivalent.

## Examples

``` r
if (FALSE) { # \dontrun{
# After fitting 2D profiles:
# fit2d <- fit_2d_functional(...)
# primary_vals <- c(V1 = 4.2, V2 = 3.8, ...)
# sof2d <- scalar_on_2d_function(primary_vals, fit2d)
} # }
```
