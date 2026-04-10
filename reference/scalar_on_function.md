# Scalar-on-Function Regression (Stage 2)

Fits a scalar-on-function regression model that relates a scalar primary
trait (e.g., yield) to the variety-specific functional profiles
estimated in Stage 1 (via
[`fit_functional_profiles()`](https://biometryhub.github.io/funcrop/reference/fit_functional_profiles.md)).
The key output is the coefficient function \\\beta(t)\\, which
identifies the time periods during which the secondary trait most
strongly predicts the primary trait.

## Usage

``` r
scalar_on_function(
  primary_trait,
  functional_profiles,
  basis = NULL,
  genomic_matrix = NULL,
  pedigree_matrix = NULL,
  environment_col = NULL,
  engine = "auto",
  ...
)
```

## Arguments

- primary_trait:

  One of:

  - A named numeric vector with names corresponding to variety
    identifiers.

  - An `fda_data` object containing a `primary_trait` column (extracted
    automatically).

  - A data.frame / data.table with columns `variety` and `yield` (or the
    name specified via the first column).

- functional_profiles:

  One of:

  - An `fda_model` object from
    [`fit_functional_profiles()`](https://biometryhub.github.io/funcrop/reference/fit_functional_profiles.md)
    (preferred).

  - A numeric matrix of B-spline coefficients (varieties x basis
    functions), in which case `basis` must also be supplied.

- basis:

  An `fda_basis` object (from
  [`bspline_basis()`](https://biometryhub.github.io/funcrop/reference/bspline_basis.md)).
  Required when `functional_profiles` is a raw matrix; ignored when it
  is an `fda_model` (the basis is extracted from the model object).

- genomic_matrix:

  Optional numeric matrix (symmetric, PSD) of genomic relationships
  among varieties. Row/column names must match variety identifiers. Used
  as `vm(variety, source = G)` in the mixed model.

- pedigree_matrix:

  Optional numeric matrix (symmetric, PSD) of pedigree-based
  relationships. Same naming convention as `genomic_matrix`.

- environment_col:

  Character or `NULL`; name of the environment column in `primary_trait`
  if multi-environment data is supplied. Default `NULL`
  (single-environment analysis).

- engine:

  Character; estimation engine. One of `"auto"` (default), `"asreml"`,
  or `"bayesreml"`. See
  [`funcrop_engines()`](https://biometryhub.github.io/funcrop/reference/funcrop_engines.md).

- ...:

  Additional arguments passed to the backend fitting function.

## Value

An `fda_model` object containing:

- fitted_curves:

  Empty data.table (not applicable for Stage 2).

- coefficient_function:

  List with `time`, `beta`, `se`, `ci_lower`, `ci_upper` – the estimated
  coefficient function \\\hat{\beta}(t)\\ evaluated on a fine grid.

- variance_components:

  data.table of estimated variance components.

- predictions:

  data.table of predicted primary trait values per variety.

- residuals:

  Model residuals.

- basis:

  The `fda_basis` object used.

- data:

  The input `primary_trait` data.

- engine:

  Character string identifying the engine used.

- call:

  The matched function call.

- extras:

  List with `beta_coefs` (estimated B-spline coefficients for
  \\\beta(t)\\), `J_matrix` (inner product matrix of B-spline basis),
  `C_matrix` (functional covariate matrix), `raw_model`.

## Details

### Model formulation

\$\$y_i = \alpha + \int \beta(t) \hat{f}\_i(t) \\ dt + u_i +
\varepsilon_i\$\$

where:

- \\\hat{f}\_i(t)\\ is the Stage 1 estimated functional profile for
  variety \\i\\.

- \\\beta(t) = \sum_k b_k B_k(t)\\ is the coefficient function, expanded
  in the same B-spline basis.

- The integral is approximated using the inner product matrix:
  \$\$J\_{jk} = \int B_j(t) B_k(t) \\ dt\$\$ computed via numerical
  quadrature (composite Simpson's rule).

- The functional covariate matrix is: \$\$C = \hat{A} J\$\$ where
  \\\hat{A}\\ is the variety x basis coefficient matrix from Stage 1.
  Each row of \\C\\ is the functional covariate for one variety.

- \\u_i\\ is an optional genomic/pedigree random effect: \\u \sim N(0,
  \sigma_g^2 G)\\ where \\G\\ is the relationship matrix.

- \\\varepsilon_i \sim N(0, \sigma^2)\\ is the residual.

### Smoothness of beta(t)

The coefficient function is regularised via a P-spline penalty on the
B-spline coefficients \\b\\: the penalty matrix from the basis is
incorporated as a random effect precision structure (mixed model
representation of P-splines).

## References

Ramsay, J.O. and Silverman, B.W. (2005). *Functional Data Analysis* (2nd
ed.). Springer.

Reiss, P.T. and Ogden, R.T. (2007). Functional principal component
regression and functional partial least squares. *Journal of the
American Statistical Association*, 102(479), 984–996.

Marx, B.D. and Eilers, P.H.C. (1999). Generalized linear regression on
sampled signals and curves: a P-spline approach. *Technometrics*, 41(1),
1–13.

## See also

[`fit_functional_profiles()`](https://biometryhub.github.io/funcrop/reference/fit_functional_profiles.md)
for Stage 1,
[`compare_methods()`](https://biometryhub.github.io/funcrop/reference/compare_methods.md)
for model comparison,
[`bspline_basis()`](https://biometryhub.github.io/funcrop/reference/bspline_basis.md)
for basis construction.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming `profiles` is an fda_model from fit_functional_profiles()
data(sim_grain_fill)

# Named yield vector
yield_vec <- setNames(
  sim_grain_fill[, .(yield = mean(yield)), by = variety]$yield,
  sim_grain_fill[, .(yield = mean(yield)), by = variety]$variety
)

# Stage 2: scalar-on-function regression
sof_fit <- scalar_on_function(
  primary_trait        = yield_vec,
  functional_profiles  = profiles,
  engine               = "auto"
)

# Inspect the coefficient function
plot(sof_fit, which = "coef")
coef(sof_fit)
} # }
```
