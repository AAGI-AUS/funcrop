# Fit a Joint Model for Primary and Functional Secondary Traits

Fits a single-stage model that simultaneously models a scalar primary
trait (e.g., yield) and a functional secondary trait (e.g., grain-fill
over time), allowing the functional trait to inform estimation of the
primary trait. This avoids the error propagation issues inherent in the
two-stage approach where functional profiles are first estimated, then
used as predictors.

## Usage

``` r
fit_fda_joint(
  data,
  primary_col = "yield",
  time_col = "time",
  secondary_col = NULL,
  id_col = "id",
  group_col = "variety",
  spatial_row_col = NULL,
  spatial_col_col = NULL,
  basis = NULL,
  n_knots = 10L,
  degree = 3L,
  penalty_order = 2L,
  spatial = c("none", "ar1", "ar1ar1"),
  genomic_matrix = NULL,
  pedigree_matrix = NULL,
  engine = "auto",
  ...
)
```

## Arguments

- data:

  A `data.frame` or
  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  in long format containing both the functional secondary trait
  observations (one row per id x time combination) and the primary trait
  (constant within each id). May also be an `fda_data` object.

- primary_col:

  Character; column name for the scalar primary trait (e.g., yield).
  Default `"yield"`.

- time_col:

  Character; column name for the time/index variable of the functional
  trait. Default `"time"`.

- secondary_col:

  Character or `NULL`; column name for the functional secondary trait
  values. If `NULL` (default), inferred as the first numeric column that
  is not `primary_col`, `time_col`, or a spatial column.

- id_col:

  Character; column name identifying observational units (e.g., plot).
  Default `"id"`.

- group_col:

  Character; column name for the grouping factor (e.g., variety,
  genotype). Default `"variety"`.

- spatial_row_col:

  Character or `NULL`; column name for spatial row coordinates. Default
  `NULL` (no spatial modelling).

- spatial_col_col:

  Character or `NULL`; column name for spatial column coordinates.
  Default `NULL`.

- basis:

  An `fda_basis` object (from
  [`bspline_basis()`](https://AAGI-AUS.github.io/funcrop/reference/bspline_basis.md)).
  If `NULL` (default), a basis is constructed from the data using
  `n_knots`, `degree`, and `penalty_order`.

- n_knots:

  Integer; number of internal B-spline knots if `basis` is `NULL`.
  Default 10.

- degree:

  Integer; B-spline polynomial degree. Default 3 (cubic).

- penalty_order:

  Integer; order of the difference penalty. Default 2.

- spatial:

  Character; spatial correlation structure for residuals. One of
  `"none"` (default), `"ar1"` (AR1 in rows), or `"ar1ar1"` (separable
  AR1 x AR1 in rows and columns).

- genomic_matrix:

  Optional numeric matrix (square, symmetric) of genomic relationships
  among groups. Row/column names must match levels of `group_col`.
  Default `NULL`.

- pedigree_matrix:

  Optional numeric matrix of pedigree-based relationships. Same naming
  requirements as `genomic_matrix`. Default `NULL`.

- engine:

  Character; estimation engine. One of `"auto"` (default), `"asreml"`,
  or `"bayesreml"`. See
  [`funcrop_engines()`](https://AAGI-AUS.github.io/funcrop/reference/funcrop_engines.md).

- ...:

  Additional arguments passed to the backend fitting function (e.g.,
  `maxiter`, `workspace` for ASReml; `mcmc_control` for bayesreml).

## Value

An `fda_model` object containing:

- fitted_curves:

  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  of variety-specific fitted functional curves for the secondary trait.

- coefficient_function:

  List with the implied relationship \\\beta(t)\\ between the secondary
  functional trait and the primary scalar trait.

- variance_components:

  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  of estimated variance components including the cross-trait
  correlations.

- predictions:

  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  of predicted primary trait values per group.

- basis:

  The `fda_basis` object used for the functional trait.

- engine:

  Character string identifying the backend used.

## Statistical Model

The joint model treats the secondary trait at each time point and the
primary trait as a multivariate response. In stacked form:

\$\$ \begin{bmatrix} y\_{\text{secondary}}(t) \\ y\_{\text{primary}}
\end{bmatrix} = \begin{bmatrix} f_v(t) + \text{spatial} +
\varepsilon_s(t) \\ \alpha + \int \beta(t) f_v(t)\\ dt + \varepsilon_p
\end{bmatrix} \$\$

where \\f_v(t) = \sum_k \alpha\_{vk} B_k(t)\\ is the variety-specific
B-spline random effect for the secondary trait, and the primary trait is
linked through the functional integral \\\int \beta(t) f_v(t)\\ dt\\.

The model is fitted as a bivariate mixed model with:

- A trait indicator variable distinguishing secondary vs primary
  observations

- B-spline basis x variety interactions for the secondary trait

- Numerical integration weights for the functional covariate in the
  primary trait equation

- Correlated random effects across traits via `us(trait)` structures

## Implementation

Data are stacked with a `trait_type` factor (`"secondary"`,
`"primary"`). For secondary-trait rows, the design matrix uses B-spline
basis x variety interactions. For primary-trait rows, the design matrix
uses the numerical integral of the B-spline basis functions (trapezoidal
rule), creating a functional covariate. Variety-level random effects are
correlated across both trait types.

## References

De Boor, C. (2001). *A Practical Guide to Splines* (Revised ed.).
Springer.

Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines
and penalties. *Statistical Science*, 11(2), 89–121.

## See also

[`fit_functional_profiles()`](https://AAGI-AUS.github.io/funcrop/reference/fit_functional_profiles.md)
for Stage 1 of the two-stage approach,
[`scalar_on_function()`](https://AAGI-AUS.github.io/funcrop/reference/scalar_on_function.md)
for Stage 2,
[`bspline_basis()`](https://AAGI-AUS.github.io/funcrop/reference/bspline_basis.md)
for basis construction.

## Examples

``` r
if (FALSE) { # \dontrun{
data(sim_grain_fill)

# Fit the joint model
joint_fit <- fit_fda_joint(
  data          = sim_grain_fill,
  primary_col   = "yield",
  time_col      = "time",
  secondary_col = "grain_weight",
  id_col        = "plot_id",
  group_col     = "variety",
  spatial_row_col = "row",
  spatial_col_col = "col",
  n_knots       = 6,
  spatial       = "ar1ar1",
  engine        = "asreml"
)
print(joint_fit)
summary(joint_fit)
} # }
```
