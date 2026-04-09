# Fit Variety-Specific Functional Profiles (Stage 1)

Fits a spatial + temporal linear mixed model per trial to extract
variety-specific functional profiles (smooth curves) from repeated
measurements of a secondary trait. This is the first stage of the
two-stage FDA approach: Stage 1 recovers the shape of each variety's
functional response; Stage 2 (via
[`scalar_on_function()`](https://AAGI-AUS.github.io/funcrop/reference/scalar_on_function.md))
regresses a primary trait (e.g., yield) onto these profiles.

## Usage

``` r
fit_functional_profiles(
  data,
  time_col = "time",
  value_col = NULL,
  id_col = "id",
  group_col = "variety",
  spatial_row_col = NULL,
  spatial_col_col = NULL,
  basis = NULL,
  n_knots = 10L,
  degree = 3L,
  penalty_order = 2L,
  spatial = c("none", "ar1", "ar1ar1", "spline"),
  engine = "auto",
  ...
)
```

## Arguments

- data:

  An `fda_data` object (from
  [`fda_data()`](https://AAGI-AUS.github.io/funcrop/reference/fda_data.md))
  or a data.frame / data.table containing the observations. If a raw
  data.frame, the column mapping arguments below must be specified.

- time_col:

  Character; name of the time column in `data`. Default `"time"`.
  Ignored when `data` is an `fda_data` object.

- value_col:

  Character or `NULL`; name of the response (value) column. Default
  `NULL`, which uses `"value"` for `fda_data` or requires explicit
  specification for raw data.

- id_col:

  Character; name of the observational unit identifier column (e.g.,
  plot). Default `"id"`.

- group_col:

  Character; name of the variety/genotype grouping column. Default
  `"variety"`.

- spatial_row_col:

  Character or `NULL`; name of the spatial row column. Default `NULL`
  (no spatial modelling unless detected from `fda_data`).

- spatial_col_col:

  Character or `NULL`; name of the spatial column column. Default
  `NULL`.

- basis:

  An `fda_basis` object (from
  [`bspline_basis()`](https://AAGI-AUS.github.io/funcrop/reference/bspline_basis.md)).
  If `NULL` (default), a basis is constructed automatically from the
  unique time points in `data` using `n_knots`, `degree`, and
  `penalty_order`.

- n_knots:

  Integer; number of internal knots for automatic basis construction.
  Default 10. Ignored if `basis` is supplied.

- degree:

  Integer; B-spline polynomial degree. Default 3 (cubic). Ignored if
  `basis` is supplied.

- penalty_order:

  Integer; order of the difference penalty for the P-spline. Default 2
  (penalises curvature). Ignored if `basis` is supplied.

- spatial:

  Character; type of spatial error model. One of `"none"` (default),
  `"ar1"` (first-order autoregressive in rows), `"ar1ar1"` (separable
  AR1 x AR1 in rows and columns), or `"spline"` (2D spatial spline).
  Requires `spatial_row_col` and `spatial_col_col` to be specified
  (except `"none"`).

- engine:

  Character; estimation engine. One of `"auto"` (default), `"asreml"`,
  or `"bayesreml"`. See
  [`funcrop_engines()`](https://AAGI-AUS.github.io/funcrop/reference/funcrop_engines.md).

- ...:

  Additional arguments passed to the backend fitting function (e.g.,
  `control` for ASReml, `mcmc_control` for bayesreml).

## Value

An `fda_model` object (see `new_fda_model()`) containing:

- fitted_curves:

  data.table of variety-specific fitted curves on a fine time grid, with
  columns `id`, `time`, `fitted`, `se`.

- coefficient_function:

  Empty list (not applicable for Stage 1; the coefficient function
  beta(t) is estimated in Stage 2).

- variance_components:

  data.table of estimated variance components.

- predictions:

  data.table of variety-level predicted values.

- residuals:

  Model residuals.

- basis:

  The `fda_basis` object used.

- data:

  The input data (as `fda_data`).

- engine:

  Character string identifying the engine used.

- call:

  The matched function call.

- extras:

  List with `spline_blups` (variety x basis coefficient matrix),
  `spline_decomposition` (from
  [`make_Zspline()`](https://AAGI-AUS.github.io/funcrop/reference/make_Zspline.md)),
  `raw_model` (backend model object), and `convergence` diagnostics.

## Details

### Model formulation

For each trial / environment, the model is: \$\$y\_{ijk}(t) = f_v(t) +
b_j + s(r, c) + \varepsilon\_{ijk}(t)\$\$

where:

- \\f_v(t) = \sum_k \alpha\_{vk} B_k(t)\\ is the variety-specific curve,
  with \\\alpha_v \sim N(0, \sigma_v^2 I)\\ as random B-spline
  coefficients (P-spline mixed model representation).

- \\b_j\\ is a block effect (if `block` column is present).

- \\s(r, c)\\ is an optional spatial error model (AR1, AR1xAR1, or
  spline).

- \\\varepsilon\_{ijk}(t) \sim N(0, \sigma^2)\\ is the residual.

The variety-specific B-spline random effects are structured as a
Kronecker product: \\Z = I_V \otimes Z\_{spline}\\, where
\\Z\_{spline}\\ comes from
[`make_Zspline()`](https://AAGI-AUS.github.io/funcrop/reference/make_Zspline.md)
and \\I_V\\ is the identity over varieties. This yields \\V \times K\\
random effect levels (V varieties, K spline coefficients).

### Backend dispatch

For ASReml: uses [`str()`](https://rdrr.io/r/utils/str.html)
specification from `.asreml_build_spline_str()` with penalised variance
structure. For bayesreml: uses known matrices from
`.bayesreml_build_spline_str()` with penalty-as-precision prior.

## References

Verbyla, A.P., Cavanagh, C.R. and Verbyla, K.L. (2012). Whole-genome
analysis of multienvironment or multitrait QTL in MAGIC. *G3: Genes,
Genomes, Genetics*, 2(9), 1085–1093.

De Faveri, J., Verbyla, A.P., Pitchford, W.S., Venkatanagappa, S. and
Cullis, B.R. (2015). Statistical methods for analysis of multi-harvest
data from perennial pasture variety selection trials. *Crop and Pasture
Science*, 66(9), 947–962.

## See also

[`scalar_on_function()`](https://AAGI-AUS.github.io/funcrop/reference/scalar_on_function.md)
for Stage 2,
[`bspline_basis()`](https://AAGI-AUS.github.io/funcrop/reference/bspline_basis.md)
for basis construction,
[`make_Zspline()`](https://AAGI-AUS.github.io/funcrop/reference/make_Zspline.md)
for mixed model reparameterisation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load simulated grain-fill data
data(sim_grain_fill)

# Create fda_data object
fd <- fda_data(
  time  = sim_grain_fill$time,
  value = sim_grain_fill$grain_weight,
  id    = sim_grain_fill$plot_id,
  group = sim_grain_fill$variety,
  spatial_row = sim_grain_fill$row,
  spatial_col = sim_grain_fill$col,
  primary_trait = sim_grain_fill$yield,
  primary_trait_name = "yield_t_ha"
)

# Fit functional profiles (Stage 1)
profiles <- fit_functional_profiles(
  data    = fd,
  n_knots = 6,
  spatial = "ar1ar1",
  engine  = "auto"
)

print(profiles)
plot(profiles, which = "fitted")

# Extract BLUPs for Stage 2
blup_matrix <- profiles$extras$spline_blups
} # }
```
