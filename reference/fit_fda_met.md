# Fit Functional Data Analysis Model for Multi-Environment Trials

Fits FDA models across multiple environments, accounting for
genotype-by-environment (GxE) interaction on the functional trait and
its relationship to the primary trait of interest. Supports both
two-stage (per-environment profiles then MET model) and single-stage
(joint) approaches, with factor-analytic (FA) variance structures for
modelling GxE.

## Usage

``` r
fit_fda_met(
  data,
  environment_col = "environment",
  time_col = "time",
  value_col = NULL,
  id_col = "id",
  group_col = "variety",
  primary_col = NULL,
  spatial_row_col = NULL,
  spatial_col_col = NULL,
  basis = NULL,
  n_knots = 10L,
  degree = 3L,
  penalty_order = 2L,
  gxe_structure = c("fa1", "fa2", "us", "diag", "compound_symmetry"),
  spatial = c("none", "ar1", "ar1ar1"),
  genomic_matrix = NULL,
  pedigree_matrix = NULL,
  two_stage = TRUE,
  engine = "auto",
  ...
)
```

## Arguments

- data:

  An `fda_data` object or a data.frame/data.table containing
  observations from multiple environments. Must include columns for
  environment, time, value, observational unit ID, and variety/genotype.

- environment_col:

  Character; name of the environment/trial column. Default
  `"environment"`.

- time_col:

  Character; name of the time column. Default `"time"`.

- value_col:

  Character or `NULL`; name of the response column. Default `NULL`
  (auto-detected from `fda_data` or uses `"value"`).

- id_col:

  Character; name of the observational unit column (e.g., plot). Default
  `"id"`.

- group_col:

  Character; name of the variety/genotype column. Default `"variety"`.

- primary_col:

  Character or `NULL`; name of the primary trait column (e.g., yield).
  Default `NULL`.

- spatial_row_col:

  Character or `NULL`; name of the spatial row column. Default `NULL`.

- spatial_col_col:

  Character or `NULL`; name of the spatial column column. Default
  `NULL`.

- basis:

  An `fda_basis` object (from
  [`bspline_basis()`](https://AAGI-AUS.github.io/funcrop/reference/bspline_basis.md)).
  If `NULL` (default), a basis is constructed automatically. A single
  basis is used across all environments (ensuring coefficient
  comparability).

- n_knots:

  Integer; number of internal knots. Default 10.

- degree:

  Integer; B-spline polynomial degree. Default 3.

- penalty_order:

  Integer; difference penalty order. Default 2.

- gxe_structure:

  Character; variance structure for GxE on B-spline coefficients. One of
  `"fa1"` (default), `"fa2"`, `"us"` (unstructured), `"diag"`
  (diagonal/independent environments), or `"compound_symmetry"`. The FA
  structures model \\\mathrm{var}(\alpha_v) = \Lambda \Lambda' + \Psi\\,
  with \\k\\ factors.

- spatial:

  Character; spatial error model per environment. One of `"none"`
  (default), `"ar1"`, or `"ar1ar1"`.

- genomic_matrix:

  Optional numeric matrix (symmetric, PSD) of genomic relationships
  among varieties. Row/colnames must match variety IDs. Passed as
  `vm(variety, source = G)` in the mixed model.

- pedigree_matrix:

  Optional numeric matrix of pedigree-based relationships. Same usage as
  `genomic_matrix`.

- two_stage:

  Logical; if `TRUE` (default), use two-stage approach. If `FALSE`, fit
  a single-stage joint model across all environments.

- engine:

  Character; estimation engine. One of `"auto"` (default), `"asreml"`,
  or `"bayesreml"`.

- ...:

  Additional arguments passed to backend fitting functions or to
  [`fit_functional_profiles()`](https://AAGI-AUS.github.io/funcrop/reference/fit_functional_profiles.md)
  (in two-stage mode).

## Value

An `fda_model` object with the following extras:

- `environment_loadings`:

  Numeric matrix of FA loadings (n_environments x k), where k is the
  number of FA factors.

- `variety_scores`:

  Numeric matrix of FA scores (n_varieties x k).

- `environment_curves`:

  data.table of environment-specific fitted curves per variety, with
  columns `environment`, `variety`, `time`, `fitted`.

- `gxe_variance`:

  data.table decomposing variance into main effect, GxE interaction, and
  residual components per environment.

- `stage1_models`:

  List of per-environment `fda_model` objects (two-stage only).

- `stage1_blups`:

  data.table of per-environment variety spline coefficient BLUPs
  (two-stage only).

## Details

### Two-stage approach

**Stage 1**: For each environment,
[`fit_functional_profiles()`](https://AAGI-AUS.github.io/funcrop/reference/fit_functional_profiles.md)
extracts variety-specific B-spline coefficient BLUPs
\\\hat{\alpha}\_{ve}\\.

**Stage 2**: The coefficients are modelled across environments:
\$\$\hat{\alpha}\_{vek} = \mu_k + g\_{vk} + (ge)\_{vek} +
\varepsilon\_{vek}\$\$

where:

- \\g\_{vk}\\ is the variety main effect on coefficient \\k\\,

- \\(ge)\_{vek}\\ is the GxE interaction, modelled with an FA structure:
  \\\mathrm{var}(\alpha\_{v\cdot}) = \Lambda \Lambda' + \Psi\\,

- \\\Lambda\\ is the loadings matrix (environments x factors),

- and variety scores are the latent factors.

### Single-stage approach

All environments are fitted simultaneously in a single mixed model:

- Fixed: intercept + environment + null-space spline terms

- Random: `variety:spline_coef` with `fa(environment, k)` variance

- Spatial: `at(environment):ar1(row):ar1(col)` (if requested)

### Genomic prediction

When `genomic_matrix` or `pedigree_matrix` is supplied, variety random
effects are modelled as `vm(variety, source = G)`, enabling genomic
prediction of unobserved varieties.

## References

Smith, A.B., Cullis, B.R. and Thompson, R. (2001). Analyzing variety by
environment data using multiplicative mixed models and adjustments for
spatial field trend. *Biometrics*, 57(4), 1138–1147.

Kelly, A.M., Smith, A.B., Eccleston, J.A. and Cullis, B.R. (2007). The
accuracy of varietal selection using factor analytic models for
multi-environment plant breeding trials. *Crop Science*, 47(3),
1063–1070.

De Faveri, J., Verbyla, A.P., Pitchford, W.S., Venkatanagappa, S. and
Cullis, B.R. (2015). Statistical methods for analysis of multi-harvest
data from perennial pasture variety selection trials. *Crop and Pasture
Science*, 66(9), 947–962.

## See also

[`fit_functional_profiles()`](https://AAGI-AUS.github.io/funcrop/reference/fit_functional_profiles.md)
for per-environment profiling,
[`scalar_on_function()`](https://AAGI-AUS.github.io/funcrop/reference/scalar_on_function.md)
for scalar-on-function regression,
[`predict_new_env()`](https://AAGI-AUS.github.io/funcrop/reference/predict_new_env.md)
for predicting into new environments,
[`make_genomic_matrix()`](https://AAGI-AUS.github.io/funcrop/reference/make_genomic_matrix.md),
[`make_pedigree_matrix()`](https://AAGI-AUS.github.io/funcrop/reference/make_pedigree_matrix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulated MET grain-fill data
data(sim_met_grain_fill)

# Two-stage MET-FDA with FA1 GxE structure
met_model <- fit_fda_met(
  data            = sim_met_grain_fill,
  environment_col = "site",
  time_col        = "das",
  value_col       = "ndvi",
  id_col          = "plot",
  group_col       = "variety",
  n_knots         = 8,
  gxe_structure   = "fa1",
  spatial         = "ar1ar1",
  spatial_row_col = "row",
  spatial_col_col = "col",
  two_stage       = TRUE,
  engine          = "auto"
)

# Examine GxE decomposition
met_model$extras$gxe_variance
met_model$extras$environment_loadings
} # }
```
