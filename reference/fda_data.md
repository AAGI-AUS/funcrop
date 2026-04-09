# Create a functional data object

Constructs a validated `fda_data` object wrapping functional
observations (e.g., repeated measurements of a secondary trait over
time) into a structured
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with metadata. Designed for crop variety trial data where secondary
functional traits (NDVI, canopy temperature, grain-fill rate) are
related to a primary trait (e.g., yield).

## Usage

``` r
fda_data(
  time,
  value,
  id,
  group = NULL,
  spatial_row = NULL,
  spatial_col = NULL,
  primary_trait = NULL,
  primary_trait_name = NULL,
  ...
)
```

## Arguments

- time:

  Numeric vector of time points (or other continuous dimension such as
  thermal time, depth, wavelength).

- value:

  Numeric vector of measured values (the secondary / functional trait,
  e.g., NDVI).

- id:

  Character or factor identifying the observational unit (e.g., plot,
  plant).

- group:

  Optional character or factor for grouping (e.g., variety, genotype).
  Default `NULL`.

- spatial_row:

  Optional integer vector of spatial row positions. Default `NULL`.

- spatial_col:

  Optional integer vector of spatial column positions. Default `NULL`.

- primary_trait:

  Optional numeric vector of primary trait values (e.g., yield). Must
  contain exactly one value per unique `id`. Recycled internally to
  match the observation-level data. Default `NULL`.

- primary_trait_name:

  Optional character string naming the primary trait (e.g.,
  `"yield_t_ha"`). Default `NULL`.

- ...:

  Additional named vectors of the same length as `time`, stored as extra
  columns in the underlying data.table.

## Value

An object of class `fda_data`, inheriting from `data.table` and
`data.frame`. Contains an attribute `"fda_meta"` with metadata (trait
names, id column, etc.).

## Validation

- `time`, `value`, and `id` must be non-NA and of equal length.

- `time` and `value` must be numeric and finite.

- If supplied, `primary_trait` must have exactly one value per unique
  `id`.

- Spatial columns must be integer-valued.

## Examples

``` r
# Simulate simple functional data
set.seed(42)
ids <- rep(paste0("plot_", 1:5), each = 10)
times <- rep(seq(100, 190, by = 10), times = 5)
vals <- sin(times / 50) + rnorm(50, sd = 0.1)
grps <- rep(c("var_A", "var_B", "var_A", "var_B", "var_A"), each = 10)
yields <- rep(c(3.2, 4.1, 3.8, 4.5, 3.0), each = 10)

fd <- fda_data(
  time = times,
  value = vals,
  id = ids,
  group = grps,
  primary_trait = yields,
  primary_trait_name = "yield_t_ha"
)
print(fd)
#> -- fda_data ---------------------------------
#>   Observations : 50 rows, 5 unique ids
#>   Time range   : [100.00, 190.00] (10 unique points)
#>   Value range  : [-0.6866, 1.0464]
#>   Groups       : 2 (var_A, var_B)
#>   Primary trait: yield_t_ha
#> ---------------------------------------------
#> Key: <id, time>
#>         id  time       value  group primary_trait
#>     <fctr> <num>       <num> <fctr>         <num>
#>  1: plot_1   100  1.04639327  var_A           3.2
#>  2: plot_1   110  0.75202659  var_A           3.2
#>  3: plot_1   120  0.71177602  var_A           3.2
#>  4: plot_1   130  0.57878763  var_A           3.2
#>  5: plot_1   140  0.37541498  var_A           3.2
#>  6: plot_1   150  0.13050756  var_A           3.2
#>  7: plot_1   160  0.09277806  var_A           3.2
#>  8: plot_1   170 -0.26500701  var_A           3.2
#>  9: plot_1   180 -0.24067807  var_A           3.2
#> 10: plot_1   190 -0.61812930  var_A           3.2
#> 11: plot_2   100  1.03978439  var_B           4.1
#> 12: plot_2   110  1.03716094  var_B           4.1
#> 13: plot_2   120  0.53657711  var_B           4.1
#> 14: plot_2   130  0.48762250  var_B           4.1
#> 15: plot_2   140  0.32165602  var_B           4.1
#> 16: plot_2   150  0.20471505  var_B           4.1
#> 17: plot_2   160 -0.08679944  var_B           4.1
#> 18: plot_2   170 -0.52118664  var_B           4.1
#> 19: plot_2   180 -0.68656714  var_B           4.1
#> 20: plot_2   190 -0.47984656  var_B           4.1
#> 21: plot_3   100  0.87863357  var_A           3.8
#> 22: plot_3   110  0.63036556  var_A           3.8
#> 23: plot_3   120  0.65827144  var_A           3.8
#> 24: plot_3   130  0.63696884  var_A           3.8
#> 25: plot_3   140  0.52450750  var_A           3.8
#> 26: plot_3   150  0.09807309  var_A           3.8
#> 27: plot_3   160 -0.08410108  var_A           3.8
#> 28: plot_3   170 -0.43185741  var_A           3.8
#> 29: plot_3   180 -0.39651071  var_A           3.8
#> 30: plot_3   190 -0.67585738  var_A           3.8
#> 31: plot_4   100  0.95484244  var_B           4.5
#> 32: plot_4   110  0.87898014  var_B           4.5
#> 33: plot_4   120  0.77897353  var_B           4.5
#> 34: plot_4   130  0.45460873  var_B           4.5
#> 35: plot_4   140  0.38548366  var_B           4.5
#> 36: plot_4   150 -0.03058086  var_B           4.5
#> 37: plot_4   160 -0.13682004  var_B           4.5
#> 38: plot_4   170 -0.34063186  var_B           4.5
#> 39: plot_4   180 -0.68394121  var_B           4.5
#> 40: plot_4   190 -0.60824563  var_B           4.5
#> 41: plot_5   100  0.92989729  var_A           3.0
#> 42: plot_5   110  0.77239067  var_A           3.0
#> 43: plot_5   120  0.75127950  var_A           3.0
#> 44: plot_5   130  0.44283089  var_A           3.0
#> 45: plot_5   140  0.19816005  var_A           3.0
#> 46: plot_5   150  0.18440181  var_A           3.0
#> 47: plot_5   160 -0.13951346  var_A           3.0
#> 48: plot_5   170 -0.11113098  var_A           3.0
#> 49: plot_5   180 -0.48566506  var_A           3.0
#> 50: plot_5   190 -0.54629310  var_A           3.0
#>         id  time       value  group primary_trait
#>     <fctr> <num>       <num> <fctr>         <num>
summary(fd)
#> == fda_data summary ==========================
#> 
#> Observations : 50
#> Unique ids   : 5
#> Time points  : 10 unique in [100.00, 190.00]
#> Obs per id   : median = 10, range = [10, 10]
#> 
#> Value (secondary trait):
#>   Mean = 0.1981, SD = 0.5460
#>   Range = [-0.6866, 1.0464]
#> 
#> Groups: 2 levels
#>   Per-group summary:
#>   group n_ids mean_value  sd_value
#>  <fctr> <int>      <num>     <num>
#>   var_A     3  0.2132907 0.5178221
#>   var_B     2  0.1752893 0.5988588
#> 
#> Primary trait (yield_t_ha):
#>   Mean = 3.7200, SD = 0.6221, Range = [3.0000, 4.5000]
#> 
#> ==============================================
```
