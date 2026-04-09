# Simulated Multi-Environment Trial with Functional NDVI Trait

A simulated multi-environment trial (MET) with 30 varieties across 4
environments, using an alpha-lattice design (3 replicates, 5 incomplete
blocks per replicate). NDVI (normalised difference vegetation index) is
measured at 6 time points as a functional stay-green trait. Includes
genotype-by-environment (GxE) interaction via a factor-analytic (FA1)
structure, AR1 x AR1 spatial correlation, and environment-specific
residual variances.

## Usage

``` r
sim_met_fda
```

## Format

A `data.table` with 2160 rows (30 varieties x 4 environments x 3 reps x
6 time points) and 10 columns:

- plot_id:

  Character. Unique plot identifier within each environment (e.g.,
  `"E1_P001"`).

- variety:

  Character. Variety identifier (`"G01"` to `"G30"`).

- environment:

  Character. Environment identifier (`"E1"` to `"E4"`).

- rep:

  Character. Replicate identifier (`"R1"` to `"R3"`).

- iblock:

  Character. Incomplete block identifier (`"IB1"` to `"IB5"`).

- row:

  Integer. Row position in the spatial layout (1–15).

- col:

  Integer. Column position in the spatial layout (1–6).

- time:

  Numeric. Days after sowing (70, 80, 90, 100, 110, 120).

- ndvi:

  Numeric. Observed NDVI value at each time point, clamped to the range
  0 to 1.

- yield:

  Numeric. Plot-level yield (t/ha), incorporating variety main effects,
  GxE interaction, and environment-specific noise.

## Source

Simulated data. See `data-raw/simulate_data.R` for the full generation
code. Seed: `set.seed(20250416)`.

## Details

The true NDVI curves follow a Gaussian decay model: \$\$\text{NDVI}(t) =
\text{NDVI}\_{\max} \exp\left(-\frac{\text{decay} \times (t -
t\_{\text{onset}})^2}{1000}\right)\$\$ for \\t \> t\_{\text{onset}}\\,
and \\\text{NDVI}\_{\max}\\ otherwise.

The GxE structure uses a factor-analytic model with 1 factor:
environment loadings \\(0.8, 0.5, -0.3, -0.7)\\ and variety scores
\\\sim N(0, 1)\\. Spatial correlation is AR1 x AR1 with
\\\phi\_{\text{row}} = 0.4\\, \\\phi\_{\text{col}} = 0.3\\ per
environment. Residual standard deviations vary by environment: \\(0.03,
0.04, 0.035, 0.05)\\.

The spatial layout per environment is a 15 x 6 grid (90 plots per
environment, 360 plots total).

## Examples

``` r
data(sim_met_fda)
str(sim_met_fda)
#> Classes ‘data.table’ and 'data.frame':   2160 obs. of  10 variables:
#>  $ variety    : chr  "G29" "G29" "G29" "G29" ...
#>  $ environment: chr  "E1" "E1" "E1" "E1" ...
#>  $ plot_id    : chr  "E1_P001" "E1_P001" "E1_P001" "E1_P001" ...
#>  $ rep        : chr  "R1" "R1" "R1" "R1" ...
#>  $ iblock     : chr  "IB1" "IB1" "IB1" "IB1" ...
#>  $ row        : num  1 1 1 1 1 1 1 1 1 1 ...
#>  $ col        : num  1 1 1 1 1 1 2 2 2 2 ...
#>  $ time       : num  70 80 90 100 110 120 70 80 90 100 ...
#>  $ ndvi       : num  0.761 0.737 0.696 0.568 0.503 ...
#>  $ yield      : num  1.61 1.61 1.61 1.61 1.61 ...
#>  - attr(*, ".internal.selfref")=<externalptr> 

# Plots per environment
sim_met_fda[, data.table::uniqueN(plot_id), by = environment]
#>    environment    V1
#>         <char> <int>
#> 1:          E1    90
#> 2:          E2    90
#> 3:          E3    90
#> 4:          E4    90

# Time points
sim_met_fda[, unique(time)]
#> [1]  70  80  90 100 110 120
```
