# Simulated Grain-Fill Trial Data

A simulated randomised complete block design (RCBD) trial with 20
varieties across 3 blocks (60 plots), measuring grain weight over 8 time
points during grain filling. Includes spatially correlated residuals
(AR1 in rows) and a primary yield trait related to grain-fill
parameters.

## Usage

``` r
sim_grain_fill
```

## Format

A `data.table` with 480 rows (60 plots x 8 time points) and 8 columns:

- plot_id:

  Character. Unique plot identifier (e.g., `"P001"`).

- variety:

  Character. Variety identifier (`"V01"` to `"V20"`).

- block:

  Character. Block identifier (`"B1"` to `"B3"`).

- row:

  Integer. Row position in the spatial layout (1–12).

- col:

  Integer. Column position in the spatial layout (1–5).

- time:

  Numeric. Days after anthesis (10, 15, 20, 25, 30, 35, 40, 45).

- grain_weight:

  Numeric. Observed grain weight (mg) at each time point, including
  block effects and spatially correlated noise.

- yield:

  Numeric. Plot-level yield (t/ha), related to underlying grain-fill
  parameters via \\5 + 0.1 W\_{\max} + 0.5 \times \text{rate} \times
  100 + \epsilon\\.

## Source

Simulated data. See `data-raw/simulate_data.R` for the full generation
code. Seed: `set.seed(20250415)`.

## Details

The true grain-fill curves follow a logistic growth model: \$\$w(t) =
W\_{\max} / (1 + \exp(-\text{rate} \times (t - t\_{\text{mid}})))\$\$

The spatial layout is a 12 x 5 grid. Residual noise has AR1 correlation
along rows (\\\phi = 0.3\\) with marginal standard deviation 1.5. Block
effects are drawn from \\N(0, 2^2)\\.

## Examples

``` r
data(sim_grain_fill)
str(sim_grain_fill)
#> Classes ‘data.table’ and 'data.frame':   480 obs. of  8 variables:
#>  $ variety     : chr  "V06" "V06" "V06" "V06" ...
#>  $ plot_id     : chr  "P001" "P001" "P001" "P001" ...
#>  $ block       : chr  "B1" "B1" "B1" "B1" ...
#>  $ row         : int  1 1 1 1 1 1 1 1 2 2 ...
#>  $ col         : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ time        : num  10 15 20 25 30 35 40 45 10 15 ...
#>  $ grain_weight: num  6.34 8.56 10.07 17.55 23.97 ...
#>  $ yield       : num  15.5 15.5 15.5 15.5 15.5 ...
#>  - attr(*, ".internal.selfref")=<externalptr> 

# Number of plots and time points
sim_grain_fill[, .N, by = plot_id]
#>     plot_id     N
#>      <char> <int>
#>  1:    P001     8
#>  2:    P002     8
#>  3:    P003     8
#>  4:    P004     8
#>  5:    P005     8
#>  6:    P006     8
#>  7:    P007     8
#>  8:    P008     8
#>  9:    P009     8
#> 10:    P010     8
#> 11:    P011     8
#> 12:    P012     8
#> 13:    P013     8
#> 14:    P014     8
#> 15:    P015     8
#> 16:    P016     8
#> 17:    P017     8
#> 18:    P018     8
#> 19:    P019     8
#> 20:    P020     8
#> 21:    P021     8
#> 22:    P022     8
#> 23:    P023     8
#> 24:    P024     8
#> 25:    P025     8
#> 26:    P026     8
#> 27:    P027     8
#> 28:    P028     8
#> 29:    P029     8
#> 30:    P030     8
#> 31:    P031     8
#> 32:    P032     8
#> 33:    P033     8
#> 34:    P034     8
#> 35:    P035     8
#> 36:    P036     8
#> 37:    P037     8
#> 38:    P038     8
#> 39:    P039     8
#> 40:    P040     8
#> 41:    P041     8
#> 42:    P042     8
#> 43:    P043     8
#> 44:    P044     8
#> 45:    P045     8
#> 46:    P046     8
#> 47:    P047     8
#> 48:    P048     8
#> 49:    P049     8
#> 50:    P050     8
#> 51:    P051     8
#> 52:    P052     8
#> 53:    P053     8
#> 54:    P054     8
#> 55:    P055     8
#> 56:    P056     8
#> 57:    P057     8
#> 58:    P058     8
#> 59:    P059     8
#> 60:    P060     8
#>     plot_id     N
#>      <char> <int>
sim_grain_fill[, unique(time)]
#> [1] 10 15 20 25 30 35 40 45
```
