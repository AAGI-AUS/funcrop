# funcrop <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/AAGI-AUS/funcrop/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AAGI-AUS/funcrop/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**Functional Data Analysis for Crop Variety Trials**

funcrop integrates B-spline basis functions with linear mixed models to
analyse functional traits (grain-fill curves, NDVI profiles, stay-green
trajectories) in crop variety trials. It supports scalar-on-function
regression, 2D functional surfaces, and multi-environment trial (MET)
analyses with genotype-by-environment (GxE) interaction.

## Key Features

- **Scalar-on-function regression** -- relate yield to the entire
  grain-fill trajectory via a coefficient function beta(t)
- **Multi-environment FDA** -- FA1/FA2/US/diagonal GxE structures on
  functional coefficients
- **Genomic prediction** -- VanRaden G-matrices, pedigree A-matrices,
  blended H-matrices
- **4 estimation backends** -- mgcv, lme4, ASReml-R, bayesreml

## Installation

```r
# From GitHub
remotes::install_github("AAGI-AUS/funcrop")
```

At least one backend must be installed:

| Backend | Install | Notes |
|---------|---------|-------|
| mgcv | Ships with R | Default. Open-source. |
| lme4 | `install.packages("lme4")` | Open-source. GLMM support. |
| ASReml-R | [vsni.co.uk](https://vsni.co.uk/software/asreml-r) | Commercial licence. FA, AR1, genomic vm(). |
| bayesreml | `install.packages("bayesreml")` | Open-source. Full posterior inference. |

## Quick Start

```r
library(funcrop)
data(sim_grain_fill)

# Stage 1: Fit variety-specific grain-fill profiles
profiles <- fit_functional_profiles(
  data      = sim_grain_fill,
  time_col  = "time",
  value_col = "grain_weight",
  id_col    = "plot_id",
  group_col = "variety",
  n_knots   = 6,
  engine    = "auto"
)

# Stage 2: Scalar-on-function regression (yield ~ grain-fill curve)
yield_vec <- setNames(
  sim_grain_fill[, .(y = mean(yield)), by = variety]$y,
  sim_grain_fill[, .(y = mean(yield)), by = variety]$variety
)

sof <- scalar_on_function(
  primary_trait       = yield_vec,
  functional_profiles = profiles,
  engine              = "auto"
)

# Inspect the coefficient function beta(t)
plot(sof, which = "coef")
```

## Backend Selection

funcrop auto-detects the best available backend (mgcv > lme4 > ASReml > bayesreml).
Override with the `engine` argument or set a session default:

```r
# Check available backends
funcrop_engines()

# Set default for the session
funcrop_default_engine("asreml")

# Or specify per call
profiles <- fit_functional_profiles(..., engine = "lme4")
```

See `vignette("choosing-backends")` for a detailed comparison.

## Core Functions

| Function | Purpose |
|----------|---------|
| `fit_functional_profiles()` | Stage 1: variety-specific functional curves |
| `scalar_on_function()` | Stage 2: coefficient function beta(t) |
| `fit_fda_joint()` | Single-stage joint scalar-functional model |
| `fit_2D_functional()` | 2D functional surfaces (time x depth) |
| `fit_fda_met()` | Multi-environment FDA with GxE |
| `bspline_basis()` | B-spline basis construction |
| `make_Zspline()` | Mixed model reparameterisation |
| `make_genomic_matrix()` | VanRaden genomic relationship matrix |

## Mathematical Foundation

The package implements P-spline mixed models where the smoothing penalty
is estimated as a variance ratio. For variety g at time t:

f_g(t) = X_0(t)' beta + Z(t)' alpha_g,  alpha_g ~ N(0, sigma^2_u I)

where X_0 spans the unpenalised null space (population mean curve) and Z
spans the penalised range space (variety-specific smooth deviations).

A full mathematical specification (24 pages, LaTeX) is available in the
project repository as `funcrop_mathematical_specification.pdf`.

## Citation

```r
citation("funcrop")
```

## Licence

GPL (>= 3)

## Acknowledgements

Funded by the [Analytics for the Australian Grains Industry (AAGI)](https://www.aagi.com.au/) project.
Developed at the University of Adelaide.
