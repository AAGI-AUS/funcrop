# Choosing and Using Backends

## Overview

funcrop v0.2.0 supports four estimation backends:

| Feature            |     mgcv     |      lme4      |    ASReml    |   bayesreml    |
|--------------------|:------------:|:--------------:|:------------:|:--------------:|
| Estimation         |  REML (GAM)  |   REML (LMM)   |     REML     | Bayesian MCMC  |
| Speed              |     Fast     |      Fast      |     Fast     | Slow (minutes) |
| Licence            |  Open (GPL)  |   Open (GPL)   |  Commercial  |   Open (GPL)   |
| Availability       | Ships with R |      CRAN      | VSNi licence |  CRAN + greta  |
| Uncertainty        |  Approx. SE  |   Approx. SE   |  Approx. SE  | Full posterior |
| AR1 residuals      |  Via gamm()  |       No       |     Yes      |      Yes       |
| FA covariance      |      No      |       No       |     Yes      |    Partial     |
| Genomic G-matrix   |      No      | Cholesky trick |     vm()     |      vm()      |
| Tensor products    | Native te()  |       No       |    Manual    |     Manual     |
| GLMM               |     Yes      |      Yes       |      No      |      Yes       |
| Large data (\>50K) |    bam()     |      Yes       |     Yes      |      Slow      |

## Check Available Engines

``` r
library(funcrop)
cat("Available engines:", paste(funcrop_engines(), collapse = ", "), "\n")
#> Available engines: mgcv, lme4
cat("Default engine:", funcrop_default_engine(), "\n")
#> Default engine: mgcv
```

## Example: Baseline RCBD Model

We compare both backends on a simple variety trial:

``` r
library(data.table)
data(sim_grain_fill)
dt <- copy(sim_grain_fill)
dt[, variety := factor(variety)]
dt[, block := factor(block)]

# Add plot-level yield noise
set.seed(42)
pids <- unique(dt$plot_id)
dt <- merge(dt,
            data.table(plot_id = pids,
                       yield_noise = rnorm(length(pids), 0, 0.8)),
            by = "plot_id", sort = FALSE)
dt[, yield_plot := yield + yield_noise]
yp <- unique(dt[, .(plot_id, variety, block, yield_plot)])
```

### ASReml

``` r
library(asreml)
m_asr <- asreml(yield_plot ~ block, random = ~ variety,
                data = yp, trace = FALSE)
cat("Converged:", m_asr$converge, "\n\n")
cat("Variance components:\n")
print(summary(m_asr)$varcomp)
```

    #> ASReml not available. Install from https://vsni.co.uk/software/asreml-r

### bayesreml

``` r
library(bayesreml)
m_bay <- bayesreml(yield_plot ~ block, random = ~ variety,
                   data = copy(yp),  # copy() prevents by-ref modification
                   n_samples = 500, warmup = 250, chains = 2,
                   verbose = FALSE, mcmc_verbose = FALSE)
cat("Variance components (posterior):\n")
print(m_bay$extras$variance_comps)
```

    #> bayesreml not available. See the bayesreml package documentation.

## Comparison

``` r
# ASReml: variances
asr_sigv <- summary(m_asr)$varcomp["variety", "component"]
asr_sige <- summary(m_asr)$varcomp["units!R", "component"]

# bayesreml: SDs (square for variance comparison)
bay_vc <- m_bay$extras$variance_comps
bay_sigv_sd <- bay_vc$estimate[grep("variety", bay_vc$component)][1]
bay_sige_sd <- bay_vc$estimate[grep("sigma_e", bay_vc$component)][1]

comp <- data.frame(
  Parameter = c("sigma_variety^2", "sigma_residual^2"),
  ASReml = round(c(asr_sigv, asr_sige), 4),
  bayesreml_var = round(c(bay_sigv_sd^2, bay_sige_sd^2), 4),
  pct_diff = round(100 * (c(bay_sigv_sd^2, bay_sige_sd^2) -
                           c(asr_sigv, asr_sige)) /
                     c(asr_sigv, asr_sige), 1)
)
cat("\nVariance component comparison:\n")
print(comp, row.names = FALSE)
```

## Expected Differences

1.  **Variance components**: bayesreml reports SDs with lognormal
    priors; REML gives point estimates that can be zero (boundary).
    After squaring bayesreml SDs, expect 5–15% differences for small
    datasets.

2.  **BLUPs vs posterior means**: Rankings are typically identical (r \>
    0.99). Bayesian estimates are slightly more shrunk due to prior.

3.  **Convergence**: ASReml may fail for complex models (AI
    singularity); bayesreml with priors often succeeds where REML fails.

4.  **Computational cost**: ASReml is 10–100x faster. Use ASReml for
    exploration, bayesreml for final inference with full uncertainty.

5.  **bayesreml modifies data.table by reference**: Always use
    [`copy()`](https://rdrr.io/pkg/data.table/man/copy.html) when
    passing data.table objects to bayesreml.

## Recommendation

- **Routine analysis**: Use ASReml (fast, well-tested)
- **Publication-quality inference**: Run bayesreml for credible
  intervals
- **Model diagnostics**: Compare both – divergence signals weak
  identification or prior sensitivity
- **No ASReml licence**: bayesreml provides full functionality (slower)
