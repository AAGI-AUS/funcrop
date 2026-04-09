# List Available Estimation Engines

Returns a character vector of estimation engines currently installed and
available for use with funcrop. At least one engine must be installed
for model fitting.

## Usage

``` r
funcrop_engines()
```

## Value

Character vector of available engine names. Possible values are
`"asreml"` (REML via ASReml-R v4.2+) and `"bayesreml"` (Bayesian MCMC
via greta). Returns `character(0)` if no engines are installed.

## Details

funcrop supports two estimation backends:

- `asreml`:

  ASReml-R v4.2+ – restricted maximum likelihood (REML) estimation.
  Commercial licence required. Faster for large models.

- `bayesreml`:

  bayesreml – Bayesian MCMC estimation via greta/TensorFlow.
  Open-source. Provides full posterior distributions and credible
  intervals.

Both packages are listed in Suggests; neither is required at install
time. However, at least one must be installed before fitting any model.

## Examples

``` r
funcrop_engines()
#> [1] "mgcv" "lme4"
```
