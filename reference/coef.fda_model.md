# Extract coefficient function beta(t) from fda_model

Extract coefficient function beta(t) from fda_model

## Usage

``` r
# S3 method for class 'fda_model'
coef(object, ...)
```

## Arguments

- object:

  An `fda_model` object.

- ...:

  Additional arguments (ignored).

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with columns `time`, `beta`, `se`, `ci_lower`, `ci_upper`.
