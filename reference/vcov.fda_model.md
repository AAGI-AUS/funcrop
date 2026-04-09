# Extract variance-covariance information from fda_model

Extract variance-covariance information from fda_model

## Usage

``` r
# S3 method for class 'fda_model'
vcov(object, ...)
```

## Arguments

- object:

  An `fda_model` object.

- ...:

  Additional arguments (ignored).

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
of variance components, or a matrix if available in
`extras$vcov_matrix`.
