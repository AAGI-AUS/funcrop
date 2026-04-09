# Subset method for fda_data

Subsets an `fda_data` object while preserving the class and updating
metadata. Accepts the same arguments as
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
subsetting (row expressions, column selections, etc.).

## Usage

``` r
# S3 method for class 'fda_data'
x[...]
```

## Arguments

- x:

  An `fda_data` object.

- ...:

  Arguments passed to `[.data.table` (e.g., row filter expressions,
  column selections, `by` grouping).

## Value

An `fda_data` object (if essential columns remain) or a `data.table`.
