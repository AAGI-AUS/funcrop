# Coerce fda_data to data.table

Coerce fda_data to data.table

## Usage

``` r
as.data.table.fda_data(x, keep.rownames = FALSE, ...)
```

## Arguments

- x:

  An `fda_data` object.

- keep.rownames:

  Ignored (for compatibility).

- ...:

  Additional arguments (ignored).

## Value

A plain
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
without `fda_data` class or metadata.
