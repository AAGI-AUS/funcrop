# Print Method for fda_basis Objects

Displays a compact summary of a B-spline basis object.

## Usage

``` r
# S3 method for class 'fda_basis'
print(x, ...)
```

## Arguments

- x:

  An object of class `"fda_basis"`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the input object `x`.

## Examples

``` r
basis <- bspline_basis(seq(0, 1, length.out = 50), n_knots = 8)
print(basis)
#> B-spline basis (fda_basis)
#>   Degree:            3 
#>   Internal knots:   8 (equally_spaced)
#>   Basis functions:   12 
#>   Eval. points:      50 
#>   Boundary:         [-0.01, 1.01]
#>   Penalty order:     2 
#>   Basis matrix dim: 50 x 12
```
