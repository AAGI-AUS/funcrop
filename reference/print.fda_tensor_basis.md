# Print Method for fda_tensor_basis Objects

Displays a compact summary of a tensor product B-spline basis.

## Usage

``` r
# S3 method for class 'fda_tensor_basis'
print(x, ...)
```

## Arguments

- x:

  An object of class `"fda_tensor_basis"`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the input object `x`.

## Examples

``` r
set.seed(42)
n <- 50
tb <- tensor_bspline_basis(runif(n), runif(n), n_knots1 = 5, n_knots2 = 4)
print(tb)
#> Tensor product B-spline basis (fda_tensor_basis)
#>   Dimension 1:
#>     Degree:          3 
#>     Internal knots: 5 (equally_spaced)
#>     Basis functions: 9 
#>     Boundary:       [-0.005901, 0.9987]
#>     Penalty order:   2 
#>   Dimension 2:
#>     Degree:          3 
#>     Internal knots: 4 (equally_spaced)
#>     Basis functions: 8 
#>     Boundary:       [-0.009587, 0.9926]
#>     Penalty order:   2 
#>   Total basis functions:72 (9 x 8)
#>   Observations:         50 
#>   Basis matrix dim:    50 x 72
```
