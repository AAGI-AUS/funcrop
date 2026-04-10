# Construct a Difference Penalty Matrix

Builds the penalty matrix \\P = D_d^T D_d\\ where \\D_d\\ is the d-th
order difference matrix of dimension `(n_basis - d) x n_basis`. This is
the standard discrete penalty used in P-spline smoothing (Eilers & Marx,
1996).

## Usage

``` r
make_penalty(n_basis, order = 2L)
```

## Arguments

- n_basis:

  Integer; number of basis functions (columns of the B-spline basis
  matrix). Must be \>= 2.

- order:

  Integer; order of the difference penalty. Default is 2 (second-order,
  penalising curvature). `order = 1` corresponds to a random walk
  penalty. Must satisfy `1 <= order < n_basis`.

## Value

A sparse penalty matrix of class `dgCMatrix` (from the Matrix package),
dimension `n_basis x n_basis`. The matrix is symmetric and positive
semi-definite with rank `n_basis - order`.

## Details

For `order = 1`, the difference matrix \\D_1\\ has entries \\(D_1)\_{ij}
= 1\\ if \\j = i+1\\ and \\-1\\ if \\j = i\\. Higher-order difference
matrices are obtained by recursive application: \\D_d = D_1^{(n-d+1)}
D\_{d-1}\\, where the superscript denotes the dimension of the
first-order matrix.

## References

Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines
and penalties. *Statistical Science*, 11(2), 89–121.

## See also

[`bspline_basis()`](https://biometryhub.github.io/funcrop/reference/bspline_basis.md)
for basis construction,
[`make_Zspline()`](https://biometryhub.github.io/funcrop/reference/make_Zspline.md)
for mixed model reparameterisation.

## Examples

``` r
# Second-order penalty for 12 basis functions
P <- make_penalty(12, order = 2)
dim(P)
#> [1] 12 12
Matrix::rankMatrix(P)
#> [1] 10
#> attr(,"method")
#> [1] "tolNorm2"
#> attr(,"useGrad")
#> [1] FALSE
#> attr(,"tol")
#> [1] 2.664535e-15

# First-order penalty (random walk)
P1 <- make_penalty(12, order = 1)
```
