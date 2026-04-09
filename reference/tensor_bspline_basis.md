# Construct a Tensor Product B-Spline Basis

Builds a two-dimensional tensor product basis from two marginal B-spline
bases, suitable for modelling functional data observed over two indices
(e.g., time x depth, time x wavelength). The penalty is the standard
additive Kronecker sum: \\P = P_1 \otimes I_2 + I_1 \otimes P_2\\.

## Usage

``` r
tensor_bspline_basis(
  x1,
  x2,
  n_knots1 = 10L,
  n_knots2 = 10L,
  degree1 = 3L,
  degree2 = 3L,
  knot_type1 = "equally_spaced",
  knot_type2 = "equally_spaced",
  penalty_order1 = 2L,
  penalty_order2 = 2L
)
```

## Arguments

- x1:

  Numeric vector of evaluation points for the first dimension.

- x2:

  Numeric vector of evaluation points for the second dimension. Must
  have the same length as `x1` (i.e., each observation has coordinates
  `(x1[i], x2[i])`).

- n_knots1, n_knots2:

  Integer; number of internal knots for each marginal basis. Defaults
  are 10.

- degree1, degree2:

  Integer; polynomial degree for each marginal basis. Defaults are 3
  (cubic).

- knot_type1, knot_type2:

  Character; knot placement for each marginal basis. One of
  `"equally_spaced"` (default), `"quantile"`, or `"custom"`.

- penalty_order1, penalty_order2:

  Integer; order of the difference penalty for each marginal basis.
  Defaults are 2.

## Value

An S3 object of class `c("fda_tensor_basis", "fda_basis")` containing:

- B:

  Tensor product basis matrix (n x p1p2), constructed as the row-wise
  Kronecker product of the two marginal basis matrices.

- basis1:

  The marginal `fda_basis` object for the first dimension.

- basis2:

  The marginal `fda_basis` object for the second dimension.

- P:

  Combined additive penalty matrix (p1p2 x p1p2), sparse. Computed as
  \\P_1 \otimes I_2 + I_1 \otimes P_2\\.

- P1:

  Marginal penalty in tensor product form: \\P_1 \otimes I_2\\ (sparse).

- P2:

  Marginal penalty in tensor product form: \\I_1 \otimes P_2\\ (sparse).

- n_basis:

  Integer; total number of tensor product basis functions (p1 times p2).

## Details

The tensor product basis matrix at observation i is:
\$\$B\_{\text{tensor}}\[i,\] = B_1\[i,\] \otimes B_2\[i,\]\$\$ where
\\\otimes\\ denotes the Kronecker product of two row vectors.

The additive penalty: \$\$P = P_1 \otimes I_2 + I_1 \otimes P_2\$\$
penalises roughness in each dimension separately, which is the standard
approach for tensor product P-splines (Marx & Eilers, 2005; Wood, 2006).

Sparse matrix representations (Matrix package) are used throughout as
tensor product bases can be large (p1 \* p2 basis functions).

## References

Marx, B.D. and Eilers, P.H.C. (2005). Multidimensional penalized signal
regression. *Technometrics*, 47(1), 13–22.

Wood, S.N. (2006). Low-rank scale-invariant tensor product smooths for
generalized additive mixed models. *Biometrics*, 62(4), 1025–1036.

## See also

[`bspline_basis()`](https://AAGI-AUS.github.io/funcrop/reference/bspline_basis.md)
for marginal basis construction,
[`make_penalty()`](https://AAGI-AUS.github.io/funcrop/reference/make_penalty.md)
for penalty matrices,
[`make_Zspline()`](https://AAGI-AUS.github.io/funcrop/reference/make_Zspline.md)
for mixed model reparameterisation.

## Examples

``` r
# 2D basis: time x depth with 100 observations
set.seed(123)
n <- 100
time_pts <- runif(n, 0, 10)
depth_pts <- runif(n, 0, 1)
tb <- tensor_bspline_basis(time_pts, depth_pts,
                           n_knots1 = 6, n_knots2 = 5)
print(tb)
#> Tensor product B-spline basis (fda_tensor_basis)
#>   Dimension 1:
#>     Degree:          3 
#>     Internal knots: 6 (equally_spaced)
#>     Basis functions: 10 
#>     Boundary:       [-0.09312, 10.04]
#>     Penalty order:   2 
#>   Dimension 2:
#>     Degree:          3 
#>     Internal knots: 5 (equally_spaced)
#>     Basis functions: 9 
#>     Boundary:       [0.0007154, 0.9954]
#>     Penalty order:   2 
#>   Total basis functions:90 (10 x 9)
#>   Observations:         100 
#>   Basis matrix dim:    100 x 90
dim(tb$B)
#> [1] 100  90

# \donttest{
# Larger tensor product basis (may be slow)
tb_large <- tensor_bspline_basis(time_pts, depth_pts,
                                 n_knots1 = 15, n_knots2 = 10)
dim(tb_large$B)
#> [1] 100 266
# }
```
