# Check and Repair a Relationship Matrix

Ensures a relationship matrix is symmetric, positive definite (bending
if necessary), and properly scaled. Useful as a pre-processing step
before passing relationship matrices to mixed model fitting functions.

## Usage

``` r
check_relationship_matrix(K, tol = 1e-06, bend = TRUE)
```

## Arguments

- K:

  Numeric square matrix (relationship matrix).

- tol:

  Numeric: eigenvalue tolerance. Eigenvalues below this value are raised
  to `tol` when bending. Default `1e-6`.

- bend:

  Logical: if `TRUE` (default), bend the matrix to ensure positive
  definiteness by replacing negative/near-zero eigenvalues with `tol`.
  If `FALSE`, only checks and reports issues without modifying.

## Value

If `bend = TRUE`, returns the repaired symmetric positive definite
matrix (with row/colnames preserved). If `bend = FALSE`, returns the
symmetrised matrix (or raises an error if not PD).

## Details

The bending procedure uses spectral decomposition:

1.  Decompose \\K = U \Lambda U'\\.

2.  Replace any \\\lambda_i \< \mathrm{tol}\\ with \\\mathrm{tol}\\.

3.  Reconstruct \\K^\* = U \Lambda^\* U'\\.

This is the simplest form of matrix bending (Jorjani et al., 2003). For
more sophisticated approaches (e.g., preserving diagonal elements), see
the `Matrix` or `nadiv` packages.

## References

Jorjani, H., Klei, L. and Emanuelson, U. (2003). A simple method for
weighted bending of genetic (co)variance matrices. *Journal of Dairy
Science*, 86(2), 677–679.
[doi:10.3168/jds.S0022-0302(03)73646-7](https://doi.org/10.3168/jds.S0022-0302%2803%2973646-7)

## See also

[`make_genomic_matrix()`](https://AAGI-AUS.github.io/funcrop/reference/make_genomic_matrix.md),
[`make_pedigree_matrix()`](https://AAGI-AUS.github.io/funcrop/reference/make_pedigree_matrix.md),
[`make_H_matrix()`](https://AAGI-AUS.github.io/funcrop/reference/make_H_matrix.md)

## Examples

``` r
# Create a near-singular matrix
set.seed(42)
K <- crossprod(matrix(rnorm(30), 6, 5))
K[1, 1] <- 0  # Make it singular
K_fixed <- check_relationship_matrix(K)
#> Matrix has 1 eigenvalue(s) below tolerance (min = -2.4775e+00).
#> Matrix bent to positive definiteness.
eigen(K_fixed)$values  # All positive
#> [1] 22.052715 15.165346  6.413134  1.759685  0.000001
```
