# Reparameterise B-Spline Basis for Mixed Model Representation

Decomposes a B-spline basis into fixed-effect (X) and random-effect (Z)
components, enabling the equivalence between penalised splines and
linear mixed models (Wand, 2003; Wood, 2017). This is essential for
fitting P-splines via REML in mixed model software (e.g., ASReml-R,
lme4).

## Usage

``` r
make_Zspline(basis, constraint = c("decompose", "absorb"))
```

## Arguments

- basis:

  An object of class `"fda_basis"` (as returned by
  [`bspline_basis()`](https://AAGI-AUS.github.io/funcrop/reference/bspline_basis.md)).

- constraint:

  Character; reparameterisation method. One of:

  `"decompose"`

  :   (Default) Spectral decomposition of the penalty matrix.
      Eigenvectors corresponding to zero eigenvalues form the fixed
      effects (X, the null space), whilst those with positive
      eigenvalues form the random effects (Z, the range space), scaled
      by the inverse square root of their eigenvalues.

  `"absorb"`

  :   Absorb the penalty into Z via Cholesky decomposition of the
      Moore-Penrose pseudoinverse. The resulting Z has an identity
      penalty structure, which is convenient for standard mixed model
      software.

## Value

A list with components:

- X:

  Fixed-effect design matrix (n x d), where d = `penalty_order` is the
  dimension of the null space of the penalty. Typically contains
  polynomial terms up to degree `penalty_order - 1`.

- Z:

  Random-effect design matrix (n x r), where r = `n_basis - d`. In the
  `"decompose"` method, columns are orthogonal; in the `"absorb"`
  method, the penalty is absorbed so that the implied penalty on the
  random effect coefficients is the identity.

- penalty:

  Sparse penalty matrix for the Z component. Identity matrix (dimension
  r x r) for the `"absorb"` method; identity for `"decompose"` when
  eigenvalues are appropriately absorbed into Z.

- rank:

  Integer; rank of the original penalty matrix
  (`n_basis - penalty_order`).

- null_dim:

  Integer; dimension of the null space (= `penalty_order`).

- constraint:

  Character; the method used.

## Details

The key identity exploited is: a penalised spline fit minimising
\$\$\\y - B\beta\\^2 + \lambda \beta^T P \beta\$\$ is equivalent to the
BLUP from a mixed model \$\$y = X\alpha + Z u + \varepsilon\$\$ where
\\u \sim N(0, \sigma_u^2 I)\\ and the smoothing parameter \\\lambda =
\sigma^2 / \sigma_u^2\\ is estimated via REML.

For `"decompose"`: the spectral decomposition \\P = U \Lambda U^T\\
separates the null space (zero eigenvalues, giving X) from the range
space (positive eigenvalues, giving Z scaled by \\\Lambda^{-1/2}\\).

For `"absorb"`: uses the pseudoinverse \\P^+\\ and its Cholesky factor
\\L\\ such that \\P^+ = LL^T\\. Then \\Z = B\_{\text{range}} L\\ absorbs
the penalty, yielding identity covariance.

## References

Wand, M.P. (2003). Smoothing and mixed models. *Computational
Statistics*, 18(2), 223–249.

Wood, S.N. (2017). *Generalized Additive Models: An Introduction with R*
(2nd ed.). Chapman & Hall/CRC.

## See also

[`bspline_basis()`](https://AAGI-AUS.github.io/funcrop/reference/bspline_basis.md),
[`make_penalty()`](https://AAGI-AUS.github.io/funcrop/reference/make_penalty.md).

## Examples

``` r
x <- seq(0, 1, length.out = 100)
basis <- bspline_basis(x, n_knots = 10, degree = 3)

# Decompose method (default)
mm <- make_Zspline(basis, constraint = "decompose")
dim(mm$X)
#> [1] 100   2
dim(mm$Z)
#> [1] 100  12

# Absorb method -- penalty becomes identity
mm2 <- make_Zspline(basis, constraint = "absorb")
dim(mm2$Z)
#> [1] 100  12
```
