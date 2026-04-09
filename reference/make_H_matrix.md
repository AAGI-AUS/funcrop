# Blend Genomic and Pedigree Relationship Matrices (H Matrix)

Creates the H-matrix for single-step GBLUP by blending a genomic
relationship matrix (G) with a pedigree-based numerator relationship
matrix (A), following the Legarra et al. (2009) formulation.

## Usage

``` r
make_H_matrix(G, A, genotyped, tau = 1, omega = 1)
```

## Arguments

- G:

  Genomic relationship matrix (square, symmetric) for genotyped
  individuals. Must have row/colnames identifying genotyped varieties.

- A:

  Pedigree relationship matrix (square, symmetric) for all individuals
  (genotyped and non-genotyped). Must have row/colnames identifying all
  varieties.

- genotyped:

  Character vector of genotyped individual IDs. Must be a subset of
  `rownames(A)` and must match `rownames(G)`.

- tau:

  Numeric: weight for \\A\_{22}^{-1}\\ in blending. Default 1.

- omega:

  Numeric: weight for \\G^{-1}\\ in blending. Default 1.

## Value

The H matrix of dimension n_all x n_all, with row/colnames from A.

## Details

The inverse of H is constructed as: \$\$H^{-1} = A^{-1} +
\begin{pmatrix} 0 & 0 \\ 0 & \omega G^{-1} - \tau A\_{22}^{-1}
\end{pmatrix}\$\$

where \\A\_{22}\\ is the submatrix of A corresponding to genotyped
individuals. H itself is then obtained by inverting \\H^{-1}\\.

Setting \\\tau = \omega = 1\\ gives the standard Legarra formulation.
Values different from 1 can be used for regularisation or when G and A
are on different scales.

## References

Legarra, A., Aguilar, I. and Misztal, I. (2009). A relationship matrix
including full pedigree and genomic information. *Journal of Dairy
Science*, 92(9), 4656–4663.
[doi:10.3168/jds.2009-2061](https://doi.org/10.3168/jds.2009-2061)

Aguilar, I., Misztal, I., Johnson, D.L., Legarra, A., Tsuruta, S. and
Lawlor, T.J. (2010). Hot topic: A unified approach to utilize
phenotypic, full pedigree, and genomic information for genetic
evaluation of Holstein final score. *Journal of Dairy Science*, 93(2),
743–752.

## See also

[`make_genomic_matrix()`](https://AAGI-AUS.github.io/funcrop/reference/make_genomic_matrix.md),
[`make_pedigree_matrix()`](https://AAGI-AUS.github.io/funcrop/reference/make_pedigree_matrix.md),
[`check_relationship_matrix()`](https://AAGI-AUS.github.io/funcrop/reference/check_relationship_matrix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Suppose A is a 50x50 pedigree matrix and G is 20x20 genomic matrix
# for the first 20 individuals
H <- make_H_matrix(G, A, genotyped = rownames(G))
dim(H)  # 50 x 50
} # }
```
