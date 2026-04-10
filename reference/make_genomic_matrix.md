# Construct a Genomic Relationship Matrix

Computes the genomic relationship matrix (G) from marker data using the
VanRaden (2008) Method 1 or Method 2. This is a convenience function;
users may also supply pre-computed matrices to fitting functions
directly.

## Usage

``` r
make_genomic_matrix(
  markers,
  method = c("vanraden1", "vanraden2"),
  min_maf = 0.01,
  scale = TRUE
)
```

## Arguments

- markers:

  Numeric matrix of dimension n_varieties x n_markers. Coded as 0, 1, 2
  (minor allele counts) or -1, 0, 1 (centred). Row names should identify
  varieties; column names (if present) identify markers.

- method:

  Character: `"vanraden1"` (default) or `"vanraden2"`.

  - **VanRaden Method 1**: \\G = ZZ' / (2 \sum p_j(1-p_j))\\, where
    \\Z\_{ij} = M\_{ij} - 2p_j\\ and \\p_j\\ is the allele frequency at
    marker \\j\\.

  - **VanRaden Method 2**: \\G = ZDZ' / n_m\\, where \\D =
    \mathrm{diag}(1/(2 p_j (1-p_j)))\\.

- min_maf:

  Numeric: minimum minor allele frequency filter. Markers with MAF below
  this threshold are removed. Default 0.01.

- scale:

  Logical: if `TRUE` (default), scale G so that the mean diagonal equals
  1.

## Value

A symmetric positive (semi-)definite matrix of dimension n_varieties x
n_varieties, with row/colnames from `rownames(markers)`.

## Details

Missing marker values (`NA`) are imputed to the column mean (i.e., twice
the allele frequency) prior to computation. If many markers are missing,
consider using a dedicated imputation tool first.

## References

VanRaden, P.M. (2008). Efficient methods to compute genomic predictions.
*Journal of Dairy Science*, 91(11), 4414–4423.
[doi:10.3168/jds.2007-0980](https://doi.org/10.3168/jds.2007-0980)

## See also

[`make_pedigree_matrix()`](https://biometryhub.github.io/funcrop/reference/make_pedigree_matrix.md),
[`make_H_matrix()`](https://biometryhub.github.io/funcrop/reference/make_H_matrix.md),
[`check_relationship_matrix()`](https://biometryhub.github.io/funcrop/reference/check_relationship_matrix.md)

## Examples

``` r
set.seed(123)
# Simulate 20 varieties, 500 markers (coded 0/1/2)
M <- matrix(sample(0:2, 20 * 500, replace = TRUE,
                    prob = c(0.25, 0.50, 0.25)),
            nrow = 20, ncol = 500)
rownames(M) <- paste0("var_", seq_len(20))

G <- make_genomic_matrix(M, method = "vanraden1")
dim(G)
#> [1] 20 20
isSymmetric(G)
#> [1] TRUE
```
