# Construct a Pedigree-Based Relationship Matrix (A Matrix)

Computes the numerator relationship matrix (A) from a pedigree
data.frame using Henderson's tabular method. The A matrix encodes the
expected proportion of alleles shared identical by descent between all
pairs of individuals.

## Usage

``` r
make_pedigree_matrix(pedigree)
```

## Arguments

- pedigree:

  A data.frame (or data.table) with exactly three columns in order:
  `id`, `sire`, `dam`. Unknown parents should be coded as `NA`, `"0"`,
  or `0`. The pedigree must be sorted such that parents appear before
  their offspring (chronological order).

## Value

A symmetric positive definite matrix of dimension n x n, where n is the
number of individuals. Row and column names are set to the individual
IDs.

## Details

Implements Henderson's (1976) tabular method without external
dependencies. Time complexity is O(n^2) in the number of individuals,
which is adequate for typical crop trial pedigrees (up to ~5000
individuals). For larger pedigrees, consider using the `pedigreemm` or
`nadiv` packages.

Inbreeding coefficients are computed on the diagonal: \\A\_{ii} = 1 +
F_i\\, where \\F_i\\ is the inbreeding coefficient of individual \\i\\.

## References

Henderson, C.R. (1976). A simple method for computing the inverse of a
numerator relationship matrix used in prediction of breeding values.
*Biometrics*, 32(1), 69–83.
[doi:10.2307/2529339](https://doi.org/10.2307/2529339)

## See also

[`make_genomic_matrix()`](https://biometryhub.github.io/funcrop/reference/make_genomic_matrix.md),
[`make_H_matrix()`](https://biometryhub.github.io/funcrop/reference/make_H_matrix.md),
[`check_relationship_matrix()`](https://biometryhub.github.io/funcrop/reference/check_relationship_matrix.md)

## Examples

``` r
ped <- data.frame(
  id   = c("A", "B", "C", "D", "E"),
  sire = c(NA,  NA,  "A", "A", "C"),
  dam  = c(NA,  NA,  "B", "B", "D")
)
A <- make_pedigree_matrix(ped)
print(round(A, 3))
#>     A   B    C    D    E
#> A 1.0 0.0 0.50 0.50 0.50
#> B 0.0 1.0 0.50 0.50 0.50
#> C 0.5 0.5 1.00 0.50 0.75
#> D 0.5 0.5 0.50 1.00 0.75
#> E 0.5 0.5 0.75 0.75 1.25
```
