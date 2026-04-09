# funcrop: Functional Data Analysis for Crop Variety Trials

Functional data analysis methods for crop variety trials using B-spline
basis functions integrated with linear mixed models. Supports
scalar-on-function regression relating primary traits (e.g., yield) to
secondary functional traits (e.g., grain-fill rate, stay-green) measured
over time. Handles two-dimensional functional traits (time x depth, time
x wavelength) via tensor product B-splines, and multi-environment trial
(MET) analyses with genotype-by-environment (GxE) interaction. Provides
a dual-backend architecture: ASReml-R (REML estimation) and bayesreml
(Bayesian MCMC via greta), with pedigree and genomic relationship matrix
support. Includes publication-ready visualisation tools for functional
profiles, coefficient functions, and GxE patterns.

## See also

Useful links:

- <https://github.com/AAGI-AUS/funcrop>

- Report bugs at <https://github.com/AAGI-AUS/funcrop/issues>

## Author

**Maintainer**: Max Moldovan <max.moldovan@adelaide.edu.au>

Authors:

- Joanne De Faveri

- Ari Verbyla

Other contributors:

- AAGI \[funder\]
