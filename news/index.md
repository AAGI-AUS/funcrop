# Changelog

## funcrop 0.2.0

### Breaking changes

- [`fit_functional_profiles()`](https://biometryhub.github.io/funcrop/reference/fit_functional_profiles.md):
  Model now uses the correct P-spline mixed model reparameterisation
  (null-space fixed effects + range-space random effects). Previously
  used raw B-spline columns without penalisation. Fitted curves and
  BLUPs will differ from v0.1.0 — the new results are statistically
  correct.
- [`fit_fda_joint()`](https://biometryhub.github.io/funcrop/reference/fit_fda_joint.md):
  Random formula now includes shared spline coefficients linking
  secondary and primary traits. Previously used only trait-by-variety
  random intercepts, defeating the functional link.
- Engine auto-detection priority changed: mgcv \> lme4 \> asreml \>
  bayesreml (open-source first). Set `funcrop_default_engine("asreml")`
  to restore old behaviour.

### New features

- **lme4 backend**: `engine = "lme4"` for `lmer()`/`glmer()`-based
  fitting. Enables GLMM extensions for non-Gaussian functional traits.
- **mgcv backend**: `engine = "mgcv"` for `gam()`/`bam()`/`gamm()`-based
  fitting. Native tensor products, automatic smoothing selection, and
  `bam()` for large datasets.
- **Engine registry**: Extensible registration system for custom
  backends via `.register_engine()`.
- [`funcrop_engines()`](https://biometryhub.github.io/funcrop/reference/funcrop_engines.md)
  now lists all 4 available backends.

### Bug fixes

- Fixed: Fitted curves now include the fixed-effect population mean
  curve (was random-effect deviations only, centred at zero).
- Fixed: `bayesreml` backend now filters invalid arguments before
  dispatch, preventing `unused argument (nsamples = 2000)` errors.
- Fixed: MET Stage 1 subscript-out-of-bounds when environments have
  different variety subsets — now uses defensive indexing.
- Fixed: `.find_spline_terms()` now recognises `Zrange_` and `Zspline_`
  column patterns from the decomposed basis.
- Improved: bayesreml default MCMC settings increased from 2000/1000 to
  4000/2000 samples/warmup for better convergence on FDA models.
- Added: Convergence quality gate between MET Stage 1 and Stage 2 —
  warns about unconverged environments before passing BLUPs.

### Documentation

- Added `inst/REVISION_LOG.md` mapping every feedback item to source
  lines and fixes.

------------------------------------------------------------------------

## funcrop 0.1.0

### New features

- Initial release.
- B-spline basis construction with penalty matrices
  ([`bspline_basis()`](https://biometryhub.github.io/funcrop/reference/bspline_basis.md)).
- Tensor product B-splines for 2D functional data
  ([`tensor_bspline_basis()`](https://biometryhub.github.io/funcrop/reference/tensor_bspline_basis.md)).
- Mixed model reparameterisation
  ([`make_Zspline()`](https://biometryhub.github.io/funcrop/reference/make_Zspline.md)).
- Dual-backend architecture: ASReml-R v4.2 (REML) and bayesreml
  (Bayesian MCMC).
- S3 classes: `fda_basis`, `fda_data`, `fda_model`.
