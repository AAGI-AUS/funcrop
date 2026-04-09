# Revision Log — funcrop v0.2.0

**Purpose**: Tracks how each feedback item was addressed, what advanced practice
was applied, and provides code-to-specification traceability.

**Revision started**: 2026-04-09
**Feedback sources**: Jules (reviewer), R_Console_review.txt, R_Console_review2.txt,
4 reverse-engineered mathematical framework PDFs.

---

## Phase 0: Audit & Bug-to-Source Mapping

### Feedback Item F01: Spline Z columns not in random formula (fit_functional_profiles)

- **Source**: `feedback/functional_model_note.pdf`
- **Issue**: `make_Zspline()` eigendecomposition computed (line 279) producing
  null-space X and range-space Z, but the actual model formula uses raw B-spline
  columns `Bsp_1, ..., Bsp_K` as both fixed AND random effects (lines 384-404).
  The decomposed X/Z matrices are used only to build a Kronecker Z_full (line 344)
  that is never referenced in the model formula.
- **Consequence**: No P-spline penalisation. Population mean curve is overfit
  (K fixed effects instead of d null-space terms). Random effects are IID on
  raw basis coefficients, not penalised through the smoothing penalty.
- **Before**: `fixed = value ~ 1 + block + Bsp_1 + ... + Bsp_K`,
  `random = ~ variety_f:Bsp_1 + ... + variety_f:Bsp_K`
- **Required fix**: `fixed = value ~ 1 + block + Xnull_1 + ... + Xnull_d`,
  `random = ~ variety_f:Zrange_1 + ... + variety_f:Zrange_q`
  where d = dim(null(P)) and q = rank(P).
- **Files**: `R/fit_profiles.R:278-448`
- **Phase**: 2.1 (Critical)

### Feedback Item F02: Fitted curves = random BLUPs only, missing fixed mean

- **Source**: `feedback/functional_model_note.pdf`
- **Issue**: `.extract_fitted_curves()` (backend_common.R:126-281) reconstructs
  curves as `fitted_v(t) = B_grid %*% alpha_v` using only random-effect BLUPs.
  The fixed-effect null-space component mu(t) = X(t)'beta is never added back.
- **Consequence**: Fitted curves are centred around zero (variety deviations)
  instead of absolute predicted profiles.
- **Before**: `fitted = B_grid %*% alpha_v` (random only)
- **Required fix**: `fitted = X_grid %*% beta_hat + Z_grid %*% alpha_v` (full)
- **Files**: `R/backend_common.R:126-281`
- **Phase**: 1.4 (Critical)

### Feedback Item F03: Tensor basis not in formula (fit_2D_functional)

- **Source**: `feedback/fit_2d_functional_framework.pdf`
- **Issue**: `fit_2D_functional()` constructs tensor product basis B_tensor and
  penalty matrices P, P_1, P_2, but the model formula passed to the engine is
  approximately `y ~ 1, random = ~ group` with optional AR1xAR1 residuals.
  The tensor basis columns are stored in extras but not in the formula.
- **Consequence**: Model fits a random-intercept model instead of a 2D functional
  mixed model. "The code builds a full 2D spline cathedral, but the formula that
  walks into the engine looks more like a modest random-intercept hut."
- **Files**: `R/fit_joint.R` (scalar_on_2d_function section)
- **Phase**: 2.2 (Critical)

### Feedback Item F04: Tensor basis — same as F03 (duplicate)

Consolidated into F03.

### Feedback Item F05: Integrated basis c_k not in primary-trait random term (fit_fda_joint)

- **Source**: `feedback/fit_fda_joint_framework.pdf`
- **Issue**: In `fit_fda_joint()`, the z_spl columns are computed and added to
  the stacked data (lines 388-429), with secondary rows getting Z(t) and primary
  rows getting the integrated c_k = integral(z_k(t))dt. However, the random
  formula at lines 432-446 is:
  `~ us(trait_type):variety + at(trait_type, 'secondary'):variety`
  This does NOT reference z_spl columns. The integrated basis vectors for the
  primary trait rows are computed but never used in the model formula.
- **Consequence**: The joint model degrades to a bivariate random-intercept
  model (variety main effect per trait) without the spline-mediated functional
  link between secondary and primary traits. The core innovation of the joint
  model — that the scalar primary trait depends on the integrated latent
  functional process — is missing from the actual fit.
- **Before**: `random = ~ us(trait_type):variety + at(trait_type, 'secondary'):variety`
- **Required fix**: `random = ~ variety:z_spl_1 + ... + variety:z_spl_q`
  where z_spl columns contain Z(t) for secondary rows and c_k for primary rows,
  creating the shared latent coefficient structure alpha_g that drives both traits.
- **Files**: `R/fit_joint.R:375-483`
- **Phase**: 2.3 (Critical)

### Feedback Item F06: Joint model formula — consolidated into F05

### Feedback Item F07: FA extraction returns NA (fit_fda_met)

- **Source**: `feedback/fit_fda_met_framework_fixed.pdf`
- **Issue**: `.extract_fa_results()` attempts to recover FA loadings and variety
  scores from the fitted model but falls back to NA placeholders when the VC
  table format doesn't match expected patterns.
- **Files**: `R/fit_met.R` (`.extract_fa_results()` function)
- **Phase**: 1.5 (High)

### Feedback Item F08: Two-stage MET Stage 2 on scalar BLUPs not coefficient vectors

- **Source**: `feedback/fit_fda_met_framework_fixed.pdf`
- **Issue**: Stage 2 model specification at line 567 uses `blup_value` as
  response. Looking at the data assembly (lines 500-520), Stage 1 BLUPs are
  correctly unpacked as (environment, variety, coef_idx, blup_value). The model
  then treats `coef_f` (coefficient index as factor) as a fixed effect.
  However, the GxE structure is on `variety` across `environment`, not on
  `variety:coef_f`. This means all spline coefficients share the same GxE
  covariance — the FA structure applies to variety main effects, not to the
  full coefficient vector.
- **Required fix**: Either (a) fit separate Stage 2 models per coefficient
  index with shared FA structure, or (b) stack coefficients and model with
  `coef_idx` interacting with the FA structure.
- **Files**: `R/fit_met.R:540-591`
- **Phase**: 2.5 (High)

### Feedback Item F09: ASReml v3 syntax

- **Source**: `feedback/feedback_notes.txt`
- **Issue**: Jules notes the code uses "an old version of asreml (v3) for its
  syntax." DESCRIPTION specifies asreml >= 4.2.0 and the code uses v4 syntax
  (`residual =` not `rcov =`, `vm()`, etc.). However, `.asreml_build_spline_str()`
  generates `str()` terms that may not work in v4.2's current API.
- **Resolution**: Add explicit v3 backend. Current code is v4-targeted.
- **Files**: `R/backend_asreml.R`
- **Phase**: 4 (High)

### Feedback Item F10: 60+ hidden functions hard to audit

- **Source**: `feedback/feedback_notes.txt`
- **Issue**: Code complexity. 17 source files, ~12,000 lines, 60+ internal
  functions. No mathematical specification document.
- **Resolution**: Produce `funcrop_mathematical_specification.pdf` (Phase 2.5)
  and improve internal documentation.
- **Phase**: 2.5 + 8 (Medium)

### Feedback Item F11: `object 'spline_coef' not found`

- **Source**: `R_Console_review.txt:144`
- **Issue**: `.asreml_build_spline_str()` (backend_asreml.R:409,448,457)
  generates formula strings referencing a column `spline_coef` that is never
  created in any fitting function. The fitting functions use `Bsp_k` columns
  instead. This function appears to be dead code — it builds str() formulas
  but is never called from the main fitting path. The error likely occurs if
  a code path attempts to use the str() approach.
- **Root cause**: Variable name mismatch between helper function and callers.
- **Files**: `R/backend_asreml.R:365-471`
- **Phase**: 1.1 (Critical)

### Feedback Item F12: `unused argument (nsamples = 2000)`

- **Source**: `R_Console_review.txt:152`
- **Issue**: `backend_bayesreml.R:74` passes `n_samples` (with underscore) to
  `bayesreml::bayesreml()`. The error message says `nsamples` (no underscore)
  is the unused argument. This suggests either (a) bayesreml's actual parameter
  name differs from what funcrop passes, or (b) the parameter was renamed in
  bayesreml and funcrop wasn't updated.
- **Root cause**: Parameter name mismatch with bayesreml API.
- **Files**: `R/backend_bayesreml.R:62-77`
- **Phase**: 1.2 (Critical)

### Feedback Item F13: `subscript out of bounds` in MET Stage 1

- **Source**: `R_Console_review.txt:185`
- **Issue**: In `.fit_met_two_stage()` (fit_met.R:446-486), the per-environment
  loop extracts `stage1_models[[e]]$extras$spline_blups`. If the spline_blups
  matrix has fewer rows than `n_varieties` (because some varieties are absent
  in an environment), downstream code that indexes by `variety_levels` will
  subscript out of bounds.
- **Root cause**: Missing defensive indexing for environments with incomplete
  variety sets.
- **Files**: `R/fit_met.R:446-520`
- **Phase**: 1.3 (Critical)

### Feedback Item F14: bayesreml Rhat up to 3.1, ESS as low as 11

- **Source**: `R_Console_review2.txt:276`
- **Issue**: bayesreml MCMC convergence is poor. Default settings
  (n_samples=2000, warmup=1000, chains=4) insufficient for complex FDA models.
- **Root cause**: FDA models with many spline coefficients per variety need
  more iterations and better initialisation.
- **Files**: `R/backend_bayesreml.R:62-67`
- **Phase**: 1.6 (High)

### Feedback Item F15: Rhat = Inf, ESS = 0 in MET Stage 2

- **Source**: `R_Console_review2.txt:880`
- **Issue**: Complete convergence failure in Stage 2 MET model. Likely caused
  by Stage 1 failures propagating bad BLUPs into Stage 2.
- **Root cause**: No convergence pre-check between stages; no quality filter
  on Stage 1 results before passing to Stage 2.
- **Files**: `R/fit_met.R:488-591`
- **Phase**: 1.6 (High)

### Feedback Item F16: "No spline random effect terms found"

- **Source**: `R_Console_review2.txt:279`
- **Issue**: `.find_spline_terms()` (backend_common.R:427-451) searches
  `model$G.param` (ASReml) or `model$ranef` (bayesreml) for patterns matching
  "spline|basis|bspline|Bsp|...". Term names in ASReml's G.param for
  `variety_f:Bsp_1` interactions may be stored differently (e.g., as a single
  compound term `variety_f:Bsp_1` or under a different naming convention).
- **Root cause**: Pattern matching in `.find_spline_terms()` doesn't match the
  actual term names produced by the engine.
- **Files**: `R/backend_common.R:427-451`
- **Phase**: 1.4 (High)

---

## Summary: Critical Implementation Gaps

The most significant finding from the audit is that **three of the four major
fitting functions have specification gaps** where computational infrastructure
(basis matrices, penalty decompositions, integration weights) is built but
not connected to the model formula passed to the engine:

| Function | Infrastructure Built | Formula Actually Used |
|----------|---------------------|----------------------|
| `fit_functional_profiles()` | Eigendecomposition → X_null + Z_range | Raw B-spline columns Bsp_k (no penalisation) |
| `fit_fda_joint()` | z_spl columns with integrated basis for primary | us(trait):variety (no spline link) |
| `fit_2D_functional()` | Tensor product B_tensor, P_1, P_2 | Random intercept ~ group |
| `fit_fda_met()` | Stage 2 on coefficient vectors | FA structure on variety, not on variety:coef |

These are not cosmetic issues — they mean the fitted models are mathematically
different from the documented/intended models.
