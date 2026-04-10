# funcrop — Project Log

## Compressed Project Context

- **Objective**: R package for Functional Data Analysis in crop variety
  trials
- **Core concept**: B-spline basis functions integrated with linear
  mixed models; scalar-on-function regression relating primary traits
  (yield) to secondary functional traits (grain-fill, stay-green)
- **Technical approach**: Multi-backend (ASReml-R v4.2 REML, bayesreml
  Bayesian MCMC, lme4, mgcv); P-splines with difference penalties;
  tensor products for 2D; FA models for GxE
- **Key assumptions**: B-splines sufficient (no need for
  wavelets/Fourier); ASReml v4.2 API stable; bayesreml v0.1.0 structures
  adequate for FDA
- **Major adaptations**: Extensible engine registry; correct P-spline
  mixed-model reparameterisation throughout
- **Status**: v0.2.0 — Major revision addressing reviewer feedback. 4
  backends (mgcv, lme4, ASReml, bayesreml). 218 tests passing. R CMD
  check: 0 errors, 0 warnings.
- **Open issues**: ASReml v3 backend (needs v3 installation for
  testing); mathematical specification document (in progress); real-data
  validation pending

------------------------------------------------------------------------

## Log Entries

### 2026-04-09 — v0.2.0 — Major Revision (Reviewer Feedback + Backend Expansion)

**Summary of Changes**: - Complete revision addressing Jules’ reviewer
feedback (7 feedback documents) - Fixed critical
implementation–specification gaps in 3 of 4 fitting functions - Added
lme4 and mgcv backends (total: 4 backends) - Created extensible engine
registration infrastructure - Produced comprehensive REVISION_LOG.md
mapping all feedback to source fixes

**Key Decisions & Rationale**: - Correct P-spline reparameterisation:
null-space (fixed) + range-space (random) — previous code used raw
B-spline columns without penalisation - Engine priority changed to mgcv
\> lme4 \> asreml \> bayesreml (open-source first) - Registration-based
engine architecture for future extensibility - bayesreml MCMC defaults
doubled (4000/2000 from 2000/1000) for FDA convergence

**Technical Updates**: - `fit_profiles.R`: Xnull_k (fixed) + Zrange_k
(random) replace Bsp_k columns - `fit_joint.R`: z_spl columns now in
random formula (shared latent link) - `.reconstruct_variety_curves()`:
includes fixed-effect mean curve - `backend.R`: dispatch handles 4
engines - New files: `engine_registry.R`, `engine_lme4.R`,
`engine_mgcv.R` - Tests: 218 passing (up from 214) - R CMD check: 0
errors, 0 warnings, 1 NOTE (CRAN new submission)

**Next Steps**: 1. Complete mathematical specification document (LaTeX,
in progress) 2. ASReml v3 backend (Phase 4 — needs v3 installation) 3.
Testing overhaul — argument combination tests, cross-backend validation
4. Documentation — 7 vignettes, update README

### 2026-04-04 — v0.0.1 — Project Initialisation

**Summary of Changes**: - Created full R package scaffold: DESCRIPTION,
LICENSE, CI workflows, pkgdown config - Reviewed AAGI-AU-RD-FDA proposal
documents and bayesreml project - Produced comprehensive implementation
plan (IMPLEMENTATION_PLAN.md v2) - Established dual-backend architecture
(ASReml-R v4.2 + bayesreml)

**Key Decisions & Rationale**: - Package name: `funcrop` (concise,
descriptive) - GPL-3 licence (ASReml compatibility, bayesreml
alignment) - S3 classes (simpler than S4/R6, fits R ecosystem) - Both
ASReml and bayesreml in Suggests (not Imports) — at least one required
for fitting - Auto-detect engine: ASReml first, bayesreml fallback

**Technical Updates**: - Package scaffold: DESCRIPTION, NAMESPACE
imports, .gitignore, .Rbuildignore - CI: GitHub Actions for R-CMD-check
(4 platforms), pkgdown, test-coverage - Repository target:
github.com/max578/funcrop

**Next Steps**: 1. Implement B-spline engine (basis.R, penalty.R,
tensor.R) 2. Implement S3 classes (classes.R) 3. Build backend
abstraction layer (backend\*.R) 4. Generate simulated test data 5. Write
unit tests for basis functions

### 2026-04-04 — v0.1.0 — Full Implementation Sprint

**Summary of Changes**: - Implemented Phases 2–7 in a single session via
parallel agent architecture - Total: ~12,000 lines of R code across 16
source files, 3 test files, 1 data-raw script

**Key Decisions & Rationale**: - `row_kronecker` consolidated into
utils.R (`.row_kronecker`) — single source of truth - Two-stage MET-FDA
uses per-environment
[`fit_functional_profiles()`](https://biometryhub.github.io/funcrop/reference/fit_functional_profiles.md)
then FA on coefficients - Single-stage MET-FDA uses stacked data with
`at(environment)` structures - Genomic matrices: VanRaden Method 1/2,
Henderson tabular A, Legarra H-matrix - All fitting functions accept
`engine = c("auto", "asreml", "bayesreml")` - Visualisation: ggplot2
primary with base R fallback

**Technical Updates**: - Phase 2: B-spline engine (basis.R, penalty.R,
tensor.R) — basis construction, P-spline penalties, mixed model
reparameterisation, tensor products - Phase 3: Backend layer (backend.R,
backend_asreml.R, backend_bayesreml.R, backend_common.R) — engine
dispatch, ASReml v4.2 str() structures, bayesreml formula translation,
result standardisation - Phase 4: Fitting (fit_profiles.R,
scalar_on_function.R, fit_joint.R) — two-stage, single-stage, 2D FDA -
Phase 5: MET-FDA (fit_met.R) — FA/US/diag GxE on B-spline coefficients,
predict_new_env() - Phase 6: Genomics (genomic.R) — GRM, A-matrix,
H-matrix, matrix repair - Phase 7: Visualisation (plot.R) — 5 exported
plot functions, backend comparison - Testing: 141+ tests in 3 files
(test-basis.R, test-classes.R, test-integration.R) - Data:
sim_grain_fill (480 rows), sim_met_fda (2160 rows) — simulated with
known parameters - Fix: .datatable.aware + missing data.table imports in
funcrop-package.R

**Next Steps**: 1. Run `roxygen2::roxygenise()` to generate NAMESPACE
and man/ pages 2. Run `R CMD check --as-cran` and fix any issues 3.
Generate simulated datasets via `data-raw/simulate_data.R` 4. Write
vignettes (introduction, grain-fill case study, bayesian-fda) 5. Create
GitHub repo (max578/funcrop) and push 6. Set up pkgdown site
