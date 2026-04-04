# funcrop — Project Log

## Compressed Project Context

- **Objective**: R package for Functional Data Analysis in crop variety trials
- **Core concept**: B-spline basis functions integrated with linear mixed models; scalar-on-function regression relating primary traits (yield) to secondary functional traits (grain-fill, stay-green)
- **Technical approach**: Dual-backend (ASReml-R v4.2 REML + bayesreml Bayesian MCMC); P-splines with difference penalties; tensor products for 2D; FA models for GxE
- **Key assumptions**: B-splines sufficient (no need for wavelets/Fourier); ASReml v4.2 API stable; bayesreml v0.1.0 structures adequate for FDA
- **Major adaptations**: Extending bayesreml's formula parser to handle funcrop's specialised B-spline structures
- **Status**: Phases 0–7 implemented (~12,000 lines). All R files parse. 141+ unit/integration tests written. Awaiting R CMD check and real-data validation.
- **Open issues**: Package name `funcrop` is provisional; need real datasets from CSIRO/QAAFI; NAMESPACE needs roxygen2 generation; vignettes not yet written

---

## Log Entries

### 2026-04-04 — v0.0.1 — Project Initialisation

**Summary of Changes**:
- Created full R package scaffold: DESCRIPTION, LICENSE, CI workflows, pkgdown config
- Reviewed AAGI-AU-RD-FDA proposal documents and bayesreml project
- Produced comprehensive implementation plan (IMPLEMENTATION_PLAN.md v2)
- Established dual-backend architecture (ASReml-R v4.2 + bayesreml)

**Key Decisions & Rationale**:
- Package name: `funcrop` (concise, descriptive)
- GPL-3 licence (ASReml compatibility, bayesreml alignment)
- S3 classes (simpler than S4/R6, fits R ecosystem)
- Both ASReml and bayesreml in Suggests (not Imports) — at least one required for fitting
- Auto-detect engine: ASReml first, bayesreml fallback

**Technical Updates**:
- Package scaffold: DESCRIPTION, NAMESPACE imports, .gitignore, .Rbuildignore
- CI: GitHub Actions for R-CMD-check (4 platforms), pkgdown, test-coverage
- Repository target: github.com/max578/funcrop

**Next Steps**:
1. Implement B-spline engine (basis.R, penalty.R, tensor.R)
2. Implement S3 classes (classes.R)
3. Build backend abstraction layer (backend*.R)
4. Generate simulated test data
5. Write unit tests for basis functions

### 2026-04-04 — v0.1.0 — Full Implementation Sprint

**Summary of Changes**:
- Implemented Phases 2–7 in a single session via parallel agent architecture
- Total: ~12,000 lines of R code across 16 source files, 3 test files, 1 data-raw script

**Key Decisions & Rationale**:
- `row_kronecker` consolidated into utils.R (`.row_kronecker`) — single source of truth
- Two-stage MET-FDA uses per-environment `fit_functional_profiles()` then FA on coefficients
- Single-stage MET-FDA uses stacked data with `at(environment)` structures
- Genomic matrices: VanRaden Method 1/2, Henderson tabular A, Legarra H-matrix
- All fitting functions accept `engine = c("auto", "asreml", "bayesreml")`
- Visualisation: ggplot2 primary with base R fallback

**Technical Updates**:
- Phase 2: B-spline engine (basis.R, penalty.R, tensor.R) — basis construction, P-spline penalties, mixed model reparameterisation, tensor products
- Phase 3: Backend layer (backend.R, backend_asreml.R, backend_bayesreml.R, backend_common.R) — engine dispatch, ASReml v4.2 str() structures, bayesreml formula translation, result standardisation
- Phase 4: Fitting (fit_profiles.R, scalar_on_function.R, fit_joint.R) — two-stage, single-stage, 2D FDA
- Phase 5: MET-FDA (fit_met.R) — FA/US/diag GxE on B-spline coefficients, predict_new_env()
- Phase 6: Genomics (genomic.R) — GRM, A-matrix, H-matrix, matrix repair
- Phase 7: Visualisation (plot.R) — 5 exported plot functions, backend comparison
- Testing: 141+ tests in 3 files (test-basis.R, test-classes.R, test-integration.R)
- Data: sim_grain_fill (480 rows), sim_met_fda (2160 rows) — simulated with known parameters
- Fix: .datatable.aware + missing data.table imports in funcrop-package.R

**Next Steps**:
1. Run `roxygen2::roxygenise()` to generate NAMESPACE and man/ pages
2. Run `R CMD check --as-cran` and fix any issues
3. Generate simulated datasets via `data-raw/simulate_data.R`
4. Write vignettes (introduction, grain-fill case study, bayesian-fda)
5. Create GitHub repo (max578/funcrop) and push
6. Set up pkgdown site
