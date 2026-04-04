# funcrop — Project Log

## Compressed Project Context

- **Objective**: R package for Functional Data Analysis in crop variety trials
- **Core concept**: B-spline basis functions integrated with linear mixed models; scalar-on-function regression relating primary traits (yield) to secondary functional traits (grain-fill, stay-green)
- **Technical approach**: Dual-backend (ASReml-R v4.2 REML + bayesreml Bayesian MCMC); P-splines with difference penalties; tensor products for 2D; FA models for GxE
- **Key assumptions**: B-splines sufficient (no need for wavelets/Fourier); ASReml v4.2 API stable; bayesreml v0.1.0 structures adequate for FDA
- **Major adaptations**: Extending bayesreml's formula parser to handle funcrop's specialised B-spline structures
- **Status**: Phase 0 (scaffolding) complete; Phases 2–3 in progress
- **Open issues**: Package name `funcrop` is provisional; need real datasets from CSIRO/QAAFI; bayesreml tensor product support TBD

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
