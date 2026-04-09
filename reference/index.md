# Package index

## B-Spline Basis Construction

Functions for constructing B-spline bases, penalties, and design
matrices

- [`bspline_basis()`](https://AAGI-AUS.github.io/funcrop/reference/bspline_basis.md)
  : Construct a B-Spline Basis with Difference Penalty
- [`tensor_bspline_basis()`](https://AAGI-AUS.github.io/funcrop/reference/tensor_bspline_basis.md)
  : Construct a Tensor Product B-Spline Basis
- [`make_Zspline()`](https://AAGI-AUS.github.io/funcrop/reference/make_Zspline.md)
  : Reparameterise B-Spline Basis for Mixed Model Representation
- [`make_penalty()`](https://AAGI-AUS.github.io/funcrop/reference/make_penalty.md)
  : Construct a Difference Penalty Matrix

## Data Structures & S3 Methods

S3 classes and methods for functional data objects

- [`fda_data()`](https://AAGI-AUS.github.io/funcrop/reference/fda_data.md)
  : Create a functional data object
- [`print(`*`<fda_basis>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/print.fda_basis.md)
  : Print Method for fda_basis Objects
- [`print(`*`<fda_comparison>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/print.fda_comparison.md)
  : Print method for fda_comparison
- [`print(`*`<fda_data>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/print.fda_data.md)
  : Print method for fda_data
- [`print(`*`<fda_model>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/print.fda_model.md)
  : Print method for fda_model
- [`print(`*`<fda_tensor_basis>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/print.fda_tensor_basis.md)
  : Print Method for fda_tensor_basis Objects
- [`summary(`*`<fda_comparison>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/summary.fda_comparison.md)
  : Summary method for fda_comparison
- [`summary(`*`<fda_data>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/summary.fda_data.md)
  : Summary method for fda_data
- [`summary(`*`<fda_model>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/summary.fda_model.md)
  : Summary method for fda_model
- [`plot(`*`<fda_basis>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/plot.fda_basis.md)
  : Plot Method for fda_basis Objects
- [`plot(`*`<fda_comparison>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/plot.fda_comparison.md)
  : Plot method for fda_comparison
- [`plot(`*`<fda_data>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/plot.fda_data.md)
  : Plot method for fda_data (spaghetti plot)
- [`plot(`*`<fda_model>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/plot.fda_model.md)
  : Plot method for fda_model
- [`coef(`*`<fda_model>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/coef.fda_model.md)
  : Extract coefficient function beta(t) from fda_model
- [`fitted(`*`<fda_model>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/fitted.fda_model.md)
  : Extract fitted curves from fda_model
- [`residuals(`*`<fda_model>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/residuals.fda_model.md)
  : Extract residuals from fda_model
- [`vcov(`*`<fda_model>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/vcov.fda_model.md)
  : Extract variance-covariance information from fda_model
- [`as.data.table.fda_data()`](https://AAGI-AUS.github.io/funcrop/reference/as.data.table.fda_data.md)
  : Coerce fda_data to data.table
- [`` `[`( ``*`<fda_data>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/sub-.fda_data.md)
  : Subset method for fda_data

## Model Fitting

Core model fitting functions for FDA in crop trials

- [`fit_functional_profiles()`](https://AAGI-AUS.github.io/funcrop/reference/fit_functional_profiles.md)
  : Fit Variety-Specific Functional Profiles (Stage 1)
- [`scalar_on_function()`](https://AAGI-AUS.github.io/funcrop/reference/scalar_on_function.md)
  : Scalar-on-Function Regression (Stage 2)
- [`fit_fda_joint()`](https://AAGI-AUS.github.io/funcrop/reference/fit_fda_joint.md)
  : Fit a Joint Model for Primary and Functional Secondary Traits
- [`fit_2d_functional()`](https://AAGI-AUS.github.io/funcrop/reference/fit_2d_functional.md)
  : Fit 2D Functional Profiles via Tensor Product B-Splines
- [`scalar_on_2d_function()`](https://AAGI-AUS.github.io/funcrop/reference/scalar_on_2d_function.md)
  : Scalar-on-2D-Function Regression
- [`fit_fda_met()`](https://AAGI-AUS.github.io/funcrop/reference/fit_fda_met.md)
  : Fit Functional Data Analysis Model for Multi-Environment Trials

## Prediction & Comparison

Prediction, comparison, and inference tools

- [`predict(`*`<fda_model>`*`)`](https://AAGI-AUS.github.io/funcrop/reference/predict.fda_model.md)
  : Predict method for fda_model
- [`predict_new_env()`](https://AAGI-AUS.github.io/funcrop/reference/predict_new_env.md)
  : Predict Variety Performance in New Environments
- [`compare_methods()`](https://AAGI-AUS.github.io/funcrop/reference/compare_methods.md)
  : Compare Multiple FDA Model Fits

## Visualisation

Publication-ready plots for FDA results

- [`plot_functional_profiles()`](https://AAGI-AUS.github.io/funcrop/reference/plot_functional_profiles.md)
  : Plot Fitted Functional Profiles
- [`plot_coefficient_function()`](https://AAGI-AUS.github.io/funcrop/reference/plot_coefficient_function.md)
  : Plot the Coefficient Function
- [`plot_gxe_heatmap()`](https://AAGI-AUS.github.io/funcrop/reference/plot_gxe_heatmap.md)
  : Heatmap of Genotype-by-Environment Interaction
- [`plot_fa_biplot()`](https://AAGI-AUS.github.io/funcrop/reference/plot_fa_biplot.md)
  : Factor-Analytic Biplot for GxE Interaction
- [`plot_backend_comparison()`](https://AAGI-AUS.github.io/funcrop/reference/plot_backend_comparison.md)
  : Compare REML and Bayesian Model Fits

## Backend Management

Engine selection and backend utilities

- [`funcrop_engines()`](https://AAGI-AUS.github.io/funcrop/reference/funcrop_engines.md)
  : List Available Estimation Engines
- [`funcrop_default_engine()`](https://AAGI-AUS.github.io/funcrop/reference/funcrop_default_engine.md)
  : Get or Set the Default Estimation Engine

## Genomic Integration

Relationship matrix construction

- [`make_genomic_matrix()`](https://AAGI-AUS.github.io/funcrop/reference/make_genomic_matrix.md)
  : Construct a Genomic Relationship Matrix
- [`make_pedigree_matrix()`](https://AAGI-AUS.github.io/funcrop/reference/make_pedigree_matrix.md)
  : Construct a Pedigree-Based Relationship Matrix (A Matrix)
- [`make_H_matrix()`](https://AAGI-AUS.github.io/funcrop/reference/make_H_matrix.md)
  : Blend Genomic and Pedigree Relationship Matrices (H Matrix)
- [`check_relationship_matrix()`](https://AAGI-AUS.github.io/funcrop/reference/check_relationship_matrix.md)
  : Check and Repair a Relationship Matrix

## Datasets

Simulated datasets for tutorials and testing

- [`sim_grain_fill`](https://AAGI-AUS.github.io/funcrop/reference/sim_grain_fill.md)
  : Simulated Grain-Fill Trial Data
- [`sim_met_fda`](https://AAGI-AUS.github.io/funcrop/reference/sim_met_fda.md)
  : Simulated Multi-Environment Trial with Functional NDVI Trait
