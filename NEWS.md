# SurvFitSummary 1.0.0

## Initial Release (2025-11-30)

This is the first official release of SurvFitSummary, a comprehensive R package for parametric survival model fitting and diagnostic visualization.

### Major Features

#### Data Processing
* `processing_dataset()`: Flexible data processing function supporting multiple formats (CSV, Excel, SAS)
  - Automatic column mapping for ARM, SURVTIME, EVENT/CNSR
  - Support for both event (EVENT=1) and censoring (CNSR=1) indicators
  - Input validation and error handling

#### Dummy Data Generation
* `generate_dummy_survival_data()`: Generate simulated survival data for testing
  - Supports multiple treatment arms
  - Exponential event time generation with specified hazard rates
  - Administrative censoring and dropout handling
  - Reproducible results with seed parameter

#### Model Fitting Functions
* `fitting_surv_mod()`: Fit single parametric distribution
  - Seven parametric distributions: exponential, Weibull, lognormal, log-logistic, Gompertz, generalized gamma, gamma
  - Dependent (covariate) and independent (stratified) model options
  - Parameter estimates with confidence intervals
  - Variance-covariance and Cholesky decomposition matrices
  - Goodness-of-fit statistics (AIC, BIC)
  - Landmark survival probabilities
  - Survival metrics (median, mean, RMST)
  - Time scale conversion support

* `fitting_surv_mod_multiple()`: Fit and compare multiple distributions
  - Batch fitting of multiple distributions
  - Support for both dependent and independent models simultaneously
  - Standardized time scale handling (day, week, month, year)
  - Cox proportional hazards model results for comparison
  - Comprehensive output combining all model results

#### Kaplan-Meier Visualization
* `plot_km_and_parametric()`: KM curves with single parametric overlay
  - Kaplan-Meier step curves
  - Parametric distribution smooth curves
  - Optional 95% confidence intervals
  - Number at risk table
  - Automatic axis scaling

* `plot_km_and_parametric_multiple()`: KM curves with multiple parametric overlays
  - Faceted plots by distribution
  - Consistent color scheme across panels
  - Shared axes for easy comparison
  - Support for 2+ treatment arms

#### Hazard Function Visualization
* `plot_hazard_and_parametric()`: Hazard plots with single distribution
  - Non-parametric hazard estimation (bshazard)
  - Parametric hazard curves with 95% CI
  - Support for multiple treatment arms

* `plot_hazard_and_parametric_multiple()`: Hazard plots with multiple distributions
  - Faceted comparison of hazard functions
  - Consistent visualization across distributions
  - Easy identification of hazard patterns

#### Diagnostic Plots
* `plot_schoenfeld()`: Scaled Schoenfeld residuals
  - Test proportional hazards assumption
  - LOESS smoothing curves
  - Beta coefficient and p-value annotations
  - Support for multiple treatment comparisons
  - Faceted display for 3+ arms

* `plot_test_ph_assumption()`: Complementary log-log and related plots
  - Four diagnostic plots: log(-log(S)) vs log(t), log(-log(S)) vs t, Φ⁻¹(S) vs log(t), Φ⁻¹(S) vs t
  - Parallel curves indicate proportional hazards
  - Faceted layout for easy comparison

* `plot_qq_survival()`: Q-Q plots for survival distributions
  - Compare survival time distributions between arms
  - Identity line and regression fit with 95% CI
  - Faceted plots for 3+ arms vs control
  - Quantile-based comparison

### Technical Features
* Full roxygen2 documentation for all functions
* Comprehensive examples in function documentation
* Input validation and informative error messages
* Consistent API design across functions
* Support for factor and character ARM variables
* Automatic handling of reference group specification
* Non-ASCII character compliance for CRAN compatibility
* No global variable warnings (proper NSE handling)

### Dependencies
* Core: R (>= 4.0.0)
* Survival Analysis: survival, flexsurv, bshazard
* Data Manipulation: dplyr, tidyr, tibble, purrr
* Visualization: ggplot2, patchwork
* Data Import: readr, readxl, haven
* Utilities: broom, rlang, magrittr, stats

### Supported Distributions
1. Exponential (`"exp"`)
2. Weibull (`"weibull"`)
3. Lognormal (`"lnorm"`)
4. Log-logistic (`"llogis"`)
5. Gompertz (`"gompertz"`)
6. Generalized Gamma (`"gengamma"`)
7. Gamma (`"gamma"`)

### Output Components
All model fitting functions return comprehensive results including:
* Parameter estimates with standard errors and confidence intervals
* Variance-covariance matrices (full and Cholesky decomposition)
* Goodness-of-fit statistics for model comparison
* Landmark survival probabilities at specified time points
* Survival metrics (median, mean, RMST) with time scale conversion
* Cox PH model results (in `fitting_surv_mod_multiple()`)

### Quality Assurance
* Passes `R CMD check` with 0 errors, 0 warnings, 0 notes
* All functions fully documented with roxygen2
* Consistent coding style and naming conventions
* Proper handling of edge cases (single arm, missing data, etc.)
* Informative messages and warnings for users

### Known Limitations
* Generalized gamma distribution may require larger sample sizes for stable parameter estimation
* Time scale conversion assumes linear relationship (appropriate for survival time units)
* Non-parametric hazard estimation (bshazard) may be sensitive to bandwidth selection

### Future Development
Planned features for future releases:
* Additional parametric distributions (e.g., piecewise exponential, spline-based)
* Model averaging and ensemble methods
* Extrapolation diagnostics
* Interactive Shiny application
* Bayesian parameter estimation options
* Additional goodness-of-fit tests

---

For bug reports and feature requests, please visit:
https://github.com/gosukehommaEX/SurvFitSummary/issues
