# SurvFitSummary 2.0.0

Version 2.0.0 redesigns the package (development in progress). The functions
of version 1.0.0 are replaced step by step.

* `processing_dataset()` keeps covariate and stratification columns
  (`column_covariates`, `column_strata`), returns `ARM` as a factor with the
  control arm (`control_arm`) as the first level, checks that the times are
  positive and that the indicator contains only 0 and 1, removes rows with
  missing values with a warning, and also reads SAS transport (`.xpt`)
  files. CSV files are read with `utils::read.csv()`.
* `generate_dummy_survival_data()` restores the random number generator
  state of the session (the simulated values are unchanged).
* New `summarise_data()`: numbers of patients, events and censored
  observations, and the median with its confidence interval, by arm.
* New `km_table()`: Kaplan-Meier estimates with standard errors and
  confidence intervals, and the numbers at risk, events and censored
  observations at regular times.
* New `estimate_hazard_ratio()`: hazard ratios from unadjusted, adjusted,
  unstratified and stratified Cox models, with the control arm as the
  reference (in version 1.0.0 the reference arm of the Cox model did not
  follow `control_arm`).
* New `test_proportional_hazards()`: Grambsch-Therneau test of
  proportional hazards.
* New `estimate_empirical_hazard()`: smoothed empirical hazard by arm with
  `bshazard::bshazard()`.
* New `fit_survival_models()` replaces `fitting_surv_mod()` and
  `fitting_surv_mod_multiple()`: fits all the parametric and Royston-Parmar
  spline models once, for each arm and for both arms with the arm as a
  covariate, and returns an object of class `survfitsummary` used by all the
  tables. Twelve distributions are available (exponential, Weibull in the
  accelerated failure time and proportional hazards forms, Rayleigh,
  lognormal, log-logistic, Gompertz, gamma, generalized gamma and
  generalized F, each also in its original parameterisation); the default
  is the seven standard models of NICE DSU TSD 14. Spline models on the
  hazard, odds and normal scales with 1 to 3 knots at quantiles of the log
  event times or at times given by the user. A model that fails or warns is
  recorded with its message and the other models are still fitted; the
  optimiser uses a smaller tolerance (`reltol = 1e-10`). Print-ready model
  names (for example "Log-logistic", "Generalized Gamma", "Hazards, 2
  knots").
* New `parameter_table()`: estimates, standard errors and confidence
  intervals on the natural and transformed scales.
* New `covariance_table()`: covariance matrices and their lower triangular
  Cholesky factors on the transformed and natural (delta method) scales (in
  version 1.0.0 the Cholesky factor was upper triangular).
* New `gof_table()`: log-likelihood, AIC and BIC sorted by AIC; the fixed
  shape of the Rayleigh model is not counted.
* New `predict_survival()`: survival, hazard or cumulative hazard at given
  times for each arm, with simulated confidence limits (`seed`, without
  changing the random number generator of the session) and the
  Kaplan-Meier estimate, `NA` after the largest observed time of the arm.
* New `survival_summary()`: median, mean and restricted mean survival.
* New `spline_knot_table()`: placement of the knots with the number at risk
  and the Kaplan-Meier estimate at each knot.
* New plots, replacing `plot_km_and_parametric()`,
  `plot_km_and_parametric_multiple()`, `plot_hazard_and_parametric()`,
  `plot_hazard_and_parametric_multiple()`, `plot_test_ph_assumption()`,
  `plot_schoenfeld()` and `plot_qq_survival()` of version 1.0.0. They draw
  from the fitted object, so plots with and without confidence bands do not
  refit the models, use print-ready model names with the AIC in the legend,
  and start the axes at (0, 0).
  * `plot_km()`: Kaplan-Meier curves with numbers at risk; each curve ends
    at the largest observed time of its arm (in version 1.0.0 the curves
    were extended to the longest arm).
  * `plot_empirical_hazard()`: smoothed empirical hazard without fitted
    models.
  * `plot_survival_fits()` and `plot_hazard_fits()`: Kaplan-Meier curve or
    smoothed empirical hazard with the fitted parametric and spline models,
    extrapolated to a time horizon, for one arm with several models or one
    model with several arms (independent or dependent fits).
  * `plot_distribution_diagnostics()`: log(-log S(t)), log-odds and probit
    of the Kaplan-Meier estimate against log time, and the log of the
    smoothed hazard against time (in version 1.0.0 the Gompertz plot used
    log(-log S(t)) against time, which is not linear for the Gompertz).
  * `plot_schoenfeld()`: the reference line is the estimated log hazard
    ratio, around which the scaled residuals are centred (it was 0 in
    version 1.0.0), and the regression on time is a natural spline with its
    confidence interval as in `survival::plot.cox.zph()`.
  * `plot_qq_survival()`: the axes are set from the data (they were fixed
    at 0 to 100 in version 1.0.0).
* New `write_survival_workbook()`: Excel workbook with one table per sheet
  (data summary, Kaplan-Meier table, hazard ratios, proportional hazards
  test, model status, fit statistics, parameters on the natural and
  transformed scales with the rows of their covariance matrix, covariance
  matrices and Cholesky factors, survival summary, predicted survival and
  knots), with print-ready model names.
* New `write_survival_report()`: Word report with the structure of the
  survival analysis report template (overview, methods, data overview,
  diagnostics, standard parametric models and spline-based extrapolation
  for each arm and for both arms), with figures without and with confidence
  intervals, smoothed hazard plots without fitted models, and print-ready
  model names. The styles come from a template in the package
  (`inst/templates/report_template.docx`: A4, Arial, numbered headings and
  tables in the style of the report template), or from a document given in
  `template`. Replaces the Word export of the Shiny application of version
  1.0.0.
* New example dataset `inst/extdata/fit_example.csv` (two arms of 300
  patients) and independent reference values for the tests
  (`inst/validation/python/expected_fits.py`, numpy and scipy).
* The license is GPL (>= 2) instead of MIT.
* broom, dplyr, magrittr, purrr, readr, rlang, tibble and tidyr are no
  longer used; flexsurv (>= 2.3), ggplot2 (>= 3.4.0), officer, openxlsx,
  patchwork and splines are imported; haven and readxl are suggested
  packages.

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
