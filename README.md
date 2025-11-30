# SurvFitSummary

<!-- badges: start -->
[![R-CMD-check](https://github.com/gosukehommaEX/SurvFitSummary/workflows/R-CMD-check/badge.svg)](https://github.com/gosukehommaEX/SurvFitSummary/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Overview

**SurvFitSummary** provides comprehensive tools for fitting parametric survival models and creating diagnostic visualizations for time-to-event data. The package supports multiple parametric distributions, automated model fitting with various configurations, and extensive diagnostic plotting capabilities.

### Key Features

- **Multiple Parametric Distributions**: Exponential, Weibull, lognormal, log-logistic, Gompertz, generalized gamma, and gamma distributions
- **Flexible Model Fitting**: Support for dependent (covariate-based) and independent (stratified) models
- **Rich Visualizations**:
  - Kaplan-Meier curves with parametric overlays
  - Hazard function estimation and plotting
  - Proportional hazards assumption testing (Schoenfeld residuals, statistical tests)
  - Q-Q plots for distribution assessment
- **Multiple Data Formats**: CSV, Excel (.xlsx, .xls), and SAS (.sas7bdat) datasets
- **Comprehensive Output**: Parameter estimates, goodness-of-fit statistics, landmark survival probabilities, and survival metrics (median, mean, RMST)

## Installation

You can install the development version of SurvFitSummary from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("gosukehommaEX/SurvFitSummary")
```

## Main Functions

### Data Processing

- `processing_dataset()`: Process and standardize survival data from various formats

### Dummy Data Generation

- `generate_dummy_survival_data()`: Generate simulated survival data for testing

### Model Fitting

- `fitting_surv_mod()`: Fit a single parametric distribution to survival data
- `fitting_surv_mod_multiple()`: Fit multiple distributions and compare models

### Diagnostic Plots

- `plot_km_and_parametric()`: Kaplan-Meier curves with parametric overlay (single distribution)
- `plot_km_and_parametric_multiple()`: Kaplan-Meier curves with multiple parametric overlays
- `plot_hazard_and_parametric()`: Hazard function plots (single distribution)
- `plot_hazard_and_parametric_multiple()`: Hazard function plots (multiple distributions)
- `plot_schoenfeld()`: Scaled Schoenfeld residuals for PH assumption testing
- `plot_test_ph_assumption()`: Complementary log-log and other PH diagnostic plots
- `plot_qq_survival()`: Q-Q plots for comparing survival distributions

## Usage Examples

### Example 1: Basic Workflow with Single Distribution

```r
library(SurvFitSummary)

# Step 1: Generate dummy survival data
dataset <- generate_dummy_survival_data(
  arm = c("Treatment", "Control"),
  n = c(100, 100),
  hazards = log(2) / c(20, 15),  # Median survival: 20 vs 15 months
  dropout_per_year = 0.1,
  seed = 123
)

# Step 2: Process the dataset
dataset_processed <- processing_dataset(
  dataset = dataset,
  column_arm = "ARM",
  column_survtime = "SURVTIME",
  column_cnsr = "CNSR",
  column_event = NULL
)

# Step 3: Fit Weibull distribution
model_results <- fitting_surv_mod(
  dataset = dataset_processed,
  distribution = "weibull",
  dependent = TRUE,
  control_arm = "Control",
  landmark_times = c(6, 12, 18, 24),
  time_horizon = 36
)

# View results
print(model_results$coef_estimates)
print(model_results$landmark_survival)
print(model_results$survival_metrics)

# Step 4: Visualize Kaplan-Meier with parametric overlay
km_plot <- plot_km_and_parametric(
  dataset = dataset_processed,
  distribution = "weibull",
  conf_int = TRUE
)

print(km_plot$plot)
print(km_plot$risktable)

# Step 5: Plot hazard functions
hazard_plot <- plot_hazard_and_parametric(
  dataset = dataset_processed,
  distribution = "weibull"
)

print(hazard_plot)
```

### Example 2: Multiple Distributions Comparison

```r
# Fit multiple distributions
multi_results <- fitting_surv_mod_multiple(
  dataset = dataset_processed,
  distribution = c("exp", "weibull", "lnorm", "gompertz"),
  dependent = "both",  # Fit both dependent and independent models
  control_arm = "Control",
  time_scale = "month",
  time_horizon = 60
)

# Compare goodness-of-fit
print(multi_results$gof_statistics)

# Visualize all distributions
km_multi_plot <- plot_km_and_parametric_multiple(
  dataset = dataset_processed,
  distribution = c("exp", "weibull", "lnorm", "gompertz"),
  conf_int = TRUE
)

print(km_multi_plot$plot)

# Compare hazard functions
hazard_multi_plot <- plot_hazard_and_parametric_multiple(
  dataset = dataset_processed,
  distribution = c("exp", "weibull", "lnorm", "gompertz")
)

print(hazard_multi_plot)
```

### Example 3: Proportional Hazards Assumption Testing

```r
# Schoenfeld residuals plot
schoenfeld_plot <- plot_schoenfeld(
  dataset = dataset_processed,
  control_arm = "Control"
)

print(schoenfeld_plot)

# Complementary log-log plot
ph_test_plot <- plot_test_ph_assumption(
  dataset = dataset_processed
)

print(ph_test_plot)
```

### Example 4: Q-Q Plot for Distribution Comparison

```r
qq_plot <- plot_qq_survival(
  dataset = dataset_processed,
  control_arm = "Control"
)

print(qq_plot)
```

### Example 5: Three-Arm Trial Analysis

```r
# Generate three-arm trial data
dataset3 <- generate_dummy_survival_data(
  arm = c("Treatment A", "Treatment B", "Control"),
  n = c(80, 90, 100),
  hazards = log(2) / c(25, 20, 15),
  dropout_per_year = 0.1,
  seed = 456
)

dataset3_processed <- processing_dataset(
  dataset = dataset3,
  column_arm = "ARM",
  column_survtime = "SURVTIME",
  column_cnsr = "CNSR",
  column_event = NULL
)

# Fit models
results3 <- fitting_surv_mod_multiple(
  dataset = dataset3_processed,
  distribution = c("weibull", "lnorm", "gompertz"),
  dependent = TRUE,
  control_arm = "Control",
  time_scale = "month",
  time_horizon = 48
)

# Visualize
km_plot3 <- plot_km_and_parametric_multiple(
  dataset = dataset3_processed,
  distribution = c("weibull", "lnorm", "gompertz"),
  conf_int = TRUE
)

print(km_plot3$plot)

# Q-Q plots (faceted by comparison)
qq_plot3 <- plot_qq_survival(
  dataset = dataset3_processed,
  control_arm = "Control"
)

print(qq_plot3)
```

### Example 6: Loading External Data

```r
# From CSV
dataset_csv <- processing_dataset(
  dataset = "path/to/survival_data.csv",
  column_arm = "TRT",
  column_survtime = "TIME",
  column_event = "STATUS"
)

# From Excel
dataset_excel <- processing_dataset(
  dataset = "path/to/survival_data.xlsx",
  column_arm = "ARM",
  column_survtime = "SURVTIME",
  column_cnsr = "CNSR"
)

# From SAS
dataset_sas <- processing_dataset(
  dataset = "path/to/adtte.sas7bdat",
  column_arm = "TRT01P",
  column_survtime = "AVAL",
  column_cnsr = "CNSR"
)

# Fit and visualize
results <- fitting_surv_mod(
  dataset = dataset_csv,
  distribution = "weibull",
  dependent = TRUE
)
```

### Example 7: Time Scale Conversion

```r
# Data in weeks, but want survival metrics in months
results_scaled <- fitting_surv_mod(
  dataset = dataset_processed,
  distribution = "weibull",
  dependent = TRUE,
  landmark_times = seq(4, 52, by = 4),  # Weekly landmarks
  time_horizon = 52,
  original_time_scale = 30.44 / 7  # Convert weeks to months
)

print(results_scaled$survival_metrics)  # Median, mean, RMST in months
```

### Example 8: Comprehensive Analysis Pipeline

```r
# Complete analysis workflow
library(SurvFitSummary)

# 1. Generate or load data
data <- generate_dummy_survival_data(
  arm = c("Drug A", "Drug B", "Placebo"),
  n = c(150, 150, 150),
  hazards = log(2) / c(24, 20, 16),
  dropout_per_year = 0.08,
  seed = 789
)

data_processed <- processing_dataset(
  dataset = data,
  column_arm = "ARM",
  column_survtime = "SURVTIME",
  column_cnsr = "CNSR"
)

# 2. Fit multiple distributions
distributions <- c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma")

fit_results <- fitting_surv_mod_multiple(
  dataset = data_processed,
  distribution = distributions,
  dependent = "both",
  control_arm = "Placebo",
  time_scale = "month",
  time_horizon = 60
)

# 3. Examine goodness-of-fit
print(fit_results$gof_statistics)

# 4. Visualize KM curves with best-fitting distribution
best_dist <- fit_results$gof_statistics %>%
  filter(model_type == "dependent") %>%
  arrange(AIC) %>%
  slice(1) %>%
  pull(distribution)

km_best <- plot_km_and_parametric(
  dataset = data_processed,
  distribution = best_dist,
  conf_int = TRUE
)

print(km_best$plot)

# 5. Compare all distributions
km_all <- plot_km_and_parametric_multiple(
  dataset = data_processed,
  distribution = distributions,
  conf_int = TRUE
)

print(km_all$plot)

# 6. Hazard function comparison
hazard_all <- plot_hazard_and_parametric_multiple(
  dataset = data_processed,
  distribution = distributions
)

print(hazard_all)

# 7. Test proportional hazards assumption
schoenfeld <- plot_schoenfeld(
  dataset = data_processed,
  control_arm = "Placebo"
)

print(schoenfeld)

ph_diagnostic <- plot_test_ph_assumption(
  dataset = data_processed
)

print(ph_diagnostic)

# 8. Q-Q plots for distribution comparison
qq <- plot_qq_survival(
  dataset = data_processed,
  control_arm = "Placebo"
)

print(qq)

# 9. Extract key results
cat("\n=== Parameter Estimates ===\n")
print(fit_results$coef_estimates %>% filter(distribution == best_dist))

cat("\n=== Landmark Survival Probabilities ===\n")
print(fit_results$landmark_survival %>% filter(distribution == best_dist))

cat("\n=== Survival Metrics ===\n")
print(fit_results$survival_metrics %>% filter(distribution == best_dist))

cat("\n=== Cox Proportional Hazards Results ===\n")
print(fit_results$cox_ph_results)
```

## Supported Distributions

| Distribution | Code | Typical Use Case |
|--------------|------|------------------|
| Exponential | `"exp"` | Constant hazard over time |
| Weibull | `"weibull"` | Monotonic increasing/decreasing hazard |
| Lognormal | `"lnorm"` | Non-monotonic hazard (peak then decline) |
| Log-logistic | `"llogis"` | Non-monotonic hazard, crossing survival curves |
| Gompertz | `"gompertz"` | Exponentially increasing/decreasing hazard |
| Generalized Gamma | `"gengamma"` | Highly flexible, nests several distributions |
| Gamma | `"gamma"` | Monotonic hazard with initial shoulder |

## Output Components

### Model Fitting Functions

- `coef_estimates`: Parameter estimates with standard errors and confidence intervals
- `varcovmat`: Variance-covariance matrix in long format
- `cholesky`: Cholesky decomposition of variance-covariance matrix
- `gof_statistics`: AIC and BIC for model comparison
- `landmark_survival`: Survival probabilities at specified time points
- `survival_metrics`: Median, mean, and restricted mean survival time (RMST)
- `cox_ph_results`: Cox proportional hazards model results (for comparison)

### Plotting Functions

- `plot`: ggplot2 object ready for display or further customization
- `risktable`: Data frame with number at risk at regular intervals

## Tips and Best Practices

1. **Model Selection**: Use AIC/BIC from `gof_statistics` to compare distributions
2. **Visual Inspection**: Always examine KM overlay plots and hazard plots
3. **PH Assumption**: Check Schoenfeld residuals before interpreting Cox model results
4. **Time Scale**: Ensure consistent time units across analysis
5. **Sample Size**: Generalized gamma requires larger sample sizes for stable estimation
6. **Censoring**: Functions handle both CNSR (1=censored) and EVENT (1=event) indicators

## Author

**Gosuke Homma**

- Email: my.name.is.gosuke@gmail.com
- ORCID: [0000-0002-6854-5627](https://orcid.org/0000-0002-6854-5627)
- GitHub: [@gosukehommaEX](https://github.com/gosukehommaEX)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Citation

If you use SurvFitSummary in your research, please cite:

```
Homma, G. (2025). SurvFitSummary: Fit and Visualize Parametric Survival Models 
with Diagnostic Plots. R package version 1.0.0. 
https://github.com/gosukehommaEX/SurvFitSummary
```

## Acknowledgments

This package builds upon the excellent work of:

- `survival` package for core survival analysis functionality
- `flexsurv` package for parametric survival model fitting
- `bshazard` package for non-parametric hazard estimation
- `ggplot2` package for visualization

## Bug Reports and Feature Requests

Please report bugs and request features at:
https://github.com/gosukehommaEX/SurvFitSummary/issues

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
