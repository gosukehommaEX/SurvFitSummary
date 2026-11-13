# Print the Survival Models Fitted by fit_survival_models()

Prints the time unit, the fits (arms), the numbers of models fitted
without problems, with warnings and failed, and the warnings and errors
of the models that had any.

## Usage

``` r
# S3 method for class 'survfitsummary'
print(x, ...)
```

## Arguments

- x:

  Object returned by
  [`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md).

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(dat, distributions = "exp",
                            spline_scales = NULL)
print(fits)
#> Survival models fitted with SurvFitSummary (time unit: months)
#> Fits: Control (independent), Treatment (independent), Both arms (dependent)
#> 3 models: 3 without problems, 0 with warnings, 0 failed
```
