# Hazard Ratios from Cox Proportional Hazards Models

Fits Cox proportional hazards models
([`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html))
with the treatment arm as the explanatory variable, and reports the
hazard ratio of each arm against the control arm (the first level of
`ARM`). The unadjusted model is always fitted; models adjusted for
`covariates`, stratified by `strata`, or both are added when these are
given.

## Usage

``` r
estimate_hazard_ratio(
  data,
  covariates = NULL,
  strata = NULL,
  ties = c("efron", "breslow", "exact"),
  conf_level = 0.95
)
```

## Arguments

- data:

  Data frame returned by
  [`processing_dataset`](https://gosukehommaEX.github.io/SurvFitSummary/reference/processing_dataset.md),
  with two or more treatment arms.

- covariates:

  Character vector with the names of covariate columns of `data` for the
  adjusted models. The default is `NULL`.

- strata:

  Character vector with the names of stratification columns of `data`
  for the stratified models. The default is `NULL`.

- ties:

  Method for tied event times: `"efron"` (the default of
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)),
  `"breslow"` or `"exact"`.

- conf_level:

  Confidence level, a single number in (0, 1). The default is 0.95.

## Value

A data frame with one row per model and comparison, and the columns
`model` (`"Unadjusted, unstratified"`, `"Unadjusted, stratified"`,
`"Adjusted, unstratified"` or `"Adjusted, stratified"`), `covariates`
and `strata` (the columns used, `"None"` if none), `comparison` (for
example `"Treatment vs Control"`), `hr`, `lcl`, `ucl` (Wald confidence
interval), `p_value` (Wald test), `n` and `events` (observations and
events used by the model).

## Examples

``` r
raw <- generate_dummy_survival_data(
  arm = c("Control", "Treatment"), n = c(100, 100),
  hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 123
)
dat <- processing_dataset(raw, column_arm = "ARM",
                          column_survtime = "SURVTIME",
                          column_cnsr = "CNSR", control_arm = "Control")
estimate_hazard_ratio(dat)
#>                      model covariates strata           comparison        hr
#> 1 Unadjusted, unstratified       None   None Treatment vs Control 0.6842886
#>        lcl       ucl    p_value   n events
#> 1 0.501417 0.9338553 0.01678775 200    167
```
