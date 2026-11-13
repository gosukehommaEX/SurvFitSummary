# Grambsch-Therneau Test of Proportional Hazards

For each treatment arm other than the control arm (the first level of
`ARM`), fits a Cox model to the data of the two arms and tests the
proportional hazards assumption with the score test of Grambsch and
Therneau (1994) as implemented in
[`survival::cox.zph()`](https://rdrr.io/pkg/survival/man/cox.zph.html).

## Usage

``` r
test_proportional_hazards(data, transform = c("km", "rank", "identity", "log"))
```

## Arguments

- data:

  Data frame returned by
  [`processing_dataset`](https://gosukehommaEX.github.io/SurvFitSummary/reference/processing_dataset.md),
  with two or more treatment arms.

- transform:

  Transformation of time for the test: `"km"` (the default of
  [`survival::cox.zph()`](https://rdrr.io/pkg/survival/man/cox.zph.html)),
  `"rank"`, `"identity"` or `"log"`.

## Value

A data frame with the columns `comparison`, `log_hr` (Cox estimate of
the log hazard ratio), `chisq`, `df` and `p_value`.

## References

Grambsch PM, Therneau TM (1994). Proportional hazards tests and
diagnostics based on weighted residuals. Biometrika 81(3):515-526.
[doi:10.1093/biomet/81.3.515](https://doi.org/10.1093/biomet/81.3.515)

## Examples

``` r
raw <- generate_dummy_survival_data(
  arm = c("Control", "Treatment"), n = c(100, 100),
  hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 123
)
dat <- processing_dataset(raw, column_arm = "ARM",
                          column_survtime = "SURVTIME",
                          column_cnsr = "CNSR", control_arm = "Control")
test_proportional_hazards(dat)
#>             comparison     log_hr     chisq df   p_value
#> 1 Treatment vs Control -0.3793755 0.0198715  1 0.8878966
```
