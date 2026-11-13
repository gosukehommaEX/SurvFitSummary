# Summarise Survival Data by Treatment Arm

Counts the patients, events and censored observations in each treatment
arm and estimates the median survival time with its confidence interval
from the Kaplan-Meier estimate
([`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)).
The median is the first time at which the Kaplan-Meier curve is at or
below 0.5, and its confidence limits are the first times at which the
pointwise lower and upper confidence limits of the curve are at or below
0.5.

## Usage

``` r
summarise_data(
  data,
  conf_level = 0.95,
  conf_type = c("log", "log-log", "plain")
)
```

## Arguments

- data:

  Data frame returned by
  [`processing_dataset`](https://gosukehommaEX.github.io/SurvFitSummary/reference/processing_dataset.md).

- conf_level:

  Confidence level, a single number in (0, 1). The default is 0.95.

- conf_type:

  Transformation for the pointwise confidence intervals of the survival
  probability: `"log"` (the default of
  [`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)),
  `"log-log"` or `"plain"`.

## Value

A data frame with one row per treatment arm and, when there are two or
more arms, a last row `"Overall"` for the pooled data. The columns are
`ARM`, `n`, `events`, `censored`, `median`, `median_lcl`, `median_ucl`
(`NA` when the curve or the confidence limit does not reach 0.5) and
`max_followup` (largest observed time).

## Examples

``` r
raw <- generate_dummy_survival_data(
  arm = c("Control", "Treatment"), n = c(100, 100),
  hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 123
)
dat <- processing_dataset(raw, column_arm = "ARM",
                          column_survtime = "SURVTIME",
                          column_cnsr = "CNSR")
summarise_data(dat)
#>         ARM   n events censored   median median_lcl median_ucl max_followup
#> 1   Control 100     87       13 14.71179   10.35018   19.30203     111.7189
#> 2 Treatment 100     80       20 18.14788   14.77999   30.69736     113.3762
#> 3   Overall 200    167       33 17.23933   14.59690   20.22997     113.3762
```
