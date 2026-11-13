# Kaplan-Meier Summary Table at Regular Times

Reports the Kaplan-Meier estimate of each treatment arm at the times
`0, interval, 2 * interval, ...` up to the largest observed time of the
arm, with the number at risk, the numbers of events and censored
observations since the previous time, the standard error (Greenwood) and
the pointwise confidence interval of
[`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html).

## Usage

``` r
km_table(
  data,
  interval,
  conf_level = 0.95,
  conf_type = c("log", "log-log", "plain")
)
```

## Arguments

- data:

  Data frame returned by
  [`processing_dataset`](https://gosukehommaEX.github.io/SurvFitSummary/reference/processing_dataset.md).

- interval:

  Positive number, the spacing of the reporting times in the time unit
  of the data.

- conf_level:

  Confidence level, a single number in (0, 1). The default is 0.95.

- conf_type:

  Transformation for the pointwise confidence intervals: `"log"` (the
  default of
  [`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)),
  `"log-log"` or `"plain"`.

## Value

A data frame with the columns `ARM`, `time`, `n_risk` (number with an
observed time at or after `time`), `n_event` and `n_censor` (numbers of
events and censored observations after the previous reporting time and
up to `time`), `surv` (Kaplan-Meier estimate), `se` (standard error of
`surv`), `lcl` and `ucl` (confidence limits; both equal to 1 while
`surv` is 1).

## Examples

``` r
raw <- generate_dummy_survival_data(
  arm = c("Control", "Treatment"), n = c(100, 100),
  hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 123
)
dat <- processing_dataset(raw, column_arm = "ARM",
                          column_survtime = "SURVTIME",
                          column_cnsr = "CNSR")
km_table(dat, interval = 6)
#>          ARM time n_risk n_event n_censor       surv         se         lcl
#> 1    Control    0    100       0        0 1.00000000 0.00000000 1.000000000
#> 2    Control    6     71      28        1 0.71946667 0.04498289 0.636489718
#> 3    Control   12     51      16        4 0.55101488 0.05051142 0.460398459
#> 4    Control   18     39      11        1 0.43032952 0.05090504 0.341279201
#> 5    Control   24     24      14        1 0.27541089 0.04647876 0.197847718
#> 6    Control   30     12      11        1 0.14829817 0.03770759 0.090095235
#> 7    Control   36      8       3        1 0.10785322 0.03389124 0.058258049
#> 8    Control   42      5       1        2 0.09437156 0.03222491 0.048326209
#> 9    Control   48      2       3        0 0.03774863 0.02436470 0.010653601
#> 10   Control   54      2       0        0 0.03774863 0.02436470 0.010653601
#> 11   Control   60      2       0        0 0.03774863 0.02436470 0.010653601
#> 12   Control   66      1       0        1 0.03774863 0.02436470 0.010653601
#> 13   Control   72      1       0        0 0.03774863 0.02436470 0.010653601
#> 14   Control   78      1       0        0 0.03774863 0.02436470 0.010653601
#> 15   Control   84      1       0        0 0.03774863 0.02436470 0.010653601
#> 16   Control   90      1       0        0 0.03774863 0.02436470 0.010653601
#> 17   Control   96      1       0        0 0.03774863 0.02436470 0.010653601
#> 18   Control  102      1       0        0 0.03774863 0.02436470 0.010653601
#> 19   Control  108      1       0        0 0.03774863 0.02436470 0.010653601
#> 20 Treatment    0    100       0        0 1.00000000 0.00000000 1.000000000
#> 21 Treatment    6     74      18        8 0.81348007 0.03975642 0.739174497
#> 22 Treatment   12     57      12        5 0.67411445 0.04928595 0.584117847
#> 23 Treatment   18     40      14        3 0.50343091 0.05399588 0.407984423
#> 24 Treatment   24     35       4        1 0.45240751 0.05422099 0.357695431
#> 25 Treatment   30     30       5        0 0.38777787 0.05362829 0.295709233
#> 26 Treatment   36     21       9        0 0.27144451 0.04961692 0.189710448
#> 27 Treatment   42     14       5        2 0.20225277 0.04568817 0.129900583
#> 28 Treatment   48      8       6        0 0.11557301 0.03737864 0.061314143
#> 29 Treatment   54      7       1        0 0.10112638 0.03538813 0.050932537
#> 30 Treatment   60      4       2        1 0.07223313 0.03061191 0.031477706
#> 31 Treatment   66      4       0        0 0.07223313 0.03061191 0.031477706
#> 32 Treatment   72      4       0        0 0.07223313 0.03061191 0.031477706
#> 33 Treatment   78      2       2        0 0.03611657 0.02367222 0.009995293
#> 34 Treatment   84      2       0        0 0.03611657 0.02367222 0.009995293
#> 35 Treatment   90      1       1        0 0.01805828 0.01741104 0.002728875
#> 36 Treatment   96      1       0        0 0.01805828 0.01741104 0.002728875
#> 37 Treatment  102      1       0        0 0.01805828 0.01741104 0.002728875
#> 38 Treatment  108      1       0        0 0.01805828 0.01741104 0.002728875
#>          ucl
#> 1  1.0000000
#> 2  0.8132610
#> 3  0.6594666
#> 4  0.5426158
#> 5  0.3833815
#> 6  0.2441011
#> 7  0.1996688
#> 8  0.1842891
#> 9  0.1337537
#> 10 0.1337537
#> 11 0.1337537
#> 12 0.1337537
#> 13 0.1337537
#> 14 0.1337537
#> 15 0.1337537
#> 16 0.1337537
#> 17 0.1337537
#> 18 0.1337537
#> 19 0.1337537
#> 20 1.0000000
#> 21 0.8952552
#> 22 0.7779771
#> 23 0.6212068
#> 24 0.5721978
#> 25 0.5085119
#> 26 0.3883925
#> 27 0.3149038
#> 28 0.2178473
#> 29 0.2007861
#> 30 0.1657562
#> 31 0.1657562
#> 32 0.1657562
#> 33 0.1305021
#> 34 0.1305021
#> 35 0.1195004
#> 36 0.1195004
#> 37 0.1195004
#> 38 0.1195004
```
