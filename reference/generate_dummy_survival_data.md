# Generate Dummy Survival Data

Simulates exponential event times for each treatment arm and exponential
dropout times, and returns the observed time and the censoring
indicator. The random number generator state of the session is restored
afterwards.

## Usage

``` r
generate_dummy_survival_data(arm, n, hazards, dropout_per_year, seed)
```

## Arguments

- arm:

  Character vector with the names of the treatment arms.

- n:

  Integer vector with the number of subjects in each arm.

- hazards:

  Numeric vector with the hazard rate of the event in each arm, per
  month.

- dropout_per_year:

  Proportion of subjects who drop out within one year, a single number
  in \[0, 1). The dropout hazard per month is
  `-log(1 - dropout_per_year) / 12`.

- seed:

  Single integer, the random seed.

## Value

A data frame with the columns `SUBJID`, `ARM`, `SURVTIME` (minimum of
the event and dropout times, in months) and `CNSR` (1 = censored by
dropout, 0 = event).

## Examples

``` r
dat <- generate_dummy_survival_data(
  arm = c("Control", "Treatment"), n = c(100, 100),
  hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 123
)
head(dat)
#>   SUBJID     ARM   SURVTIME CNSR
#> 1    ID1 Control 14.6022193    0
#> 2    ID2 Control  9.4265013    1
#> 3    ID3 Control 23.0090504    0
#> 4    ID4 Control  0.5466780    0
#> 5    ID5 Control  0.9731436    0
#> 6    ID6 Control  5.4793768    0
```
