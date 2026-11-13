# Prepare a Dataset for Survival Analysis

Takes a data frame or reads a file, renames the columns of the treatment
arm, the survival time and the event or censoring indicator to `ARM`,
`SURVTIME` and `EVENT`, keeps the requested covariate and stratification
columns, and checks the values. The treatment arm is returned as a
factor whose first level is the control arm, so that every function of
the package uses the same reference arm.

## Usage

``` r
processing_dataset(
  dataset,
  column_arm,
  column_survtime,
  column_cnsr = NULL,
  column_event = NULL,
  column_covariates = NULL,
  column_strata = NULL,
  control_arm = NULL
)
```

## Arguments

- dataset:

  A data frame, or a single character string with the path of a CSV
  (`.csv`), Excel (`.xlsx`, `.xls`), SAS (`.sas7bdat`) or SAS transport
  (`.xpt`) file. Excel files need the package readxl and SAS files the
  package haven.

- column_arm:

  Name of the column with the treatment arm.

- column_survtime:

  Name of the column with the survival time. The times must be positive.

- column_cnsr:

  Name of the column with the censoring indicator (1 = censored, 0 =
  event). The default is `NULL`.

- column_event:

  Name of the column with the event indicator (1 = event, 0 = censored).
  The default is `NULL`. When both `column_cnsr` and `column_event` are
  given, `column_event` is used.

- column_covariates:

  Character vector with the names of the columns kept as covariates for
  adjusted hazard ratios. The default is `NULL` (none).

- column_strata:

  Character vector with the names of the columns kept as stratification
  factors. The default is `NULL` (none).

- control_arm:

  Value of the treatment arm used as the control (reference) arm. The
  default is `NULL`, which uses the first level when the arm column is a
  factor and the first value in sorted order otherwise.

## Value

A data frame with the columns `ARM` (factor with the control arm as the
first level), `SURVTIME` (numeric), `EVENT` (integer, 1 = event, 0 =
censored), followed by the covariate and stratification columns with
their original names. Rows with a missing value in any of these columns
are removed with a warning.

## Examples

``` r
raw <- generate_dummy_survival_data(
  arm = c("Control", "Treatment"), n = c(50, 50),
  hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 1
)
dat <- processing_dataset(raw, column_arm = "ARM",
                          column_survtime = "SURVTIME",
                          column_cnsr = "CNSR", control_arm = "Control")
head(dat)
#>       ARM  SURVTIME EVENT
#> 1 Control 13.073965     1
#> 2 Control 20.457002     1
#> 3 Control  2.522524     1
#> 4 Control  2.420183     1
#> 5 Control  7.549369     1
#> 6 Control 50.118681     1
levels(dat$ARM)
#> [1] "Control"   "Treatment"
```
