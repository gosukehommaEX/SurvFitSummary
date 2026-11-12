test_that("columns are renamed and EVENT is derived from CNSR", {
  raw <- make_test_raw()
  out <- make_test_data()
  expect_named(out, c("ARM", "SURVTIME", "EVENT", "AGE", "SEX", "REGION"))
  expect_identical(out$SURVTIME, raw$AVAL)
  expect_identical(out$EVENT, as.integer(1 - raw$CNSR))
  expect_identical(out$AGE, raw$AGE)
  expect_identical(out$REGION, raw$REGION)
})

test_that("the control arm is the first level of ARM", {
  expect_identical(levels(make_test_data("Control")$ARM),
                   c("Control", "Drug"))
  expect_identical(levels(make_test_data("Drug")$ARM), c("Drug", "Control"))
  out <- processing_dataset(make_test_raw(), column_arm = "TRT01P",
                            column_survtime = "AVAL", column_cnsr = "CNSR")
  expect_identical(levels(out$ARM), c("Control", "Drug"))
  raw <- make_test_raw()
  raw$TRT01P <- factor(raw$TRT01P, levels = c("Drug", "Control"))
  out <- processing_dataset(raw, column_arm = "TRT01P",
                            column_survtime = "AVAL", column_cnsr = "CNSR")
  expect_identical(levels(out$ARM), c("Drug", "Control"))
})

test_that("the event column takes priority over the censoring column", {
  raw <- make_test_raw()
  raw$EVT <- 1 - raw$CNSR
  raw$CNSR <- 0
  out <- processing_dataset(raw, column_arm = "TRT01P",
                            column_survtime = "AVAL", column_cnsr = "CNSR",
                            column_event = "EVT")
  expect_identical(out$EVENT, as.integer(raw$EVT))
  raw$EVT <- raw$EVT == 1
  out <- processing_dataset(raw, column_arm = "TRT01P",
                            column_survtime = "AVAL", column_event = "EVT")
  expect_identical(out$EVENT, as.integer(raw$EVT))
})

test_that("rows with missing values are removed with a warning", {
  raw <- make_test_raw()
  raw$AVAL[3] <- NA
  raw$AGE[5] <- NA
  expect_warning(
    out <- processing_dataset(raw, column_arm = "TRT01P",
                              column_survtime = "AVAL", column_cnsr = "CNSR",
                              column_covariates = "AGE"),
    "2 row\\(s\\) with missing values"
  )
  expect_identical(nrow(out), 23L)
})

test_that("invalid inputs are rejected", {
  raw <- make_test_raw()
  expect_error(processing_dataset(raw, "TRT01P", "AVAL"),
               "Either 'column_cnsr' or 'column_event'")
  expect_error(processing_dataset(raw, "TRT", "AVAL", column_cnsr = "CNSR"),
               "not found")
  expect_error(processing_dataset(raw, "TRT01P", "AVAL", column_cnsr = "CNSR",
                                  control_arm = "Placebo"),
               "'control_arm' must be one of")
  expect_error(processing_dataset(raw, "TRT01P", "AVAL", column_cnsr = "AVAL"),
               "different columns")
  raw0 <- raw
  raw0$AVAL[1] <- 0
  expect_error(processing_dataset(raw0, "TRT01P", "AVAL",
                                  column_cnsr = "CNSR"),
               "must be positive: 1 row")
  raw2 <- raw
  raw2$CNSR[1] <- 2
  expect_error(processing_dataset(raw2, "TRT01P", "AVAL",
                                  column_cnsr = "CNSR"),
               "only 0 and 1")
  raw3 <- raw
  raw3$EVENT <- 1
  expect_error(processing_dataset(raw3, "TRT01P", "AVAL",
                                  column_cnsr = "CNSR",
                                  column_covariates = "EVENT"),
               "cannot be named")
  expect_error(processing_dataset(list(a = 1), "a", "b", column_cnsr = "c"),
               "must be a data frame")
})

test_that("CSV files are read and unsupported files are rejected", {
  raw <- make_test_raw()
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(raw, path, row.names = FALSE)
  out <- processing_dataset(path, column_arm = "TRT01P",
                            column_survtime = "AVAL", column_cnsr = "CNSR")
  expect_identical(out$SURVTIME, raw$AVAL)
  expect_identical(out$EVENT, as.integer(1 - raw$CNSR))
  txt <- tempfile(fileext = ".txt")
  on.exit(unlink(txt), add = TRUE)
  file.create(txt)
  expect_error(processing_dataset(txt, "TRT01P", "AVAL",
                                  column_cnsr = "CNSR"),
               "Unsupported file format")
  expect_error(processing_dataset(tempfile(fileext = ".csv"), "TRT01P",
                                  "AVAL", column_cnsr = "CNSR"),
               "File not found")
})

test_that("Excel files are read", {
  skip_if_not_installed("readxl")
  # inst/extdata/survival_example.xlsx holds make_test_raw(), written with
  # Python openpyxl
  raw <- make_test_raw()
  path <- system.file("extdata", "survival_example.xlsx",
                      package = "SurvFitSummary")
  out <- processing_dataset(path, column_arm = "TRT01P",
                            column_survtime = "AVAL", column_cnsr = "CNSR",
                            column_covariates = "AGE")
  expect_identical(out$SURVTIME, raw$AVAL)
  expect_identical(out$AGE, raw$AGE)
})

test_that("SAS transport files are read and labels are dropped", {
  skip_if_not_installed("haven")
  raw <- make_test_raw()
  path <- tempfile(fileext = ".xpt")
  on.exit(unlink(path), add = TRUE)
  haven::write_xpt(raw, path)
  out <- processing_dataset(path, column_arm = "TRT01P",
                            column_survtime = "AVAL", column_cnsr = "CNSR",
                            column_strata = "REGION")
  expect_identical(out$SURVTIME, raw$AVAL)
  expect_identical(out$REGION, raw$REGION)
  expect_identical(levels(out$ARM), c("Control", "Drug"))
})
