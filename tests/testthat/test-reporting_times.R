# Expected values: the largest observed time of fit_example.csv is 70.2602
# months, so the spacing of the time axis is 12 months (70.2602 / 6 > 10 and
# 70.2602 / 12 <= 10) and, read as days, 14 days (70.2602 / 7 > 10 and
# 70.2602 / 14 <= 10). The default landmarks are 0.5, 1, 2, 5 and 10 years.

test_that("default reporting times", {
  x <- get_fits("report")
  out <- reporting_times(x)
  expect_named(out, c("interval", "landmarks", "time_horizon", "time_unit"))
  expect_identical(out$interval, 12)
  expect_identical(out$landmarks$time, c(6, 12, 24, 60, 120))
  expect_identical(out$landmarks$label,
                   c("6 months", "1 year", "2 years", "5 years", "10 years"))
  expect_identical(out$time_horizon, 120)
  expect_identical(out$time_unit, "months")
})

test_that("reporting times given by the user", {
  x <- get_fits("report")
  out <- reporting_times(x, interval = 6, landmarks = c(48, 36, 48),
                         time_horizon = 60)
  expect_identical(out$interval, 6)
  expect_identical(out$landmarks$time, c(36, 48))
  expect_identical(out$landmarks$label, c("36 months", "48 months"))
  expect_identical(out$time_horizon, 60)
  # Without a time horizon, the last landmark time
  expect_identical(reporting_times(x, landmarks = c(48, 36))$time_horizon,
                   48)
})

test_that("the defaults follow the time unit", {
  x <- fit_survival_models(make_fit_data(), distributions = "exp",
                           spline_scales = NULL, fit_type = "independent",
                           time_unit = "days")
  out <- reporting_times(x)
  expect_identical(out$interval, 14)
  expect_equal(out$landmarks$time,
               c(182.625, 365.25, 730.5, 1826.25, 3652.5))
  expect_identical(out$landmarks$label,
                   c("6 months", "1 year", "2 years", "5 years", "10 years"))
  expect_equal(out$time_horizon, 3652.5)
  expect_identical(out$time_unit, "days")
  expect_identical(reporting_times(x, landmarks = 30)$landmarks$label,
                   "30 days")
})

test_that("the workbook and the report use the same reporting times", {
  x <- get_fits("report")
  out <- reporting_times(x, interval = 6, landmarks = c(48, 36))
  s <- .report_settings(x, 6, c(48, 36), NULL, NULL)
  expect_identical(s$interval, out$interval)
  expect_identical(s$landmarks, out$landmarks)
  expect_identical(s$time_horizon, out$time_horizon)
  expect_identical(s$unit, out$time_unit)
})

test_that("wrong arguments stop with a message", {
  x <- get_fits("report")
  expect_error(reporting_times(list()), "fit_survival_models")
  expect_error(reporting_times(x, interval = 0), "'interval'")
  expect_error(reporting_times(x, interval = c(6, 12)), "'interval'")
  expect_error(reporting_times(x, landmarks = c(12, -1)), "'landmarks'")
  expect_error(reporting_times(x, landmarks = numeric(0)), "'landmarks'")
  expect_error(reporting_times(x, time_horizon = NA_real_), "'time_horizon'")
})
