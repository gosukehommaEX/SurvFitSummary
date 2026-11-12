# Expected values: inst/validation/python/expected_fits.py (quantiles of the
# log event times, number at risk and Kaplan-Meier estimate with Greenwood
# standard error at each knot, comparisons on the log time scale).

test_that("knots and Kaplan-Meier estimates match the independent values", {
  x <- get_fits("hazard")
  out <- spline_knot_table(x)
  expected <- read_expected("knots")
  expect_named(out, c("fit_type", "arm", "n_knots", "knot", "quantile",
                      "log_time", "time", "n_risk", "surv_est", "surv_se"))
  expect_identical(nrow(out), nrow(expected))
  e <- lookup(out, expected, c("fit_type", "arm", "n_knots", "knot"))
  expect_equal(out$quantile, e$quantile, tolerance = 1e-10)
  expect_equal(out$log_time, e$log_time, tolerance = 1e-10)
  expect_equal(out$time, e$time, tolerance = 1e-10)
  expect_identical(out$n_risk, as.integer(e$n_risk))
  expect_equal(out$surv_est, e$surv_est, tolerance = 1e-10)
  expect_equal(out$surv_se, e$surv_se, tolerance = 1e-8)
  # The knots do not depend on the scale of the spline models
  expect_identical(spline_knot_table(get_fits("normal")), out)
})

test_that("a knot at an observed time counts that patient and event", {
  x <- get_fits("hazard")
  out <- spline_knot_table(x)
  d <- make_fit_data()
  ctrl <- d[d$ARM == "Control", ]
  median_knot <- out[out$arm == "Control" & out$n_knots == 1L, ]
  # With 271 events the median of the log event times is an observed one
  expect_true(any(log(ctrl$SURVTIME) == median_knot$log_time))
  at_knot <- log(ctrl$SURVTIME) == median_knot$log_time
  expect_identical(median_knot$n_risk,
                   sum(log(ctrl$SURVTIME) >= median_knot$log_time))
  expect_true(any(ctrl$EVENT[at_knot] == 1L))
})

test_that("no spline model gives an empty table", {
  out <- spline_knot_table(get_fits("simple"))
  expect_identical(nrow(out), 0L)
  expect_named(out, c("fit_type", "arm", "n_knots", "knot", "quantile",
                      "log_time", "time", "n_risk", "surv_est", "surv_se"))
  expect_error(spline_knot_table(list()), "fit_survival_models")
})
