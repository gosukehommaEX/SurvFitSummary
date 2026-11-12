test_that("the print method summarises the fits and returns x", {
  x <- get_fits("simple")
  out <- capture.output(res <- print(x))
  expect_identical(res, x)
  expect_match(out[1], "time unit: months")
  expect_match(out[2], "Control \\(independent\\), Treatment \\(independent\\), Both arms \\(dependent\\)")
  expect_match(out[3], "24 models: 24 without problems, 0 with warnings, 0 failed")
  expect_length(out, 3)
})

test_that("models with warnings or errors are listed", {
  d <- make_fit_data()
  d$EVENT[d$ARM == "Treatment"] <- 0L
  d$EVENT[d$ARM == "Treatment"][1] <- 1L
  x <- fit_survival_models(d, distributions = "exp",
                           spline_scales = "hazard", spline_knots = 1,
                           fit_type = "independent")
  out <- capture.output(print(x))
  expect_match(out[3], "4 models: 3 without problems, 0 with warnings, 1 failed")
  expect_true(any(grepl("Models with warnings or errors", out)))
  expect_true(any(grepl("two distinct event times", out)))
})
