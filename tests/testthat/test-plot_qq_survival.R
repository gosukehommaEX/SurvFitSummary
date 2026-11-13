# Expected values: inst/validation/python/expected_fits.py (Kaplan-Meier
# quantiles with the midpoint rule of survival::quantile.survfit,
# expected/qq.csv).

test_that("the quantiles match the independent Kaplan-Meier quantiles", {
  d <- make_fit_data()
  p <- plot_qq_survival(d)
  expect_s3_class(p, "ggplot")
  points <- layer_input(p, "GeomPoint")
  e <- read_expected("qq")
  expect_equal(points$prob, seq(0.05, 0.95, by = 0.05))
  expect_equal(points$control, e$quantile[e$arm == "Control"],
               tolerance = 1e-10)
  expect_equal(points$arm, e$quantile[e$arm == "Treatment"],
               tolerance = 1e-10)
  lines <- layer_input(p, "GeomAbline")
  fit <- stats::lm(points$arm ~ points$control)
  expect_equal(lines$intercept, c(0, unname(stats::coef(fit)[1])))
  expect_equal(lines$slope, c(1, unname(stats::coef(fit)[2])))
  # Axes from the data, the same on both axes
  limit <- 1.05 * max(e$quantile)
  r <- panel_ranges(p)
  expect_equal(r$x, c(0, limit))
  expect_equal(r$y, c(0, limit))
})

test_that("quantiles not reached are left out", {
  d <- make_fit_data()
  # The Treatment curve stays above 0.01, so its 99% quantile is not reached
  p <- plot_qq_survival(d, probs = c(0.5, 0.99))
  points <- layer_input(p, "GeomPoint")
  expect_identical(points$prob, 0.5)
  expect_identical(layer_input(p, "GeomAbline")$type,
                   "Identity line (y = x)")
})

test_that("invalid inputs are rejected", {
  d <- make_fit_data()
  expect_error(plot_qq_survival(d[d$ARM == "Control", ]),
               "At least 2 treatment arms")
  expect_error(plot_qq_survival(d, probs = 1), "'probs'")
})

test_that("the legends are stacked and placed with legend_position", {
  d <- make_fit_data()
  p <- plot_qq_survival(d)
  expect_identical(p$theme$legend.position, "bottom")
  expect_identical(p$theme$legend.box, "vertical")
  p <- plot_qq_survival(d, legend_position = "left")
  expect_identical(p$theme$legend.position, "left")
  expect_error(plot_qq_survival(d, legend_position = "middle"))
})
