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
  lines <- layer_input(p, "GeomSegment")
  fit <- stats::lm(points$arm ~ points$control)
  expect_equal(lines$intercept, c(0, unname(stats::coef(fit)[1])))
  expect_equal(lines$slope, c(1, unname(stats::coef(fit)[2])))
  # The lines are segments across the plot, so that their keys in the
  # legend are horizontal lines (those of geom_abline() are diagonal)
  limit <- 1.05 * max(e$quantile)
  expect_equal(lines$x, c(0, 0))
  expect_equal(lines$xend, c(limit, limit))
  expect_equal(lines$y, lines$intercept)
  expect_equal(lines$yend, lines$intercept + lines$slope * limit)
  # Axes from the data, the same on both axes
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
  expect_identical(layer_input(p, "GeomSegment")$type,
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

test_that("colours, line types and text sizes are set by the arguments", {
  d <- make_fit_data()
  default <- plot_qq_survival(d)
  expect_identical(group_values(default, "GeomSegment", "colour"),
                   c("grey30", "grey30"))
  expect_identical(group_values(default, "GeomSegment", "linetype"),
                   c("dashed", "solid"))
  p <- plot_qq_survival(d, colours = c(quantiles = "black",
                                       least_squares = "red"),
                        linetypes = c(identity = "dotted"),
                        text_sizes = c(axis_title = 11, legend = 10))
  expect_identical(group_values(p, "GeomSegment", "colour"),
                   c("grey30", "red"))
  expect_identical(group_values(p, "GeomSegment", "linetype"),
                   c("dotted", "solid"))
  expect_identical(group_values(p, "GeomPoint", "colour"), "black")
  expect_identical(theme_size(p, "axis.title.y"), 11)
  expect_identical(theme_size(p, "legend.text"), 10)
  # An unnamed vector gives the colours in the order of the elements
  p <- plot_qq_survival(d, colours = c("blue", "green4", "orange"))
  expect_identical(group_values(p, "GeomSegment", "colour"),
                   c("blue", "green4"))
  expect_identical(group_values(p, "GeomPoint", "colour"), "orange")
})

test_that("invalid colours, line types and text sizes are rejected", {
  d <- make_fit_data()
  expect_error(plot_qq_survival(d, linetypes = c(quantiles = "dashed")),
               "names of 'linetypes'")
  expect_error(plot_qq_survival(d, colours = c(identity = "#12")),
               "invalid colours")
})
