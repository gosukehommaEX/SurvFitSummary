# The smoothed hazards are those of estimate_empirical_hazard() (tested in
# test-estimate_empirical_hazard.R).

test_that("the curves are the smoothed empirical hazards of the arms", {
  d <- make_fit_data()
  p <- plot_empirical_hazard(d)
  expect_s3_class(p, "ggplot")
  line <- layer_input(p, "GeomLine")
  hz <- estimate_empirical_hazard(d)
  expect_equal(line$hazard, hz$hazard)
  expect_identical(as.character(line$ARM), hz$ARM)
  expect_length(layer_index(p, "GeomRibbon"), 0)
  r <- panel_ranges(p)
  expect_equal(r$x, c(0, max(d$SURVTIME)))
  expect_equal(r$y, c(0, 1.05 * max(hz$hazard)))
})

test_that("one arm with its confidence band", {
  d <- make_fit_data()
  p <- plot_empirical_hazard(d, arms = "Control", conf_int = TRUE)
  band <- layer_input(p, "GeomRibbon")
  expect_identical(unique(as.character(band$ARM)), "Control")
  expect_true(all(band$lcl <= band$hazard & band$hazard <= band$ucl))
  expect_equal(panel_ranges(p)$y, c(0, 1.05 * max(band$ucl)))
  expect_equal(panel_ranges(p)$x,
               c(0, max(d$SURVTIME[d$ARM == "Control"])))
})

test_that("an arm keeps its colour when plotted alone", {
  d <- make_fit_data()
  p <- plot_empirical_hazard(d, arms = "Treatment")
  line <- ggplot2::layer_data(p, layer_index(p, "GeomLine"))
  expect_identical(unique(line$colour),
                   unname(.plot_colours(arms = levels(d$ARM))["Treatment"]))
})

test_that("invalid inputs are rejected", {
  d <- make_fit_data()
  expect_error(plot_empirical_hazard(make_test_raw()), "processing_dataset")
  expect_error(plot_empirical_hazard(d, arms = "Placebo"), "'arms'")
  expect_error(plot_empirical_hazard(d, conf_level = 0), "'conf_level'")
})

test_that("the legend is placed with legend_position", {
  d <- make_fit_data()
  expect_identical(plot_empirical_hazard(d)$theme$legend.position, "bottom")
  p <- plot_empirical_hazard(d, legend_position = "none")
  expect_identical(p$theme$legend.position, "none")
  expect_error(plot_empirical_hazard(d, legend_position = "middle"))
})

test_that("colours, line types and text sizes are set by the arguments", {
  d <- make_fit_data()
  control <- unname(.plot_colours(arms = levels(d$ARM))["Control"])
  p <- plot_empirical_hazard(d, conf_int = TRUE,
                             colours = c(Treatment = "black"),
                             linetypes = c("dotted", "solid"),
                             text_sizes = c(axis_text = 8, title = 15))
  expect_identical(group_values(p, "GeomLine", "colour"),
                   c(control, "black"))
  expect_identical(group_values(p, "GeomLine", "linetype"),
                   c("dotted", "solid"))
  expect_identical(group_values(p, "GeomRibbon", "fill"),
                   c(control, "black"))
  expect_identical(theme_size(p, "axis.text.x"), 8)
  expect_identical(theme_size(p, "plot.title"), 15)
  default <- plot_empirical_hazard(d)
  expect_identical(group_values(default, "GeomLine", "linetype"),
                   c("solid", "solid"))
  expect_identical(theme_size(p, "legend.text"),
                   theme_size(default, "legend.text"))
})

test_that("invalid colours, line types and text sizes are rejected", {
  d <- make_fit_data()
  expect_error(plot_empirical_hazard(d, arms = "Control",
                                     colours = c(Treatment = "red")),
               "names of 'colours'")
  expect_error(plot_empirical_hazard(d, linetypes = "dashed"),
               "'linetypes' must have 2 values")
  expect_error(plot_empirical_hazard(d, text_sizes = c(title = 0)),
               "'text_sizes'")
})
