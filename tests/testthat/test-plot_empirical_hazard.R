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
