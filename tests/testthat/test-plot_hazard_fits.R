# The fitted hazards are those of predict_survival(type = "hazard") and the
# empirical hazards those of estimate_empirical_hazard() (both tested in
# their own files).

test_that("one arm with several models from time_horizon / n_points on", {
  x <- get_fits("simple")
  d <- make_fit_data()
  models <- c("Weibull", "Gompertz")
  p <- plot_hazard_fits(x, models = models, arms = "Treatment",
                        time_horizon = 100)
  expect_s3_class(p, "ggplot")
  expect_length(layer_index(p, "GeomLine"), 2)
  observed <- layer_input(p, "GeomLine", 1)
  hz <- estimate_empirical_hazard(d[d$ARM == "Treatment", ])
  expect_equal(observed$value, hz$hazard)
  expect_identical(as.character(unique(observed$curve)),
                   "Smoothed empirical hazard (Treatment)")
  fitted <- layer_input(p, "GeomLine", 2)
  times <- seq(0, 100, length.out = 201)[-1]
  pred <- predict_survival(x, times, type = "hazard")
  pred <- pred[pred$fit_type == "independent" & pred$arm == "Treatment" &
                 pred$model %in% models, ]
  expect_equal(fitted$est, pred$est)
  expect_true(min(fitted$time) > 0)
  y_max <- 1.2 * max(c(observed$value, fitted$est[fitted$time >= 5]))
  expect_equal(panel_ranges(p)$y, c(0, y_max))
  expect_equal(panel_ranges(p)$x, c(0, 100))
  expect_equal(panel_ranges(plot_hazard_fits(x, models = models,
                                             arms = "Treatment",
                                             y_max = 0.1))$y,
               c(0, 0.1))
})

test_that("one model with both arms and confidence bands", {
  x <- get_fits("simple")
  p <- plot_hazard_fits(x, models = "Log-logistic", conf_int = TRUE,
                        n_sim = 200)
  expect_length(layer_index(p, "GeomRibbon"), 2)
  observed <- layer_input(p, "GeomLine", 1)
  expect_identical(levels(observed$arm), c("Control", "Treatment"))
  band <- layer_input(p, "GeomRibbon", 2)
  expect_true(all(band$lcl < band$est & band$est < band$ucl))
  p <- plot_hazard_fits(x, models = c("Weibull", "Gamma"), arms = "Control",
                        conf_int = TRUE)
  expect_length(layer_index(p, "GeomRibbon"), 1)
})

test_that("invalid inputs are rejected", {
  x <- get_fits("simple")
  expect_error(plot_hazard_fits(list()), "fit_survival_models")
  expect_error(plot_hazard_fits(x), "either one arm")
  expect_error(plot_hazard_fits(x, models = "Weibull", y_max = -1),
               "'y_max'")
  expect_error(plot_hazard_fits(x, models = "Hazards, 1 knot"), "'models'")
})

test_that("the legend is placed with legend_position", {
  x <- get_fits("simple")
  p <- plot_hazard_fits(x, models = "Weibull", arms = "Control")
  expect_identical(p$theme$legend.position, "bottom")
  p <- plot_hazard_fits(x, models = "Weibull", arms = "Control",
                        legend_position = "top")
  expect_identical(p$theme$legend.position, "top")
  expect_error(plot_hazard_fits(x, legend_position = "middle"))
})
