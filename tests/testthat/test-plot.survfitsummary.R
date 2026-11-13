# plot() calls plot_survival_fits() or plot_hazard_fits() (tested in their
# own files); these tests check the default arm and the arguments passed on.

test_that("plot() draws the survival or hazard fits of the control arm", {
  x <- get_fits("simple")
  p <- plot(x, risk_table = FALSE)
  expected <- plot_survival_fits(x, arms = "Control", risk_table = FALSE)
  expect_identical(levels(layer_input(p, "GeomLine")$curve),
                   levels(layer_input(expected, "GeomLine")$curve))
  expect_equal(layer_input(p, "GeomLine")$est,
               layer_input(expected, "GeomLine")$est)
  expect_s3_class(plot(x), "patchwork")
  p <- plot(x, type = "hazard")
  expected <- plot_hazard_fits(x, arms = "Control")
  expect_identical(levels(layer_input(p, "GeomLine", 2)$curve),
                   levels(layer_input(expected, "GeomLine", 2)$curve))
  expect_equal(layer_input(p, "GeomLine", 2)$est,
               layer_input(expected, "GeomLine", 2)$est)
})

test_that("the arms, the models and the other arguments are passed on", {
  x <- get_fits("simple")
  # A single model: all the arms
  p <- plot(x, models = "Weibull", risk_table = FALSE)
  expect_identical(levels(layer_input(p, "GeomLine")$arm),
                   c("Control", "Treatment"))
  p <- plot(x, type = "hazard", models = "Weibull", arms = "Treatment",
            colours = c(Weibull = "red"), legend_position = "right")
  expect_identical(as.character(unique(layer_input(p, "GeomLine")$curve)),
                   "Smoothed empirical hazard (Treatment)")
  expect_identical(group_values(p, "GeomLine", "colour", 2), "red")
  expect_identical(p$theme$legend.position, "right")
  # The dependent models for the control arm
  p <- plot(x, models = c("Exponential", "Weibull"), fit_type = "dependent",
            risk_table = FALSE)
  horizon <- max(x$data$SURVTIME[x$data$ARM == "Control"])
  pred <- predict_survival(x, seq(0, horizon, length.out = 200))
  pred <- pred[pred$fit_type == "dependent" & pred$arm == "Control" &
                 pred$model %in% c("Exponential", "Weibull"), ]
  expect_equal(layer_input(p, "GeomLine")$est, pred$est)
})

test_that("invalid inputs are rejected", {
  x <- get_fits("simple")
  expect_error(plot(x, type = "density"))
  expect_error(plot.survfitsummary(list()), "fit_survival_models")
  expect_error(plot(x, arms = "Placebo"), "'arms'")
})
