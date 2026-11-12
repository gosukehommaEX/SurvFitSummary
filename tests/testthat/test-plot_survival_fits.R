# The fitted curves are those of predict_survival() and the AIC those of
# gof_table() (both compared with inst/validation/python/expected_fits.py
# in their own test files).

test_that("one arm with several models: black Kaplan-Meier and AIC labels", {
  x <- get_fits("simple")
  d <- make_fit_data()
  models <- c("Exponential", "Weibull", "Log-logistic")
  p <- plot_survival_fits(x, models = models, arms = "Control",
                          time_horizon = 120, risk_table = FALSE)
  expect_s3_class(p, "ggplot")
  lines <- layer_input(p, "GeomLine")
  gof <- gof_table(x)
  gof <- gof[gof$arm == "Control", ]
  labels <- paste0(models, " (AIC = ",
                   formatC(gof$AIC[match(models, gof$model)], format = "f",
                           digits = 1), ")")
  expect_identical(levels(lines$curve), c("Kaplan-Meier (Control)", labels))
  times <- seq(0, 120, length.out = 200)
  pred <- predict_survival(x, times)
  pred <- pred[pred$fit_type == "independent" & pred$arm == "Control" &
                 pred$model %in% models, ]
  expect_equal(lines$est, pred$est)
  expect_identical(as.character(lines$curve),
                   labels[match(pred$model, models)])
  km <- layer_input(p, "GeomStep")
  expect_identical(as.character(unique(km$curve)), "Kaplan-Meier (Control)")
  expect_equal(max(km$time), max(d$SURVTIME[d$ARM == "Control"]))
  expect_equal(panel_ranges(p)$x, c(0, 120))
  expect_equal(panel_ranges(p)$y, c(0, 1))
  step <- ggplot2::layer_data(p, layer_index(p, "GeomStep"))
  expect_identical(unique(step$colour), "black")
})

test_that("one model with both arms from the dependent fit", {
  x <- get_fits("simple")
  d <- make_fit_data()
  p <- plot_survival_fits(x, models = "Weibull", fit_type = "dependent",
                          risk_table = FALSE)
  lines <- layer_input(p, "GeomLine")
  expect_identical(levels(lines$arm), c("Control", "Treatment"))
  expect_identical(levels(lines$type), c("Kaplan-Meier", "Weibull"))
  pred <- predict_survival(x, seq(0, max(d$SURVTIME), length.out = 200))
  pred <- pred[pred$fit_type == "dependent" & pred$model == "Weibull", ]
  expect_equal(lines$est, pred$est)
  km <- layer_input(p, "GeomStep")
  for (a in c("Control", "Treatment")) {
    expect_equal(max(km$time[km$arm == a]), max(d$SURVTIME[d$ARM == a]))
  }
  p <- plot_survival_fits(x, models = "Weibull", fit_type = "dependent")
  expect_s3_class(p, "patchwork")
  expect_identical(sort(unique(as.character(p[[2]]$data$arm))),
                   c("Control", "Treatment"))
})

test_that("arms keep their colours; numbers at risk black with one arm", {
  d <- make_fit_data()
  arm <- as.character(d$ARM)
  treated <- which(arm == "Treatment")
  arm[treated[seq(1, length(treated), by = 2)]] <- "Treatment B"
  d$ARM <- factor(arm, levels = c("Control", "Treatment", "Treatment B"))
  x <- fit_survival_models(d, distributions = "exp", spline_scales = NULL,
                           fit_type = "independent")
  shown <- c("Treatment", "Treatment B")
  colours <- unname(.plot_colours(arms = x$arms)[shown])
  p <- plot_survival_fits(x, models = "Exponential", arms = shown)
  main <- p[[1]]
  step <- ggplot2::layer_data(main, layer_index(main, "GeomStep"))
  expect_setequal(unique(step$colour), colours)
  expect_setequal(unique(ggplot2::layer_data(p[[2]])$colour), colours)
  p <- plot_survival_fits(x, models = "Exponential", arms = "Treatment B")
  expect_identical(unique(ggplot2::layer_data(p[[2]])$colour), "black")
})

test_that("confidence bands: Kaplan-Meier always, a model only alone", {
  x <- get_fits("simple")
  p <- plot_survival_fits(x, models = c("Exponential", "Weibull"),
                          arms = "Control", conf_int = TRUE,
                          risk_table = FALSE)
  expect_length(layer_index(p, "GeomRibbon"), 1)
  p <- plot_survival_fits(x, models = "Weibull", arms = "Control",
                          conf_int = TRUE, n_sim = 200, risk_table = FALSE)
  expect_length(layer_index(p, "GeomRibbon"), 2)
  band <- layer_input(p, "GeomRibbon", 2)
  inside <- band$time > 0
  expect_true(all(band$lcl[inside] < band$est[inside] &
                    band$est[inside] < band$ucl[inside]))
})

test_that("invalid inputs are rejected", {
  x <- get_fits("simple")
  expect_error(plot_survival_fits(list()), "fit_survival_models")
  expect_error(plot_survival_fits(x), "either one arm")
  expect_error(plot_survival_fits(x, models = "Weibul", arms = "Control"),
               "'models'")
  expect_error(plot_survival_fits(x, models = "Weibull", arms = "Placebo"),
               "'arms'")
  expect_error(plot_survival_fits(x, models = "Weibull", time_horizon = -1),
               "'time_horizon'")
  expect_error(plot_survival_fits(x, models = "Weibull", n_points = 1),
               "'n_points'")
  expect_error(plot_survival_fits(x, models = "Weibull", conf_int = NA),
               "'conf_int'")
  expect_error(plot_survival_fits(x, models = "Weibull",
                                  fit_type = "joint"))
})
