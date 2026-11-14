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
  expect_error(plot_survival_fits(x, models = "Weibull", arms = "Control",
                                  range_times = -1), "'range_times'")
  expect_error(plot_survival_fits(x, models = "Weibull", arms = "Control",
                                  range_times = NA), "'range_times'")
  expect_error(plot_survival_fits(x, models = "Weibull",
                                  fit_type = "joint"))
})

test_that("the legend is placed with legend_position", {
  x <- get_fits("simple")
  models <- c("Exponential", "Weibull")
  p <- plot_survival_fits(x, models = models, arms = "Control",
                          risk_table = FALSE)
  expect_identical(p$theme$legend.position, "bottom")
  p <- plot_survival_fits(x, models = models, arms = "Control",
                          legend_position = "right")
  expect_identical(p[[1]]$theme$legend.position, "right")
  expect_error(plot_survival_fits(x, legend_position = "middle"))
})

test_that("colours, line types and text sizes are set by the arguments", {
  x <- get_fits("simple")
  models <- c("Exponential", "Weibull")
  exponential <- unname(.plot_colours(
    models = x$models[match("Exponential", x$models$model), ]
  ))
  # One arm: by the Kaplan-Meier curve and the models
  p <- plot_survival_fits(x, models = models, arms = "Control",
                          colours = c("Kaplan-Meier" = "grey40",
                                      Weibull = "red"),
                          linetypes = c(Exponential = "dashed"),
                          text_sizes = c(legend = 9, risk_table = 10))
  main <- p[[1]]
  expect_identical(group_values(main, "GeomStep", "colour"), "grey40")
  expect_identical(group_values(main, "GeomStep", "linetype"), "solid")
  expect_identical(group_values(main, "GeomLine", "colour"),
                   c(exponential, "red"))
  expect_identical(group_values(main, "GeomLine", "linetype"),
                   c("dashed", "solid"))
  expect_identical(theme_size(main, "legend.text"), 9)
  # The numbers at risk in the colour of the Kaplan-Meier curve
  numbers <- ggplot2::layer_data(p[[2]])
  expect_identical(unique(numbers$colour), "grey40")
  expect_equal(unique(numbers$size), 10 / ggplot2::.pt)
  # An unnamed vector gives the values in the order of the legend
  p <- plot_survival_fits(x, models = models, arms = "Control",
                          risk_table = FALSE,
                          colours = c("black", "blue", "orange"))
  expect_identical(group_values(p, "GeomLine", "colour"),
                   c("blue", "orange"))
  # Several arms: colours by arm, line types by curve
  p <- plot_survival_fits(x, models = "Weibull", fit_type = "dependent",
                          risk_table = FALSE,
                          colours = c("black", "grey50"),
                          linetypes = c(Weibull = "dotted"))
  expect_identical(group_values(p, "GeomStep", "colour"),
                   c("black", "grey50"))
  expect_identical(group_values(p, "GeomStep", "linetype"),
                   c("solid", "solid"))
  expect_identical(group_values(p, "GeomLine", "colour"),
                   c("black", "grey50"))
  expect_identical(group_values(p, "GeomLine", "linetype"),
                   c("dotted", "dotted"))
})

test_that("invalid colours, line types and text sizes are rejected", {
  x <- get_fits("simple")
  models <- c("Exponential", "Weibull")
  expect_error(plot_survival_fits(x, models = models, arms = "Control",
                                  colours = c(Gamma = "red")),
               "names of 'colours'")
  expect_error(plot_survival_fits(x, models = models, arms = "Control",
                                  colours = c("red", "blue")),
               "'colours' must have 3 values")
  expect_error(plot_survival_fits(x, models = "Weibull",
                                  fit_type = "dependent",
                                  linetypes = c(Control = "dashed")),
               "names of 'linetypes'")
  expect_error(plot_survival_fits(x, models = "Weibull", arms = "Control",
                                  text_sizes = c(size = 3)),
               "'text_sizes'")
})

test_that("range_times marks the range of the survival of the models", {
  x <- get_fits("simple")
  models <- c("Exponential", "Weibull", "Log-logistic")
  p <- plot_survival_fits(x, models = models, arms = "Control",
                          time_horizon = 120, risk_table = FALSE,
                          range_times = c(150, 60, 120))
  ranges <- layer_input(p, "GeomSegment")
  # Sorted, and the times after the time horizon left out
  expect_identical(ranges$time, c(60, 120))
  e <- read_expected("predictions")
  e <- e[e$fit_type == "independent" & e$arm == "Control", ]
  in_plot <- e[e$model %in% models, ]
  low <- vapply(c(60, 120), function(t) min(in_plot$surv[in_plot$time == t]),
                numeric(1))
  high <- vapply(c(60, 120), function(t) max(in_plot$surv[in_plot$time == t]),
                 numeric(1))
  # Absolute differences (edition 3 of testthat has no scale argument)
  expect_lte(max(abs(ranges$low - low)), TOL_SURV)
  expect_lte(max(abs(ranges$high - high)), TOL_SURV)
  # Labels in percent, from the reference values (0.04 to 2.25 percent at
  # 60 months, 0.00 to 0.56 percent at 120 months)
  labels <- layer_input(p, "GeomText")
  expect_identical(labels$label, c("0-2%", "0-1%"))
  # Right of the line, and left of it at the time horizon
  expect_identical(labels$hjust, c(-0.15, 1.15))
  expect_identical(labels$high, ranges$high)
  # In the size of the axis labels
  text_layer <- p$layers[[layer_index(p, "GeomText")]]
  expect_equal(text_layer$aes_params$size, 0.8 * 14 / ggplot2::.pt)
  p <- plot_survival_fits(x, models = models, arms = "Control",
                          time_horizon = 120, risk_table = FALSE,
                          range_times = 60, text_sizes = c(axis_text = 9))
  text_layer <- p$layers[[layer_index(p, "GeomText")]]
  expect_equal(text_layer$aes_params$size, 9 / ggplot2::.pt)
  # One model: one value (Weibull at 12 months, 0.411 in the reference)
  p <- plot_survival_fits(x, models = "Weibull", arms = "Control",
                          time_horizon = 120, risk_table = FALSE,
                          range_times = 12)
  expect_identical(layer_input(p, "GeomText")$label, "41%")
  # No line without range_times, or with all the times after the horizon
  p <- plot_survival_fits(x, models = models, arms = "Control",
                          time_horizon = 120, risk_table = FALSE)
  expect_length(layer_index(p, "GeomSegment"), 0)
  p <- plot_survival_fits(x, models = models, arms = "Control",
                          time_horizon = 120, risk_table = FALSE,
                          range_times = 200)
  expect_length(layer_index(p, "GeomSegment"), 0)
  # With the numbers at risk, the lines are in the plot
  p <- plot_survival_fits(x, models = models, arms = "Control",
                          time_horizon = 120, range_times = 60)
  expect_identical(layer_input(p[[1]], "GeomSegment")$time, 60)
})
