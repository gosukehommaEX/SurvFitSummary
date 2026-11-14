# Kaplan-Meier values: inst/validation/python/expected_fits.py (column km of
# expected/predictions.csv at 6 to 60 months).

test_that("each curve starts at (0, 1) and ends at the last time of its arm", {
  d <- make_fit_data()
  p <- plot_km(d, risk_table = FALSE)
  expect_s3_class(p, "ggplot")
  km <- layer_input(p, "GeomStep")
  expected <- read_expected("predictions")
  for (a in c("Control", "Treatment")) {
    k <- km[km$arm == a, ]
    expect_identical(c(k$time[1], k$surv[1]), c(0, 1))
    expect_equal(max(k$time), max(d$SURVTIME[d$ARM == a]))
    e <- expected[expected$fit_type == "independent" & expected$arm == a &
                    expected$model == "Exponential" & expected$time > 0 &
                    !is.na(expected$km), ]
    step <- k$surv[findInterval(e$time, k$time)]
    expect_equal(step, e$km, tolerance = 1e-10)
  }
  r <- panel_ranges(p)
  expect_equal(r$x, c(0, max(d$SURVTIME)))
  expect_equal(r$y, c(0, 1))
  expect_length(layer_index(p, "GeomRibbon"), 0)
})

test_that("the numbers at risk are given at the breaks of the time axis", {
  d <- make_fit_data()
  p <- plot_km(d)
  expect_s3_class(p, "patchwork")
  table <- p[[2]]$data
  expect_identical(sort(unique(table$time)), seq(0, 60, by = 12))
  expected_n <- mapply(function(a, t) sum(d$SURVTIME[d$ARM == a] >= t),
                       as.character(table$arm), table$time)
  expect_identical(table$n, unname(expected_n))
  # A number is shown only up to the last observed time of the arm
  table <- plot_km(d, break_by = 6)[[2]]$data
  expect_identical(max(table$time[table$arm == "Control"]), 60)
  expect_identical(max(table$time[table$arm == "Treatment"]), 66)
})

test_that("with one arm, the row of numbers at risk is labelled N at risk", {
  d <- make_fit_data()
  title_drawn <- function(p) {
    g <- ggplot2::ggplotGrob(p)
    !inherits(g$grobs[[which(g$layout$name == "title")]], "zeroGrob")
  }
  table <- plot_km(d, arms = "Treatment")[[2]]
  expect_identical(ggplot2::layer_scales(table)$y$get_labels(), "N at risk")
  expect_false(title_drawn(table))
  expect_identical(unique(as.character(table$data$arm)), "Treatment")
  # With two arms, one row per arm below the title "N at risk"
  table <- plot_km(d)[[2]]
  expect_identical(ggplot2::layer_scales(table)$y$get_labels(),
                   c("Treatment", "Control"))
  expect_true(title_drawn(table))
})

test_that("the names of the arms are clear of the numbers at time 0", {
  d <- make_fit_data()
  table <- plot_km(d)[[2]]
  first <- table$data$n[table$data$time == 0]
  # The numbers are centred at time 0: half the width of the widest one, in
  # points (a digit is about 0.6 times the font size of 0.8 * 14 points)
  half_width <- 0.5 * max(nchar(first)) * 0.6 * 0.8 * 14
  expect_gt(as.numeric(table$theme$axis.text.y$margin)[2], half_width)
})

test_that("one arm with its confidence band", {
  d <- make_fit_data()
  p <- plot_km(d, arms = "Treatment", conf_int = TRUE, risk_table = FALSE)
  expect_identical(unique(as.character(layer_input(p, "GeomStep")$arm)),
                   "Treatment")
  band <- layer_input(p, "GeomRibbon")
  expect_true(all(0 <= band$lcl & band$lcl <= band$ucl & band$ucl <= 1))
  expect_equal(panel_ranges(p)$x,
               c(0, max(d$SURVTIME[d$ARM == "Treatment"])))
  # The arm keeps the colour it has in the plot of both arms
  step <- ggplot2::layer_data(p, layer_index(p, "GeomStep"))
  expect_identical(unique(step$colour),
                   unname(.plot_colours(arms = levels(d$ARM))["Treatment"]))
})

test_that("invalid inputs are rejected", {
  d <- make_fit_data()
  expect_error(plot_km(make_test_raw()), "processing_dataset")
  expect_error(plot_km(d, arms = "Placebo"), "'arms'")
  expect_error(plot_km(d, conf_level = 2), "'conf_level'")
  expect_error(plot_km(d, break_by = -1), "'break_by'")
  expect_error(plot_km(d, time_unit = "hours"))
})

test_that("the legend is placed with legend_position", {
  d <- make_fit_data()
  p <- plot_km(d, risk_table = FALSE)
  expect_identical(p$theme$legend.position, "bottom")
  p <- plot_km(d, legend_position = "right")
  expect_identical(p[[1]]$theme$legend.position, "right")
  expect_error(plot_km(d, legend_position = "middle"))
})

test_that("colours, line types and text sizes are set by the arguments", {
  d <- make_fit_data()
  control <- unname(.plot_colours(arms = levels(d$ARM))["Control"])
  default <- plot_km(d)
  expect_identical(group_values(default[[1]], "GeomStep", "linetype"),
                   c("solid", "solid"))
  p <- plot_km(d, conf_int = TRUE, colours = c(Treatment = "red"),
               linetypes = c(Control = "dashed"),
               text_sizes = c(axis_title = 20, axis_text = 9, legend = 8,
                              risk_table = 10))
  main <- p[[1]]
  expect_identical(group_values(main, "GeomStep", "colour"),
                   c(control, "red"))
  expect_identical(group_values(main, "GeomStep", "linetype"),
                   c("dashed", "solid"))
  expect_identical(group_values(main, "GeomRibbon", "fill"),
                   c(control, "red"))
  expect_identical(theme_size(main, "axis.title.x"), 20)
  expect_identical(theme_size(main, "axis.text.y"), 9)
  expect_identical(theme_size(main, "legend.text"), 8)
  expect_identical(theme_size(main, "plot.title"),
                   theme_size(default[[1]], "plot.title"))
  # The numbers at risk in the colours of the arms and in the given size
  table <- p[[2]]
  numbers <- ggplot2::layer_data(table)
  expect_identical(unique(numbers$colour[table$data$arm == "Treatment"]),
                   "red")
  expect_equal(unique(numbers$size), 10 / ggplot2::.pt)
  expect_identical(table$theme$axis.text.y$size, 10)
  expect_equal(table$theme$plot.title$size, 11.25)
  expect_equal(default[[2]]$theme$plot.title$size, 0.9 * 14)
  # An unnamed vector gives the values in the order of the arms
  p <- plot_km(d, risk_table = FALSE, colours = c("black", "grey50"),
               linetypes = c("solid", "44"))
  expect_identical(group_values(p, "GeomStep", "colour"),
                   c("black", "grey50"))
  expect_identical(group_values(p, "GeomStep", "linetype"),
                   c("solid", "44"))
})

test_that("invalid colours, line types and text sizes are rejected", {
  d <- make_fit_data()
  expect_error(plot_km(d, colours = c(Placebo = "red")),
               "names of 'colours'")
  expect_error(plot_km(d, colours = "red"), "'colours' must have 2 values")
  expect_error(plot_km(d, colours = c("red", "no colour")),
               "invalid colours: no colour")
  expect_error(plot_km(d, colours = c(Control = NA_character_)),
               "'colours' must be NULL")
  expect_error(plot_km(d, linetypes = c(Control = "wavy")),
               "invalid line types: wavy")
  expect_error(plot_km(d, linetypes = 1:2), "'linetypes' must be NULL")
  expect_error(plot_km(d, text_sizes = c(axis = 12)), "'text_sizes'")
  expect_error(plot_km(d, text_sizes = c(legend = -1)), "'text_sizes'")
  expect_error(plot_km(d, text_sizes = 12), "'text_sizes'")
})
