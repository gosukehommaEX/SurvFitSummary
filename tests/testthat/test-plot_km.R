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
