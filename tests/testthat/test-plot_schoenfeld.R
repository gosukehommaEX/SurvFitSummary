# The log hazard ratio and p-value are those of test_proportional_hazards()
# (compared with inst/validation/python/expected_km_cox.py in its own test
# file). The regression line is checked by another computation: lm.fit() on
# the same natural spline basis, and the interval from (X'X)^-1 instead of
# the QR factor.

test_that("residuals, log hazard ratio and regression line", {
  d <- make_fit_data()
  p <- plot_schoenfeld(d)
  expect_s3_class(p, "ggplot")
  fit <- survival::coxph(survival::Surv(SURVTIME, EVENT) ~ ARM, data = d)
  zph <- survival::cox.zph(fit, transform = "identity")
  points <- layer_input(p, "GeomPoint")
  expect_equal(points$y, unname(zph$y[, 1]))
  expect_equal(points$time, unname(zph$time))
  reference <- layer_input(p, "GeomHline")
  ph <- test_proportional_hazards(d)
  expect_equal(reference$log_hr, ph$log_hr)
  expect_match(reference$label, sprintf("p = %.3f", ph$p_value),
               fixed = TRUE)
  smooth <- layer_input(p, "GeomLine", 1)
  grid <- smooth$time
  basis <- splines::ns(c(grid, points$time), df = 4, intercept = TRUE)
  pmat <- basis[seq_along(grid), ]
  xmat <- basis[-seq_along(grid), ]
  coefs <- stats::lm.fit(xmat, points$y)$coefficients
  expect_equal(smooth$fit, drop(pmat %*% coefs))
  half <- stats::qnorm(0.975) *
    sqrt(zph$var[1, 1] * rowSums((pmat %*% solve(crossprod(xmat))) * pmat))
  expect_equal(smooth$ucl - smooth$fit, half)
  expect_equal(smooth$fit - smooth$lcl, half)
})

test_that("panel title with the test and breaks of the time axis", {
  d <- make_fit_data()
  p <- plot_schoenfeld(d)
  ph <- test_proportional_hazards(d)
  reference <- layer_input(p, "GeomHline")
  expect_identical(
    as.character(reference$panel),
    sprintf(paste("Treatment vs Control (log HR = %.3f;",
                  "Grambsch-Therneau test p = %.3f)"),
            ph$log_hr, ph$p_value)
  )
  expect_identical(levels(layer_input(p, "GeomPoint")$panel),
                   as.character(reference$panel))
  expect_length(layer_index(p, "GeomText"), 0)
  breaks <- ggplot2::layer_scales(p)$x$get_breaks()
  expect_equal(breaks[!is.na(breaks)], seq(0, 60, by = 12))
  # Legends of the points and of the lines one above the other
  expect_identical(p$theme$legend.box, "vertical")
})

test_that("one panel per comparison with the control arm", {
  d <- make_fit_data()
  arm <- as.character(d$ARM)
  treated <- which(arm == "Treatment")
  arm[treated[seq(1, length(treated), by = 2)]] <- "Treatment B"
  d$ARM <- factor(arm, levels = c("Control", "Treatment", "Treatment B"))
  p <- plot_schoenfeld(d)
  reference <- layer_input(p, "GeomHline")
  expect_identical(as.character(reference$comparison),
                   c("Treatment vs Control", "Treatment B vs Control"))
  expect_equal(reference$log_hr, test_proportional_hazards(d)$log_hr)
})

test_that("invalid inputs are rejected", {
  d <- make_fit_data()
  expect_error(plot_schoenfeld(d[d$ARM == "Control", ]),
               "At least 2 treatment arms")
  expect_error(plot_schoenfeld(d, conf_level = 1), "'conf_level'")
})

test_that("the legend is placed with legend_position", {
  d <- make_fit_data()
  expect_identical(plot_schoenfeld(d)$theme$legend.position, "bottom")
  p <- plot_schoenfeld(d, legend_position = "right")
  expect_identical(p$theme$legend.position, "right")
  expect_error(plot_schoenfeld(d, legend_position = "middle"))
})

test_that("colours, line types and text sizes are set by the arguments", {
  d <- make_fit_data()
  p <- plot_schoenfeld(d, colours = c(residuals = "grey50", log_hr = "red"),
                       linetypes = c(regression = "dashed",
                                     confidence_interval = "dotdash"),
                       text_sizes = c(title = 10, axis_text = 9))
  expect_identical(group_values(p, "GeomPoint", "colour"), "grey50")
  expect_identical(group_values(p, "GeomHline", "colour"), "red")
  expect_identical(group_values(p, "GeomHline", "linetype"), "solid")
  expect_identical(group_values(p, "GeomLine", "colour", 1), "black")
  expect_identical(group_values(p, "GeomLine", "linetype", 1), "dashed")
  expect_identical(group_values(p, "GeomLine", "linetype", 2), "dotdash")
  expect_identical(theme_size(p, "strip.text.x"), 10)
  expect_identical(theme_size(p, "axis.text.x"), 9)
  # An unnamed vector gives the colours in the order of the elements
  p <- plot_schoenfeld(d, colours = c("black", "blue", "grey40", "grey60"))
  expect_identical(group_values(p, "GeomPoint", "colour"), "black")
  expect_identical(group_values(p, "GeomHline", "colour"), "blue")
  expect_identical(group_values(p, "GeomLine", "colour", 1), "grey40")
  expect_identical(group_values(p, "GeomLine", "colour", 3), "grey60")
})

test_that("invalid colours, line types and text sizes are rejected", {
  d <- make_fit_data()
  expect_error(plot_schoenfeld(d, colours = c(points = "red")),
               "names of 'colours'")
  expect_error(plot_schoenfeld(d, linetypes = c(residuals = "dashed")),
               "names of 'linetypes'")
  expect_error(plot_schoenfeld(d, colours = c("red", "blue")),
               "'colours' must have 4 values")
})
