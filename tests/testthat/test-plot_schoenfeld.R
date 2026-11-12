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
