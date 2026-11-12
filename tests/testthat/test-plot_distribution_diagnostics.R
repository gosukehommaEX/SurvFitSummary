# Each plot is checked with data simulated from the distribution it is
# meant for: the regression line has the slope given by theory (Weibull:
# shape; log-logistic: shape; lognormal: 1 / sdlog; Gompertz: shape).

test_that("each plot is a straight line with the theoretical slope", {
  withr::local_seed(42)
  u <- stats::runif(3000)
  weibull <- single_arm(10 * (-log(u))^(1 / 1.5))
  llogis <- single_arm(10 * (u / (1 - u))^(1 / 2))
  lnorm <- single_arm(exp(2 + 0.8 * stats::qnorm(u)))
  gompertz <- single_arm(log(1 - 0.1 * log(u) / 0.05) / 0.1)
  expect_equal(drawn_slope(plot_distribution_diagnostics(
    weibull, which = "weibull")), 1.5, tolerance = 0.05)
  expect_equal(drawn_slope(plot_distribution_diagnostics(
    llogis, which = "loglogistic")), 2, tolerance = 0.05)
  expect_equal(drawn_slope(plot_distribution_diagnostics(
    lnorm, which = "lognormal")), 1 / 0.8, tolerance = 0.05)
  expect_equal(drawn_slope(plot_distribution_diagnostics(
    gompertz, which = "gompertz")), 0.1, tolerance = 0.2)
})

test_that("the points are transformations of the Kaplan-Meier estimate", {
  d <- make_fit_data()
  ctrl <- d[d$ARM == "Control", ]
  fit <- survival::survfit(survival::Surv(SURVTIME, EVENT) ~ 1, data = ctrl)
  keep <- fit$n.event > 0 & fit$surv > 0 & fit$surv < 1
  s <- fit$surv[keep]
  expected <- list(weibull = log(-log(s)), loglogistic = log((1 - s) / s),
                   lognormal = stats::qnorm(1 - s))
  for (w in names(expected)) {
    p <- plot_distribution_diagnostics(d, which = w, regression = FALSE)
    expect_s3_class(p, "ggplot")
    expect_length(layer_index(p, "GeomLine"), 0)
    points <- p$data[p$data$arm == "Control", ]
    expect_equal(points$x, log(fit$time[keep]))
    expect_equal(points$y, expected[[w]])
  }
})

test_that("the Gompertz plot has the breaks of the other time axes", {
  d <- make_fit_data()
  p <- plot_distribution_diagnostics(d, which = "gompertz")
  breaks <- ggplot2::layer_scales(p)$x$get_breaks()
  expect_equal(breaks[!is.na(breaks)], seq(0, 60, by = 12))
})

test_that("all four plots in the order of 'which'", {
  d <- make_fit_data()
  p <- plot_distribution_diagnostics(d)
  expect_s3_class(p, "patchwork")
  hz <- estimate_empirical_hazard(d)
  hz <- hz[hz$hazard > 0, ]
  expect_equal(p[[4]]$data$y, log(hz$hazard))
  expect_equal(p[[4]]$data$x, hz$time)
  p <- plot_distribution_diagnostics(d, which = c("lognormal", "weibull"))
  expect_equal(p[[2]]$data$y,
               plot_distribution_diagnostics(d, which = "weibull")$data$y)
  expect_error(plot_distribution_diagnostics(d, which = "gamma"))
  expect_error(plot_distribution_diagnostics(make_test_raw()),
               "processing_dataset")
})
