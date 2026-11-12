# Expected values: inst/validation/python/expected_fits.py (knots); the
# parameters and statistics of the fits are tested in the files of the
# table functions.

test_that("each arm and both arms are fitted with print-ready labels", {
  x <- get_fits("simple")
  expect_s3_class(x, "survfitsummary")
  expect_named(x$models, c("model_id", "fit_type", "arm", "family",
                           "distribution", "scale", "n_knots", "model", "n",
                           "events", "status", "message"))
  labels <- c("Exponential", "Weibull", "Weibull (PH)", "Rayleigh",
              "Lognormal", "Log-logistic", "Gompertz", "Gamma")
  expect_identical(x$models$model, rep(labels, 3))
  expect_identical(x$models$arm,
                   rep(c("Control", "Treatment", "Both arms"), each = 8))
  expect_identical(x$models$fit_type,
                   rep(c("independent", "independent", "dependent"),
                       each = 8))
  expect_identical(x$models$status, rep("ok", 24))
  expect_true(all(is.na(x$models$message)))
  expect_identical(x$models$n, rep(c(300L, 300L, 600L), each = 8))
  expect_identical(x$models$events, rep(c(271L, 240L, 511L), each = 8))
  expect_length(x$fits, 24)
  expect_true(all(vapply(x$fits, inherits, logical(1), "flexsurvreg")))
  expect_identical(x$arms, c("Control", "Treatment"))
  expect_identical(x$time_unit, "months")
  expect_identical(nrow(x$knots), 0L)
})

test_that("spline knots are at quantiles of the log event times", {
  x <- get_fits("hazard")
  expect_identical(x$models$model,
                   rep(c("Hazards, 1 knot", "Hazards, 2 knots",
                         "Hazards, 3 knots"), 3))
  expect_identical(x$models$n_knots, rep(1:3, 3))
  expect_identical(x$models$status, rep("ok", 9))
  e <- lookup(x$knots, read_expected("knots"),
              c("fit_type", "arm", "n_knots", "knot"))
  expect_identical(nrow(x$knots), 18L)
  expect_equal(x$knots$quantile, e$quantile, tolerance = 1e-10)
  expect_equal(x$knots$log_time, e$log_time, tolerance = 1e-10)
  # The knots passed to flexsurv are the boundary and internal knots
  d <- make_fit_data()
  for (i in seq_along(x$fits)) {
    m <- x$models[i, ]
    d_fit <- if (m$fit_type == "dependent") d else d[d$ARM == m$arm, ]
    log_events <- log(d_fit$SURVTIME[d_fit$EVENT == 1L])
    internal <- x$knots$log_time[x$knots$arm == m$arm &
                                   x$knots$n_knots == m$n_knots]
    expect_equal(unname(x$fits[[i]]$knots),
                 c(min(log_events), internal, max(log_events)),
                 tolerance = 1e-12)
  }
  expect_equal(x$boundary_knots$lower,
               log(c(min(d$SURVTIME[d$EVENT == 1L & d$ARM == "Control"]),
                     min(d$SURVTIME[d$EVENT == 1L & d$ARM == "Treatment"]),
                     min(d$SURVTIME[d$EVENT == 1L]))),
               tolerance = 1e-12)
})

test_that("knots given by the user are used", {
  d <- make_fit_data()
  x <- fit_survival_models(d, distributions = NULL, spline_scales = "odds",
                           knots = list(10, c(8, 16)),
                           fit_type = "independent")
  expect_identical(x$models$model, rep(c("Odds, 1 knot", "Odds, 2 knots"),
                                       2))
  expect_equal(x$knots$log_time, log(rep(c(10, 8, 16), 2)),
               tolerance = 1e-12)
  control_events <- d$SURVTIME[d$EVENT == 1L & d$ARM == "Control"]
  expect_equal(x$knots$quantile[1:3],
               c(mean(control_events <= 10), mean(control_events <= 8),
                 mean(control_events <= 16)),
               tolerance = 1e-12)
  expect_equal(unname(x$fits[[2]]$knots[2:3]), log(c(8, 16)),
               tolerance = 1e-12)
})

test_that("warnings and errors are recorded and the other models fitted", {
  d <- make_fit_data()
  # Too few iterations: the fit is kept with a warning
  x <- fit_survival_models(d, distributions = c("exp", "gengamma"),
                           spline_scales = NULL, fit_type = "independent",
                           control = list(maxit = 1))
  gg <- x$models$model == "Generalized Gamma"
  expect_identical(x$models$status[gg], rep("warning", 2))
  expect_match(x$models$message[gg], "converge", all = TRUE)
  expect_false(any(vapply(x$fits, is.null, logical(1))))
  # A single event in one arm: its spline model fails, the others are fitted
  d$EVENT[d$ARM == "Treatment"] <- 0L
  d$EVENT[d$ARM == "Treatment"][1] <- 1L
  x <- fit_survival_models(d, distributions = "exp",
                           spline_scales = "hazard", spline_knots = 1,
                           fit_type = "independent")
  expect_identical(x$models$status, c("ok", "ok", "ok", "error"))
  expect_match(x$models$message[4], "two distinct event times")
  expect_null(x$fits[[4]])
  expect_identical(unique(x$knots$arm), "Control")
  # An error of flexsurv is returned instead of stopping
  res <- .fit_one_model(d, dependent = FALSE, distribution = "unknown",
                        conf_level = 0.95, control = list())
  expect_identical(res$status, "error")
  expect_null(res$fit)
  expect_true(nzchar(res$message))
})

test_that("a single arm gives independent fits only", {
  d <- make_fit_data()
  d <- d[d$ARM == "Control", ]
  expect_warning(
    x <- fit_survival_models(d, distributions = "exp", spline_scales = NULL),
    "two or more treatment arms"
  )
  expect_identical(x$models$fit_type, "independent")
  expect_identical(x$arms, "Control")
})

test_that("invalid inputs are rejected", {
  d <- make_fit_data()
  expect_error(fit_survival_models(make_test_raw()), "processing_dataset")
  expect_error(fit_survival_models(d, distributions = "weibul"),
               "'distributions'")
  expect_error(fit_survival_models(d, spline_scales = "logit"),
               "'spline_scales'")
  expect_error(fit_survival_models(d, distributions = NULL,
                                   spline_scales = NULL), "No model")
  expect_error(fit_survival_models(d, spline_knots = 4), "'spline_knots'")
  expect_error(fit_survival_models(d, knots = list(c(10, 5))), "'knots'")
  expect_error(fit_survival_models(d, knots = list(10, 12)),
               "different numbers")
  expect_error(fit_survival_models(d, distributions = NULL,
                                   knots = list(500)), "strictly between")
  expect_error(fit_survival_models(d, fit_type = "joint"), "'fit_type'")
  expect_error(fit_survival_models(d, time_unit = "hours"))
  expect_error(fit_survival_models(d, conf_level = 1), "'conf_level'")
  expect_error(fit_survival_models(d, control = 1), "'control'")
})
