# Expected values: inst/validation/python/expected_fits.py (independent
# maximum likelihood fits with numpy and scipy). The transformed scale is
# compared with the reference; the natural scale is checked against the
# transformed scale and against flexsurv's own back-transformation
# (flexsurvreg()$res).

check_parameters <- function(set) {
  x <- get_fits(set)
  out <- parameter_table(x)
  expected <- expected_for("parameters", x)
  expect_identical(nrow(out), nrow(expected))
  keys <- c("fit_type", "arm", "model", "parameter", "scale")
  tr <- out[out$scale == "transformed", ]
  e <- lookup(tr, expected, keys)
  expect_lte(scaled_diff(tr$est, e$est, TOL_EST * e$se), 1)
  expect_lte(scaled_diff(tr$lcl, e$lcl, TOL_CI * e$se), 1)
  expect_lte(scaled_diff(tr$ucl, e$ucl, TOL_CI * e$se), 1)
  expect_lte(scaled_diff(tr$se, e$se, TOL_SE * e$se), 1)
  # Natural scale: exponential of the log-transformed parameters
  na <- out[out$scale == "natural", ]
  id <- function(d) paste(d$fit_type, d$arm, d$model, d$parameter)
  expect_identical(id(na), id(tr))
  e_na <- lookup(na, expected, keys)
  is_log <- abs(e_na$est - e$est) > 1e-12
  expect_equal(na$est, ifelse(is_log, exp(tr$est), tr$est),
               tolerance = 1e-12)
  expect_equal(na$lcl, ifelse(is_log, exp(tr$lcl), tr$lcl),
               tolerance = 1e-12)
  expect_equal(na$ucl, ifelse(is_log, exp(tr$ucl), tr$ucl),
               tolerance = 1e-12)
  expect_equal(na$se, ifelse(is_log, na$est * tr$se, tr$se),
               tolerance = 1e-12)
  # The same as flexsurv's estimates on the natural scale
  for (i in seq_along(x$fits)) {
    rows <- na[na$fit_type == x$models$fit_type[i] &
                 na$arm == x$models$arm[i] &
                 na$model == x$models$model[i], ]
    res <- x$fits[[i]]$res[rows$parameter, 1:4, drop = FALSE]
    expect_equal(unname(as.matrix(rows[c("est", "lcl", "ucl", "se")])),
                 unname(res), tolerance = 1e-6)
  }
}

test_that("standard distributions match the independent fits", {
  check_parameters("simple")
})

test_that("generalized gamma distributions match the independent fits", {
  check_parameters("gengamma")
})

test_that("generalized F (Prentice) distributions match the independent fits", {
  check_parameters("genf")
})

test_that("original generalized F distributions match the independent fits", {
  check_parameters("genf_orig")
})

test_that("hazard spline models match the independent fits", {
  check_parameters("hazard")
})

test_that("odds spline models match the independent fits", {
  check_parameters("odds")
})

test_that("normal spline models match the independent fits", {
  check_parameters("normal")
})

test_that("the scale argument and the Rayleigh shape", {
  x <- get_fits("simple")
  out <- parameter_table(x, scale = "natural")
  expect_identical(unique(out$scale), "natural")
  expect_named(out, c("fit_type", "arm", "model", "parameter", "scale",
                      "est", "lcl", "ucl", "se"))
  expect_identical(out$parameter[out$model == "Rayleigh"],
                   c("scale", "scale", "scale", "ARMTreatment"))
  expect_identical(
    out$parameter[out$model == "Weibull" & out$arm == "Both arms"],
    c("shape", "scale", "ARMTreatment")
  )
  expect_error(parameter_table(x, scale = "log"))
  expect_error(parameter_table(list()), "fit_survival_models")
})
