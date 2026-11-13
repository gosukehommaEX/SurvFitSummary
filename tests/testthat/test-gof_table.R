# Expected values: inst/validation/python/expected_fits.py (maximum of the
# log-likelihood; AIC and BIC with n = number of patients).

check_gof <- function(set) {
  x <- get_fits(set)
  out <- gof_table(x)
  expected <- expected_for("gof", x)
  expect_identical(nrow(out), nrow(expected))
  e <- lookup(out, expected, c("fit_type", "arm", "model"))
  expect_identical(out$n, as.integer(e$n))
  expect_identical(out$events, as.integer(e$events))
  expect_identical(out$npars, as.integer(e$npars))
  expect_lte(scaled_diff(out$loglik, e$loglik, TOL_LOGLIK), 1)
  if (set %in% c("hazard", "odds", "normal")) {
    # The spline fits reach the maximum (see TOL_LOGLIK_SPLINE)
    expect_lte(scaled_diff(out$loglik, e$loglik, TOL_LOGLIK_SPLINE), 1)
  }
  expect_lte(scaled_diff(out$AIC, e$AIC, 2 * TOL_LOGLIK), 1)
  expect_lte(scaled_diff(out$BIC, e$BIC, 2 * TOL_LOGLIK), 1)
  # Sorted by AIC within each fit, the fits in the order of the models
  expect_identical(unique(out$arm), c("Control", "Treatment", "Both arms"))
  for (a in unique(out$arm)) {
    expect_false(is.unsorted(out$AIC[out$arm == a]))
  }
}

test_that("standard distributions match the independent fits", {
  check_gof("simple")
})

test_that("generalized gamma distributions match the independent fits", {
  check_gof("gengamma")
})

test_that("generalized F (Prentice) distributions match the independent fits", {
  check_gof("genf")
})

test_that("original generalized F distributions match the independent fits", {
  check_gof("genf_orig")
})

test_that("hazard spline models match the independent fits", {
  check_gof("hazard")
})

test_that("odds spline models match the independent fits", {
  check_gof("odds")
})

test_that("normal spline models match the independent fits", {
  check_gof("normal")
})

test_that("AIC and BIC agree with flexsurv and Rayleigh counts one shape", {
  x <- get_fits("simple")
  out <- gof_table(x)
  expect_named(out, c("fit_type", "arm", "model", "n", "events", "npars",
                      "loglik", "AIC", "BIC"))
  rayleigh <- out$model == "Rayleigh"
  expect_identical(out$npars[rayleigh], c(1L, 1L, 2L))
  for (i in which(x$models$model != "Rayleigh")) {
    row <- out$fit_type == x$models$fit_type[i] &
      out$arm == x$models$arm[i] & out$model == x$models$model[i]
    expect_equal(out$AIC[row], x$fits[[i]]$AIC, tolerance = 1e-10)
    expect_equal(out$BIC[row], stats::BIC(x$fits[[i]]), tolerance = 1e-10)
  }
  expect_error(gof_table(list()), "fit_survival_models")
})
