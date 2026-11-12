# Expected values: inst/validation/python/expected_fits.py (median by root
# finding; mean and restricted mean by numerical integration, the mean also
# by closed forms where they exist). The mean is compared where the
# reference flags it as reliable (mean_check = 1: closed form and integral
# agree, or a light tail); time horizon 60.

check_summary <- function(set) {
  x <- get_fits(set)
  expected <- expected_for("summary", x)
  out <- survival_summary(x, time_horizon = 60)
  expect_identical(nrow(out), nrow(expected))
  e <- lookup(out, expected, c("fit_type", "arm", "model"))
  expect_lte(scaled_diff(out$median, e$median, TOL_TIME * e$median), 1)
  expect_lte(scaled_diff(out$rmst, e$rmst, TOL_TIME * e$rmst), 1)
  checked <- e$mean_check == 1
  if (any(checked)) {
    expect_lte(scaled_diff(out$mean[checked], e$mean[checked],
                           TOL_TIME * e$mean[checked]), 1)
  }
  expect_true(all(is.na(out$median_lcl)) && all(is.na(out$rmst_ucl)))
}

test_that("standard distributions match the independent values", {
  check_summary("simple")
})

test_that("generalized gamma distributions match the independent values", {
  check_summary("gengamma")
})

test_that("generalized F (Prentice) distributions match the independent values", {
  check_summary("genf")
})

test_that("original generalized F distributions match the independent values", {
  check_summary("genf_orig")
})

test_that("hazard spline models match the independent values", {
  check_summary("hazard")
})

test_that("odds spline models match the independent values", {
  check_summary("odds")
})

test_that("normal spline models match the independent values", {
  check_summary("normal")
})

test_that("layout, time horizon and confidence limits", {
  x <- get_fits("simple")
  out <- survival_summary(x)
  expect_named(out, c("fit_type", "arm", "model", "median", "median_lcl",
                      "median_ucl", "mean", "mean_lcl", "mean_ucl", "rmst",
                      "rmst_lcl", "rmst_ucl"))
  expect_true(all(is.na(out$rmst)))
  expect_identical(out$arm[out$fit_type == "dependent" &
                             out$model == "Gamma"],
                   c("Control", "Treatment"))
  x$fits <- x$fits[2]
  x$models <- x$models[2, ]
  withr::local_seed(3)
  state <- .Random.seed
  ci <- survival_summary(x, time_horizon = 60, conf_int = TRUE,
                         n_sim = 200)
  expect_identical(.Random.seed, state)
  expect_true(ci$median_lcl < ci$median && ci$median < ci$median_ucl)
  expect_true(ci$rmst_lcl < ci$rmst && ci$rmst < ci$rmst_ucl)
  expect_identical(ci, survival_summary(x, time_horizon = 60,
                                        conf_int = TRUE, n_sim = 200))
  expect_error(survival_summary(x, time_horizon = 0), "'time_horizon'")
  expect_error(survival_summary(x, n_sim = 0), "'n_sim'")
  expect_error(survival_summary(list()), "fit_survival_models")
})
