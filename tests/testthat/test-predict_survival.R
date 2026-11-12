# Expected values: inst/validation/python/expected_fits.py (survival and
# hazard functions of the independent fits; Kaplan-Meier estimates).

times <- c(0, 6, 12, 24, 36, 60, 120)

check_predictions <- function(set) {
  x <- get_fits(set)
  expected <- expected_for("predictions", x)
  keys <- c("fit_type", "arm", "model", "time")
  s <- predict_survival(x, times)
  expect_identical(nrow(s), nrow(expected))
  e <- lookup(s, expected, keys)
  expect_lte(scaled_diff(s$est, e$surv, TOL_SURV), 1)
  expect_identical(is.na(s$km), is.na(e$km))
  expect_equal(s$km[!is.na(s$km)], e$km[!is.na(e$km)], tolerance = 1e-10)
  h <- predict_survival(x, times, type = "hazard")
  positive <- h$time > 0
  expect_lte(scaled_diff(h$est[positive], e$hazard[positive],
                         TOL_HAZ * e$hazard[positive]), 1)
  expect_true(all(is.na(h$est[!positive])))
  # Cumulative hazard against -log(S) of the reference, where S > 1e-10
  # (further in the tail flexsurv's S = 1 - F loses its relative precision)
  ch <- predict_survival(x, times, type = "cumhaz")
  used <- ch$time > 0 & e$surv > 1e-10
  expect_lte(scaled_diff(ch$est[used], -log(e$surv[used]),
                         TOL_HAZ * -log(e$surv[used])), 1)
  expect_true(all(ch$est[ch$time == 0] == 0))
}

test_that("standard distributions match the independent predictions", {
  check_predictions("simple")
})

test_that("generalized gamma distributions match the independent values", {
  check_predictions("gengamma")
})

test_that("generalized F (Prentice) distributions match the independent predictions", {
  check_predictions("genf")
})

test_that("original generalized F distributions match the independent predictions", {
  check_predictions("genf_orig")
})

test_that("hazard spline models match the independent predictions", {
  check_predictions("hazard")
})

test_that("odds spline models match the independent predictions", {
  check_predictions("odds")
})

test_that("normal spline models match the independent predictions", {
  check_predictions("normal")
})

test_that("layout, Kaplan-Meier after follow-up and values at time 0", {
  x <- get_fits("simple")
  s <- predict_survival(x, times)
  expect_named(s, c("fit_type", "arm", "model", "time", "est", "lcl", "ucl",
                    "km"))
  # Dependent fits are predicted for each arm
  dep <- s[s$fit_type == "dependent" & s$model == "Weibull", ]
  expect_identical(dep$arm, rep(c("Control", "Treatment"), each = 7))
  expect_true(all(s$est[s$time == 0] == 1))
  expect_true(all(is.na(s$lcl)) && all(is.na(s$ucl)))
  # The Kaplan-Meier estimate stops at the largest observed time of the arm
  d <- make_fit_data()
  for (a in c("Control", "Treatment")) {
    late <- s$arm == a & s$time > max(d$SURVTIME[d$ARM == a])
    expect_true(all(is.na(s$km[late])) && !anyNA(s$km[!late & s$arm == a]))
  }
  expect_true(all(predict_survival(x, 0, type = "cumhaz")$est == 0))
  expect_identical(predict_survival(x, c(12, 6))$est[1:2],
                   predict_survival(x, c(6, 12))$est[2:1])
})

test_that("simulated confidence limits are reproducible", {
  x <- get_fits("simple")
  x$fits <- x$fits[1:2]
  x$models <- x$models[1:2, ]
  withr::local_seed(5)
  state <- .Random.seed
  a <- predict_survival(x, c(0, 12, 24), conf_int = TRUE, n_sim = 200)
  expect_identical(.Random.seed, state)
  b <- predict_survival(x, c(0, 12, 24), conf_int = TRUE, n_sim = 200)
  expect_identical(a, b)
  d <- predict_survival(x, c(0, 12, 24), conf_int = TRUE, n_sim = 200,
                        seed = 2)
  expect_false(identical(a$lcl[a$time > 0], d$lcl[d$time > 0]))
  inside <- a$time > 0
  expect_true(all(a$lcl[inside] < a$est[inside] &
                    a$est[inside] < a$ucl[inside]))
  expect_true(all(a$lcl[!inside] == 1 & a$ucl[!inside] == 1))
})

test_that("invalid inputs are rejected", {
  x <- get_fits("simple")
  expect_error(predict_survival(x, -1), "'times'")
  expect_error(predict_survival(x, numeric(0)), "'times'")
  expect_error(predict_survival(x, 12, type = "density"))
  expect_error(predict_survival(x, 12, conf_int = NA), "'conf_int'")
  expect_error(predict_survival(x, 12, n_sim = 1), "'n_sim'")
  expect_error(predict_survival(x, 12, seed = 1.5), "'seed'")
  expect_error(predict_survival(list(), 12), "fit_survival_models")
})
