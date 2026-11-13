# The fit statistics are those of gof_table() (compared with
# inst/validation/python/expected_fits.py in test-gof_table.R).

test_that("the summary gives the fit statistics of all the models", {
  x <- get_fits("simple")
  s <- summary(x)
  gof <- gof_table(x)
  expect_identical(names(s), c("fit_type", "arm", "model", "status", "n",
                               "events", "npars", "loglik", "AIC", "BIC"))
  expect_identical(nrow(s), nrow(x$models))
  expect_equal(s[names(gof)], gof)
  expect_identical(unique(s$status), "ok")
})

test_that("models that failed are kept with missing statistics", {
  d <- make_fit_data()
  d$EVENT[d$ARM == "Treatment"] <- 0L
  d$EVENT[d$ARM == "Treatment"][1] <- 1L
  x <- fit_survival_models(d, distributions = "exp",
                           spline_scales = "hazard", spline_knots = 1,
                           fit_type = "independent")
  s <- summary(x)
  expect_identical(nrow(s), 4L)
  failed <- which(s$status == "error")
  expect_length(failed, 1)
  expect_true(all(is.na(s[failed, c("npars", "loglik", "AIC", "BIC")])))
  expect_identical(s$n[failed], sum(d$ARM == s$arm[failed]))
  # The failed model is the last of its arm
  expect_identical(max(which(s$arm == s$arm[failed])), failed)
  ok <- s[-failed, names(gof_table(x))]
  rownames(ok) <- NULL
  expect_equal(ok, gof_table(x))
})

test_that("invalid inputs are rejected", {
  expect_error(summary.survfitsummary(list()), "'object'")
})
