test_that("the output has the documented structure", {
  out <- generate_dummy_survival_data(arm = c("A", "B"), n = c(30, 20),
                                      hazards = c(0.1, 0.05),
                                      dropout_per_year = 0.2, seed = 1)
  expect_named(out, c("SUBJID", "ARM", "SURVTIME", "CNSR"))
  expect_identical(nrow(out), 50L)
  expect_identical(out$ARM, rep(c("A", "B"), c(30, 20)))
  expect_true(all(out$SURVTIME > 0))
  expect_true(all(out$CNSR %in% c(0, 1)))
})

test_that("the seed gives reproducible data and the session RNG is kept", {
  args <- list(arm = c("A", "B"), n = c(30, 20), hazards = c(0.1, 0.05),
               dropout_per_year = 0.2, seed = 42)
  expect_identical(do.call(generate_dummy_survival_data, args),
                   do.call(generate_dummy_survival_data, args))
  set.seed(7)
  expected <- stats::runif(3)
  set.seed(7)
  invisible(do.call(generate_dummy_survival_data, args))
  expect_identical(stats::runif(3), expected)
})

test_that("the proportion censored agrees with the competing exponentials", {
  # P(dropout before event) = r / (r + lambda) with r = -log(1 - p) / 12
  lambda <- 0.05
  p <- 0.3
  r <- -log(1 - p) / 12
  out <- generate_dummy_survival_data(arm = "A", n = 20000, hazards = lambda,
                                      dropout_per_year = p, seed = 2024)
  expect_equal(mean(out$CNSR), r / (r + lambda), tolerance = 0.04)
  # The observed time is exponential with rate lambda + r
  expect_equal(mean(out$SURVTIME), 1 / (lambda + r), tolerance = 0.03)
})

test_that("invalid inputs are rejected", {
  expect_error(generate_dummy_survival_data("A", c(1, 2), 0.1, 0.1, 1),
               "'arm' must match length of 'n'")
  expect_error(generate_dummy_survival_data("A", 1, c(0.1, 0.2), 0.1, 1),
               "'arm' must match length of 'hazards'")
  expect_error(generate_dummy_survival_data("A", 1, 0.1, 1, 1),
               "'dropout_per_year'")
})
