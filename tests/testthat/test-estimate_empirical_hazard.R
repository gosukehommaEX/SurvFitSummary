# bshazard has no closed form, so the checks use data with a known constant
# hazard and the properties of the output.

simulate_exponential <- function(rate, n, seed) {
  set.seed(seed)
  event_time <- stats::rexp(n, rate)
  censor_time <- stats::runif(n, 0, 30)
  data.frame(
    ARM = factor(rep("A", n)),
    SURVTIME = pmin(event_time, censor_time),
    EVENT = as.integer(event_time <= censor_time)
  )
}

test_that("a constant hazard is recovered", {
  d <- simulate_exponential(rate = 0.1, n = 400, seed = 11)
  out <- estimate_empirical_hazard(d)
  expect_named(out, c("ARM", "time", "hazard", "lcl", "ucl"))
  middle <- out$time > 2 & out$time < 15
  expect_equal(mean(out$hazard[middle]), 0.1, tolerance = 0.25)
  expect_true(all(out$lcl <= out$hazard & out$hazard <= out$ucl))
})

test_that("each arm is smoothed separately and the level sets the width", {
  d <- rbind(simulate_exponential(0.1, 200, 1),
             transform(simulate_exponential(0.3, 200, 2), ARM = "B"))
  d$ARM <- factor(d$ARM, levels = c("A", "B"))
  out <- estimate_empirical_hazard(d)
  expect_identical(unique(out$ARM), c("A", "B"))
  expect_gt(mean(out$hazard[out$ARM == "B"]),
            mean(out$hazard[out$ARM == "A"]))
  narrow <- estimate_empirical_hazard(d, conf_level = 0.5)
  expect_true(all(narrow$ucl - narrow$lcl < out$ucl - out$lcl))
})

test_that("invalid inputs are rejected", {
  expect_error(estimate_empirical_hazard(make_test_raw()),
               "processing_dataset")
  expect_error(estimate_empirical_hazard(make_test_data(), conf_level = 1),
               "'conf_level'")
})
