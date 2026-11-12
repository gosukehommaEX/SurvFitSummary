# Expected values: inst/validation/python/expected_km_cox.py (numpy, scipy)

expected_km <- data.frame(
  ARM = rep(c("Control", "Drug"), c(4, 5)),
  time = c(0, 5, 10, 15, 0, 5, 10, 15, 20),
  n_risk = c(13L, 10L, 5L, 2L, 12L, 10L, 7L, 4L, 2L),
  n_event = c(0L, 3L, 3L, 3L, 0L, 1L, 3L, 1L, 2L),
  n_censor = c(0L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L),
  surv = c(1.0, 0.7615384615384615, 0.5076923076923077, 0.1523076923076923,
           1.0, 0.9166666666666666, 0.6285714285714286, 0.5238095238095238,
           0.19642857142857145),
  se = c(0.0, 0.12064347275897093, 0.1441816484401799, 0.12869648641959933,
         0.0, 0.07978559231302818, 0.14925796376517905, 0.15689701101567124,
         0.16114949576533819)
)

ci_log <- cbind(
  c(1.0, 0.5582694714290025, 0.290980790491762, 0.029072084366252637,
    1.0, 0.7729009736896152, 0.3946679403865125, 0.29121412266232216,
    0.03934423425123576),
  c(1.0, 1.0, 0.88580238872242, 0.7979349827088035,
    1.0, 1.0, 1.0, 0.9421810134933385, 0.9806820340456239)
)
ci_loglog <- cbind(
  c(1.0, 0.42672762345413784, 0.21419994068863096, 0.010703498357874966,
    1.0, 0.5389771805707411, 0.282209586825932, 0.20128131147977443,
    0.012636368737650843),
  c(1.0, 0.9165482581641975, 0.7421332424388727, 0.4581671444003902,
    1.0, 0.987825565400458, 0.8433235785043969, 0.7704105847113921,
    0.5455683625996981)
)
ci_plain <- cbind(
  c(1.0, 0.5250815999610394, 0.2251014695179394, 0.0,
    1.0, 0.7602897792479356, 0.3360311951858932, 0.21629703293682406, 0.0),
  c(1.0, 0.9979953231158836, 0.7902831458666759, 0.40454817062695514,
    1.0, 1.0, 0.9211116619569639, 0.8313220146822236, 0.5122757792554242)
)

test_that("counts, estimates and standard errors match independent values", {
  out <- km_table(make_test_data(), interval = 5)
  expect_identical(out$ARM, expected_km$ARM)
  expect_equal(out$time, expected_km$time, tolerance = 1e-12)
  expect_identical(out$n_risk, expected_km$n_risk)
  expect_identical(out$n_event, expected_km$n_event)
  expect_identical(out$n_censor, expected_km$n_censor)
  expect_equal(out$surv, expected_km$surv, tolerance = 1e-10)
  expect_equal(out$se, expected_km$se, tolerance = 1e-10)
})

test_that("the three confidence intervals match independent values", {
  d <- make_test_data()
  out <- km_table(d, interval = 5)
  expect_equal(cbind(out$lcl, out$ucl), ci_log, tolerance = 1e-10)
  out <- km_table(d, interval = 5, conf_type = "log-log")
  expect_equal(cbind(out$lcl, out$ucl), ci_loglog, tolerance = 1e-10)
  out <- km_table(d, interval = 5, conf_type = "plain")
  expect_equal(cbind(out$lcl, out$ucl), ci_plain, tolerance = 1e-10)
})

test_that("the last reporting time is kept despite rounding", {
  d <- make_test_data()
  d$SURVTIME <- d$SURVTIME / 10
  out <- km_table(d, interval = 0.1)
  expect_equal(max(out$time[out$ARM == "Drug"]), 2.2, tolerance = 1e-12)
})

test_that("invalid inputs are rejected", {
  expect_error(km_table(make_test_data(), interval = 0), "'interval'")
  expect_error(km_table(make_test_raw(), interval = 5), "processing_dataset")
})
