# Expected values: inst/validation/python/expected_km_cox.py (numpy, scipy)

test_that("counts and medians match independent values", {
  out <- summarise_data(make_test_data())
  expect_identical(out$ARM, c("Control", "Drug", "Overall"))
  expect_identical(out$n, c(13L, 12L, 25L))
  expect_identical(out$events, c(9L, 7L, 16L))
  expect_identical(out$censored, c(4L, 5L, 9L))
  expect_identical(out$median, c(11, 16, 12))
  expect_identical(out$median_lcl, c(6, 10, 7))
  expect_identical(out$median_ucl, c(NA_real_, NA_real_, NA_real_))
  expect_identical(out$max_followup, c(18, 22, 22))
})

test_that("a single arm has no overall row and invalid inputs are rejected", {
  d <- make_test_data()
  d <- d[d$ARM == "Drug", ]
  out <- summarise_data(d)
  expect_identical(out$ARM, "Drug")
  expect_identical(out$median, 16)
  expect_error(summarise_data(make_test_raw()), "processing_dataset")
  expect_error(summarise_data(make_test_data(), conf_level = 2),
               "'conf_level'")
})
