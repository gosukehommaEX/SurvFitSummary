# Expected values: inst/validation/python/expected_km_cox.py (score test for
# the coefficient of ARM * g(t) with g(t) = 1 - KM(t-) centred over the
# events, at the Cox estimate, with the Efron approximation)

test_that("the test statistic matches independent values", {
  out <- test_proportional_hazards(make_test_data())
  expect_identical(out$comparison, "Drug vs Control")
  expect_equal(out$log_hr, -0.6843892641570603, tolerance = 1e-6)
  expect_equal(out$chisq, 0.007789449949722455, tolerance = 1e-5)
  expect_identical(out$df, 1)
  expect_equal(out$p_value, 0.9296717306120232, tolerance = 1e-5)
})

test_that("each arm is compared with the control arm", {
  d <- make_test_data()
  d3 <- rbind(d, transform(d[d$ARM == "Drug", ], ARM = "Other"))
  d3$ARM <- factor(as.character(d3$ARM),
                   levels = c("Control", "Drug", "Other"))
  out <- test_proportional_hazards(d3)
  expect_identical(out$comparison, c("Drug vs Control", "Other vs Control"))
  expect_equal(out$chisq[1], out$chisq[2], tolerance = 1e-10)
})

test_that("invalid inputs are rejected", {
  d <- make_test_data()
  expect_error(test_proportional_hazards(d[d$ARM == "Drug", ]),
               "At least 2 treatment arms")
  expect_error(test_proportional_hazards(d, transform = "x"))
})
