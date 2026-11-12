# Expected values: inst/validation/python/expected_km_cox.py (Cox partial
# likelihood with the Efron approximation, Newton-Raphson, numpy and scipy)

test_that("the four models match independent values", {
  out <- estimate_hazard_ratio(make_test_data(), covariates = c("AGE", "SEX"),
                               strata = "REGION")
  expect_identical(out$model, c("Unadjusted, unstratified",
                                "Unadjusted, stratified",
                                "Adjusted, unstratified",
                                "Adjusted, stratified"))
  expect_identical(out$covariates, c("None", "None", "AGE, SEX", "AGE, SEX"))
  expect_identical(out$strata, c("None", "REGION", "None", "REGION"))
  expect_identical(out$comparison, rep("Drug vs Control", 4))
  expect_equal(out$hr, c(0.5043981895776362, 0.5498916085777518,
                         0.4487517967934286, 0.5001280179135201),
               tolerance = 1e-6)
  expect_equal(out$lcl, c(0.17663440549777637, 0.19410183110499646,
                          0.1522848753364791, 0.17075670078459806),
               tolerance = 1e-6)
  expect_equal(out$ucl, c(1.4403622721871125, 1.5578461030625672,
                          1.3223780410259263, 1.4648211938554123),
               tolerance = 1e-6)
  expect_equal(out$p_value, c(0.20111724371965933, 0.26033692426768495,
                              0.14617016925558982, 0.20632471158818155),
               tolerance = 1e-6)
  expect_identical(out$n, rep(25L, 4))
  expect_identical(out$events, rep(16L, 4))
})

test_that("the control arm defines the direction of the hazard ratio", {
  out <- estimate_hazard_ratio(make_test_data(control_arm = "Drug"))
  expect_identical(out$comparison, "Control vs Drug")
  expect_equal(out$hr, 1 / 0.5043981895776362, tolerance = 1e-6)
})

test_that("invalid inputs are rejected", {
  d <- make_test_data()
  expect_error(estimate_hazard_ratio(d[d$ARM == "Drug", ]),
               "At least 2 treatment arms")
  expect_error(estimate_hazard_ratio(d, covariates = "WEIGHT"), "not found")
})
