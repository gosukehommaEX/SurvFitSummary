# Expected values: inst/validation/python/expected_fits.py (inverse of a
# numerical Hessian with Richardson extrapolation; delta method for the
# natural scale). The Cholesky factor is checked by L %*% t(L).

check_covariance <- function(set) {
  x <- get_fits(set)
  out <- covariance_table(x)
  expected <- expected_for("covariance", x)
  expect_identical(nrow(out), nrow(expected))
  e <- lookup(out, expected,
              c("fit_type", "arm", "model", "scale", "row", "col"))
  # Bound of each element: relative to sqrt(V_ii * V_jj) of the reference,
  # plus, on the natural scale, the effect of the tolerance on the estimates
  # of log-transformed parameters (the factors exp(est) of the delta method)
  key <- function(ft, arm, model, scale, r, c) {
    paste(ft, arm, model, scale, r, c, sep = "\t")
  }
  e_key <- key(expected$fit_type, expected$arm, expected$model,
               expected$scale, expected$row, expected$col)
  v_row <- expected$value[match(key(out$fit_type, out$arm, out$model,
                                    out$scale, out$row, out$row), e_key)]
  v_col <- expected$value[match(key(out$fit_type, out$arm, out$model,
                                    out$scale, out$col, out$col), e_key)]
  par <- expected_for("parameters", x)
  tr <- par[par$scale == "transformed", ]
  na <- par[par$scale == "natural", ]
  p_key <- function(d, p) paste(d$fit_type, d$arm, d$model, p, sep = "\t")
  na <- na[match(p_key(tr, tr$parameter), p_key(na, na$parameter)), ]
  log_se <- ifelse(abs(na$est - tr$est) > 1e-12, tr$se, 0)
  s_row <- log_se[match(p_key(out, out$row), p_key(tr, tr$parameter))]
  s_col <- log_se[match(p_key(out, out$col), p_key(tr, tr$parameter))]
  shift <- ifelse(out$scale == "natural", TOL_EST * (s_row + s_col), 0)
  bound <- (TOL_SE + shift) * sqrt(v_row * v_col)
  expect_lte(scaled_diff(out$value, e$value, bound), 1)
  # Cholesky factors: lower triangular with L %*% t(L) = covariance
  for (s in c("transformed", "natural")) {
    cv <- covariance_table(x, scale = s)
    ch <- covariance_table(x, type = "cholesky", scale = s)
    expect_identical(ch[c("fit_type", "arm", "model", "row", "col")],
                     cv[c("fit_type", "arm", "model", "row", "col")])
    group <- paste(cv$fit_type, cv$arm, cv$model)
    for (g in unique(group)) {
      rows <- group == g
      p <- round(sqrt(sum(rows)))
      v <- matrix(cv$value[rows], p, p, byrow = TRUE)
      l <- matrix(ch$value[rows], p, p, byrow = TRUE)
      expect_true(all(l[upper.tri(l)] == 0) && all(diag(l) > 0))
      expect_equal(l %*% t(l), v, tolerance = 1e-10)
    }
  }
}

test_that("standard distributions match the independent covariances", {
  check_covariance("simple")
})

test_that("generalized gamma distributions match the independent values", {
  check_covariance("gengamma")
})

test_that("generalized F (Prentice) distributions match the independent covariances", {
  check_covariance("genf")
})

test_that("original generalized F distributions match the independent covariances", {
  check_covariance("genf_orig")
})

test_that("hazard spline models match the independent covariances", {
  check_covariance("hazard")
})

test_that("odds spline models match the independent covariances", {
  check_covariance("odds")
})

test_that("normal spline models match the independent covariances", {
  check_covariance("normal")
})

test_that("layout, scale argument and missing covariance matrices", {
  x <- get_fits("simple")
  cv <- covariance_table(x, scale = "transformed")
  expect_named(cv, c("fit_type", "arm", "model", "scale", "row", "col",
                     "value"))
  first <- cv$fit_type == "independent" & cv$arm == "Control" &
    cv$model == "Weibull"
  expect_identical(cv$row[first], c("shape", "shape", "scale", "scale"))
  expect_identical(cv$col[first], c("shape", "scale", "shape", "scale"))
  expect_identical(unique(cv$scale), "transformed")
  # A model without covariance matrix gives NA values
  i <- which(x$models$model == "Weibull" & x$models$arm == "Control")
  x$fits[[i]]$cov <- NA
  cv <- covariance_table(x, scale = "transformed")
  ch <- covariance_table(x, type = "cholesky", scale = "natural")
  expect_true(all(is.na(cv$value[first])))
  expect_true(all(is.na(ch$value[first])))
  expect_false(anyNA(cv$value[!first]))
  pt <- parameter_table(x, scale = "transformed")
  weibull <- pt$arm == "Control" & pt$model == "Weibull"
  expect_true(all(is.na(pt$se[weibull])) && !anyNA(pt$est[weibull]))
  expect_error(covariance_table(x, type = "correlation"))
})
