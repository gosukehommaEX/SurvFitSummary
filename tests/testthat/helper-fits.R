# Dataset and reference values of the tests of the model fitting functions.
# inst/extdata/fit_example.csv was simulated by
# inst/validation/python/make_fit_example.py, and the expected values in
# inst/validation/python/expected/ were computed independently by
# inst/validation/python/expected_fits.py (numpy and scipy, no survival
# analysis library): maximum likelihood fits of the 12 distributions and of
# the Royston-Parmar models, covariance matrices from a numerical Hessian,
# predictions, medians, means and restricted means, and Kaplan-Meier
# estimates at the knots.

make_fit_data <- function() {
  processing_dataset(
    system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
    column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
    control_arm = "Control"
  )
}

read_expected <- function(name) {
  utils::read.csv(
    system.file("validation", "python", "expected", paste0(name, ".csv"),
                package = "SurvFitSummary"),
    stringsAsFactors = FALSE
  )
}

# Sets of models fitted once per run of the tests and shared by the test
# files (each set takes a few seconds). Each fits every arm separately and
# both arms with the arm as a covariate.
fit_sets <- list(
  simple = list(distributions = c("exp", "weibull", "weibullPH", "rayleigh",
                                  "lnorm", "llogis", "gompertz", "gamma"),
                spline_scales = NULL),
  gengamma = list(distributions = c("gengamma", "gengamma.orig"),
                  spline_scales = NULL),
  genf = list(distributions = "genf", spline_scales = NULL),
  genf_orig = list(distributions = "genf.orig", spline_scales = NULL),
  hazard = list(distributions = NULL, spline_scales = "hazard"),
  odds = list(distributions = NULL, spline_scales = "odds"),
  normal = list(distributions = NULL, spline_scales = "normal"),
  report = list(distributions = c("exp", "weibull"),
                spline_scales = "hazard", spline_knots = 1L)
)
fit_cache <- new.env(parent = emptyenv())
get_fits <- function(set) {
  if (is.null(fit_cache[[set]])) {
    fit_cache[[set]] <- do.call(fit_survival_models,
                                c(list(make_fit_data()), fit_sets[[set]]))
  }
  fit_cache[[set]]
}

# Tolerances of the comparisons with the reference values. flexsurv (optim
# with reltol = 1e-10) and the reference (Newton steps to a negligible
# gradient) reach the maximum to within about 1e-3 standard errors; the
# bounds leave a margin and are far smaller than the differences produced by
# a wrong parameterisation. Estimates are compared in units of their
# standard error, because the precision of an optimum is relative to it.
TOL_EST <- 0.02     # estimates on the transformed scale, in standard errors
TOL_CI <- 0.03      # confidence limits on the transformed scale, in SEs
TOL_SE <- 5e-3      # standard errors and covariances, relative
TOL_LOGLIK <- 1e-3  # log-likelihood, AIC and BIC, absolute
TOL_SURV <- 1e-3    # survival probabilities, absolute
TOL_HAZ <- 1e-2     # hazards, relative
TOL_TIME <- 2e-3    # median, mean and restricted mean, relative

# Rows of 'expected' in the order of the rows of 'actual', matched on the
# key columns; a row without a match gives NA values (and fails the tests)
lookup <- function(actual, expected, keys) {
  key <- function(d) {
    do.call(paste, c(lapply(d[keys], as.character), sep = "\t"))
  }
  expected[match(key(actual), key(expected)), , drop = FALSE]
}

# Largest difference in units of 'bound'; the tests require at most 1
scaled_diff <- function(actual, expected, bound) {
  max(abs(actual - expected) / bound)
}

# Expected values of the models of one fitted set
expected_for <- function(name, x) {
  e <- read_expected(name)
  e[e$model %in% unique(x$models$model), , drop = FALSE]
}
