# Package against the independent reference implementation in Python.
#
# inst/validation/python/expected_fits.py (numpy and scipy, no survival
# analysis library) fits the 12 distributions and the Royston-Parmar models
# on the hazard, odds and normal scales with 1 to 3 knots to
# inst/extdata/fit_example.csv (each arm, and both arms with the arm as a
# covariate) and writes the reference values to python/expected/*.csv. The
# tests compare the same values (tests/testthat/helper-fits.R); this script
# writes the verdict of each value, with the tolerances of the tests, each
# capped at half the shift of the mutation run (see add()):
#   log-likelihood, AIC, BIC   1e-3 (spline models 1e-6) and twice that
#   estimates                  0.02 standard errors of the reference
#   confidence limits          0.03 standard errors (transformed scale)
#   standard errors            0.5%, relative
#   covariances                0.5% of sqrt(V_ii * V_jj) of the reference
#                              (transformed scale; the natural scale is
#                              checked by the tests)
#   survival                   1e-3, absolute, at 0, 6, 12, 24, 36, 60 and
#                              120 months
#   hazard                     1%, relative (times after 0)
#   Kaplan-Meier column        1e-10, and NA at the same rows
#   median, mean, RMST to 60   0.2%, relative (the mean where the reference
#                              flags it as reliable, mean_check = 1)
#   knots                      1e-10 relative (quantile, log time, time,
#                              Kaplan-Meier), 1e-8 (its standard error),
#                              exact (number at risk)
#   numbers of rows            exact
# The values are taken at the rows of the reference: a row that the package
# does not give has no value and FAILs.
#
# Run from the root of the package after devtools::load_all():
#   source(file.path("inst", "validation", "validate_python_reference.R"))
# Output: inst/validation/results/validate_python_reference.csv

source(file.path("inst", "validation", "verdict.R"))
script <- "validate_python_reference"

read_reference <- function(name) {
  utils::read.csv(file.path("inst", "validation", "python", "expected",
                            paste0(name, ".csv")),
                  stringsAsFactors = FALSE)
}
# Values of 'column' of 'actual' at the rows of 'reference', matched on the
# key columns (NA where the package has no such row)
at_reference <- function(actual, reference, keys, column) {
  key <- function(d) {
    do.call(paste, c(lapply(d[keys], as.character), sep = "\t"))
  }
  actual[[column]][match(key(reference), key(actual))]
}
label_of <- function(d, keys) {
  do.call(paste, c(lapply(d[keys], as.character), sep = ", "))
}
result <- NULL
# Each tolerance is at most half the shift of the mutation run (verdict.R),
# 0.005 * max(1, |reference|), so that an error of that size FAILs. Without
# this cap the tolerances in standard errors exceeded it for parameters whose
# standard error is large compared with their value (generalized F, gamma3
# and gamma4 of the spline models with 3 knots), and with the relative
# tolerance of 1% for a hazard above 1: 71 comparisons passed in the mutation
# run of 2026-11-13. The largest difference of the normal run was then 0.18
# of its capped tolerance.
add <- function(item, reference, value, tolerance) {
  tolerance <- pmin(tolerance, 0.005 * pmax(1, abs(reference)))
  result <<- rbind(result, verdict_table(script, item, reference, value,
                                         tolerance))
}

dat <- processing_dataset(
  file.path("inst", "extdata", "fit_example.csv"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(
  dat,
  distributions = c("exp", "weibull", "weibullPH", "rayleigh", "lnorm",
                    "llogis", "gompertz", "gamma", "gengamma",
                    "gengamma.orig", "genf", "genf.orig"),
  spline_scales = c("hazard", "odds", "normal"), spline_knots = 1:3
)
times <- c(0, 6, 12, 24, 36, 60, 120)
model_keys <- c("fit_type", "arm", "model")

# Fit statistics
ref <- read_reference("gof")
out <- gof_table(fits)
add("number of rows: fit statistics", nrow(ref), nrow(out), 0)
tol_loglik <- ifelse(grepl("knot", ref$model), 1e-6, 1e-3)
prefix <- paste0(label_of(ref, model_keys), ": ")
for (column in c("n", "events", "npars")) {
  add(paste0(prefix, column), ref[[column]],
      at_reference(out, ref, model_keys, column), 0)
}
add(paste0(prefix, "log-likelihood"), ref$loglik,
    at_reference(out, ref, model_keys, "loglik"), tol_loglik)
add(paste0(prefix, "AIC"), ref$AIC,
    at_reference(out, ref, model_keys, "AIC"), 2 * tol_loglik)
add(paste0(prefix, "BIC"), ref$BIC,
    at_reference(out, ref, model_keys, "BIC"), 2 * tol_loglik)

# Parameters on both scales
ref <- read_reference("parameters")
out <- parameter_table(fits)
add("number of rows: parameters", nrow(ref), nrow(out), 0)
keys <- c(model_keys, "parameter", "scale")
prefix <- paste0(label_of(ref, keys), ": ")
add(paste0(prefix, "est"), ref$est, at_reference(out, ref, keys, "est"),
    0.02 * ref$se)
add(paste0(prefix, "se"), ref$se, at_reference(out, ref, keys, "se"),
    5e-3 * ref$se)
tr <- ref[ref$scale == "transformed", ]
for (column in c("lcl", "ucl")) {
  add(paste0(label_of(tr, keys), ": ", column), tr[[column]],
      at_reference(out, tr, keys, column), 0.03 * tr$se)
}

# Covariance matrices on the transformed scale
ref <- read_reference("covariance")
ref <- ref[ref$scale == "transformed", ]
out <- covariance_table(fits, scale = "transformed")
add("number of rows: covariance (transformed scale)", nrow(ref), nrow(out),
    0)
keys <- c(model_keys, "scale", "row", "col")
diagonal <- function(parameter) {
  d <- ref[ref$row == ref$col, ]
  d$value[match(paste(label_of(ref, model_keys), parameter),
                paste(label_of(d, model_keys), d$row))]
}
add(paste0(label_of(ref, keys), ": covariance"), ref$value,
    at_reference(out, ref, keys, "value"),
    5e-3 * sqrt(diagonal(ref$row) * diagonal(ref$col)))

# Predictions
ref <- read_reference("predictions")
keys <- c(model_keys, "time")
s <- predict_survival(fits, times)
h <- predict_survival(fits, times, type = "hazard")
add("number of rows: predictions", nrow(ref), nrow(s), 0)
prefix <- paste0(label_of(ref, model_keys), ": ")
add(paste0(prefix, "S(", ref$time, ")"), ref$surv,
    at_reference(s, ref, keys, "est"), 1e-3)
positive <- ref$time > 0
add(paste0(prefix, "h(", ref$time, ")")[positive], ref$hazard[positive],
    at_reference(h, ref, keys, "est")[positive],
    1e-2 * ref$hazard[positive])
km <- at_reference(s, ref, keys, "km")
add("Kaplan-Meier column: rows with NA that differ from the reference", 0,
    sum(is.na(km) != is.na(ref$km)), 0)
observed <- !is.na(ref$km)
add(paste0(prefix, "Kaplan-Meier at ", ref$time)[observed],
    ref$km[observed], km[observed], 1e-10)

# Median, mean and restricted mean to 60
ref <- read_reference("summary")
out <- survival_summary(fits, time_horizon = 60)
add("number of rows: survival summary", nrow(ref), nrow(out), 0)
prefix <- paste0(label_of(ref, model_keys), ": ")
add(paste0(prefix, "median"), ref$median,
    at_reference(out, ref, model_keys, "median"), 2e-3 * ref$median)
add(paste0(prefix, "RMST to 60"), ref$rmst,
    at_reference(out, ref, model_keys, "rmst"), 2e-3 * ref$rmst)
checked <- ref$mean_check == 1
add(paste0(prefix, "mean")[checked], ref$mean[checked],
    at_reference(out, ref, model_keys, "mean")[checked],
    2e-3 * ref$mean[checked])

# Knots
ref <- read_reference("knots")
out <- spline_knot_table(fits)
add("number of rows: knots", nrow(ref), nrow(out), 0)
keys <- c("fit_type", "arm", "n_knots", "knot")
prefix <- paste0(label_of(ref, keys), ": ")
for (column in c("quantile", "log_time", "time", "surv_est")) {
  add(paste0(prefix, column), ref[[column]],
      at_reference(out, ref, keys, column),
      1e-10 * pmax(1, abs(ref[[column]])))
}
add(paste0(prefix, "surv_se"), ref$surv_se,
    at_reference(out, ref, keys, "surv_se"), 1e-8)
add(paste0(prefix, "n_risk"), ref$n_risk,
    at_reference(out, ref, keys, "n_risk"), 0)

write_verdicts(result, file.path("inst", "validation", "results"))
