# Reproduction of Royston P, Parmar MKB (2002). Flexible parametric
# proportional-hazards and proportional-odds models for censored survival
# data, with application to prognostic modelling and estimation of
# treatment effects. Statistics in Medicine 21, 2175-2197.
#
# Breast cancer data of Section 2.1 (686 patients, 299 events; available as
# flexsurv::bc) with the three prognostic groups (reference
# Good). Spline models with d.f. = number of internal knots + 1 on the
# hazard scale (PH) and the odds scale (PO), with the group acting on
# gamma0 (models of all arms fitted together, fit_type = "dependent"). The
# models with 1 d.f. are the Weibull (PH) and log-logistic (PO) models.
#
# Published values reproduced (printed digits; tolerance half a unit of the
# last digit):
#   Table II: AIC and the coefficients of the Medium and Poor groups of the
#     PH and PO models with 1 to 4 d.f.
#   Section 5.1: knots of the PO(2) model (median, minimum and maximum
#     uncensored times) and its coefficients with standard errors.
#
# Time in years of 365.25 days (rectime / 365.25): the knots of Section 5.1
# (1.769, 0.197 and 6.724 years) are those of rectime / 365.25, whereas
# flexsurv::bc$recyrs is rectime / 365 (upper knot 6.729). The scale of time
# changes gamma0 but not the effects of the groups, nor the AIC of the
# article (log time).
#
# Two conversions, both exact:
#   1. The article gives the log-likelihood of log time: its AIC is the AIC
#      of the package minus 2 * (sum of the log event times).
#   2. The log-logistic model of the package has the group effect on log
#      scale (accelerated failure time); the PO coefficient is
#      -shape * (group effect). The Weibull (PH) model of the package gives
#      the PH coefficient directly.
#
# Run from the root of the package after devtools::load_all():
#   source(file.path("inst", "reproduce", "reproduce_RoystonParmar2002.R"))
# Output: inst/reproduce/results/reproduce_RoystonParmar2002.csv

source(file.path("inst", "validation", "verdict.R"))
script <- "reproduce_RoystonParmar2002"

bc <- flexsurv::bc
dat <- processing_dataset(
  data.frame(GROUP = bc$group, YEARS = bc$rectime / 365.25,
             EVENT = bc$censrec),
  column_arm = "GROUP", column_survtime = "YEARS", column_event = "EVENT",
  control_arm = "Good"
)
fits <- fit_survival_models(dat, distributions = c("weibullPH", "llogis"),
                            spline_scales = c("hazard", "odds"),
                            spline_knots = 1:3, fit_type = "dependent",
                            time_unit = "years")
par <- parameter_table(fits, scale = "natural")
gof <- gof_table(fits)
value_of <- function(model, parameter, column = "est") {
  par[[column]][par$model == model & par$parameter == parameter]
}
log_time_constant <- 2 * sum(log(dat$SURVTIME[dat$EVENT == 1]))

# Table II
models <- data.frame(
  scale = rep(c("PH", "PO"), each = 4),
  df = rep(1:4, times = 2),
  model = c("Weibull (PH)", "Hazards, 1 knot", "Hazards, 2 knots",
            "Hazards, 3 knots", "Log-logistic", "Odds, 1 knot",
            "Odds, 2 knots", "Odds, 3 knots"),
  aic = c(1284.9, 1248.8, 1248.0, 1247.0, 1259.1, 1241.0, 1241.5, 1240.6),
  beta1 = c(0.847, 0.834, 0.834, 0.836, 1.093, 1.052, 1.052, 1.054),
  beta2 = c(1.672, 1.612, 1.612, 1.613, 2.260, 2.171, 2.162, 2.159),
  stringsAsFactors = FALSE
)
coefficient <- function(model, parameter) {
  if (model == "Weibull (PH)") return(value_of(model, parameter))
  if (model == "Log-logistic") {
    return(-value_of(model, "shape") * value_of(model, parameter))
  }
  value_of(model, parameter)
}
label <- paste0("Table II, ", models$scale, "(", models$df, ")")
aic <- vapply(models$model, function(m) gof$AIC[gof$model == m],
              numeric(1)) - log_time_constant
result <- rbind(
  verdict_table(script, paste(label, "AIC"), models$aic, aic,
                printed_tolerance(1)),
  verdict_table(script, paste(label, "beta1 (Medium)"), models$beta1,
                vapply(models$model, coefficient, numeric(1),
                       parameter = "ARMMedium"),
                printed_tolerance(3)),
  verdict_table(script, paste(label, "beta2 (Poor)"), models$beta2,
                vapply(models$model, coefficient, numeric(1),
                       parameter = "ARMPoor"),
                printed_tolerance(3))
)

# Section 5.1: knots and coefficients of PO(2)
knots <- spline_knot_table(fits)
k <- knots[knots$n_knots == 1, ]
boundary <- range(log(dat$SURVTIME[dat$EVENT == 1]))
po2 <- "Odds, 1 knot"
result <- rbind(result, verdict_table(
  script,
  c("Section 5.1, PO(2) internal knot (time)",
    "Section 5.1, PO(2) lower boundary knot (time)",
    "Section 5.1, PO(2) upper boundary knot (time)",
    "Section 5.1, PO(2) lambda_1"),
  c(1.769, 0.197, 6.724, 0.378),
  c(exp(k$log_time), exp(boundary),
    (boundary[2] - k$log_time) / (boundary[2] - boundary[1])),
  c(printed_tolerance(3), printed_tolerance(3), printed_tolerance(3),
    printed_tolerance(3))
), verdict_table(
  script,
  paste("Section 5.1, PO(2)",
        rep(c("gamma0", "gamma1", "gamma2", "Medium", "Poor"), each = 2),
        rep(c("est", "se"), times = 5)),
  c(-3.451, 0.203, 2.915, 0.298, 0.191, 0.044, 1.052, 0.206, 2.171, 0.209),
  c(value_of(po2, "gamma0"), value_of(po2, "gamma0", "se"),
    value_of(po2, "gamma1"), value_of(po2, "gamma1", "se"),
    value_of(po2, "gamma2"), value_of(po2, "gamma2", "se"),
    value_of(po2, "ARMMedium"), value_of(po2, "ARMMedium", "se"),
    value_of(po2, "ARMPoor"), value_of(po2, "ARMPoor", "se")),
  printed_tolerance(3)
))
result <- rbind(result, verdict_table(
  script, "Constant 2 * sum(log event times) of the AIC conversion",
  NA, log_time_constant, NA
))

write_verdicts(result, file.path("inst", "reproduce", "results"))
