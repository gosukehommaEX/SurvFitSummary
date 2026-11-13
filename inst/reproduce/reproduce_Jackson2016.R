# Reproduction of Jackson CH (2016). flexsurv: a platform for parametric
# survival modeling in R. Journal of Statistical Software 70(8), 1-33.
#
# Breast cancer data (flexsurv::bc: 686 patients, 299 recurrences or deaths,
# time in years), with the prognostic group (Good, Medium, Poor; reference
# Good) as covariate. In SurvFitSummary, the groups are the arms and the
# models with the group as covariate are the models of all arms fitted
# together (fit_type = "dependent"): the group acts on the location
# parameter of a parametric model and on gamma0 of a spline model, as in
# the article.
#
# Published values reproduced (printed digits; tolerance half a unit of the
# last digit):
#   Section 3.1, model fs1 (Weibull): estimates, 95% CI and standard errors
#     of shape, scale and the group effects; log-likelihood and AIC.
#   Table 3: -2 log-likelihood, number of parameters and AIC of fs1
#     (Weibull), fs2 (generalized gamma), fs4 (log-logistic), sp1 (spline,
#     odds scale, 1 knot) and sp3 (spline, hazard scale, 1 knot).
#   Section 5.1: group effects (log hazard ratios) and standard errors of
#     sp3, and the 5-year survival of the Good group under sp1. The article
#     computes the latter as 1 - psurvspline(5, gamma, knots), whose default
#     is scale = "hazard", for the model sp1 on the odds scale: the printed
#     0.6896969 is exp(-exp(eta)) with the linear predictor eta of sp1 at 5
#     years, not the survival 1 / (1 + exp(eta)) of sp1. Compared are the
#     survival of the package with the survival implied by the printed
#     value, 1 / (1 - log(0.6896969)), and exp(-exp(eta)) of the package
#     with the printed value.
#
# Run from the root of the package after devtools::load_all():
#   source(file.path("inst", "reproduce", "reproduce_Jackson2016.R"))
# Output: inst/reproduce/results/reproduce_Jackson2016.csv

source(file.path("inst", "validation", "verdict.R"))
script <- "reproduce_Jackson2016"

bc <- flexsurv::bc
dat <- processing_dataset(
  data.frame(GROUP = bc$group, YEARS = bc$recyrs, EVENT = bc$censrec),
  column_arm = "GROUP", column_survtime = "YEARS", column_event = "EVENT",
  control_arm = "Good"
)
fits <- fit_survival_models(dat,
                            distributions = c("weibull", "gengamma",
                                              "llogis"),
                            spline_scales = c("hazard", "odds"),
                            spline_knots = 1L, fit_type = "dependent",
                            time_unit = "years")
par <- parameter_table(fits, scale = "natural")
gof <- gof_table(fits)
value_of <- function(model, parameter, column) {
  par[[column]][par$model == model & par$parameter == parameter]
}
fit_of <- function(model, column) gof[[column]][gof$model == model]

# Section 3.1: model fs1
fs1 <- data.frame(
  parameter = rep(c("shape", "scale", "ARMMedium", "ARMPoor"), each = 4),
  column = rep(c("est", "lcl", "ucl", "se"), times = 4),
  published = c(1.3797, 1.2548, 1.5170, 0.0668,
                11.4229, 9.1818, 14.2110, 1.2728,
                -0.6136, -0.8623, -0.3649, 0.1269,
                -1.2122, -1.4583, -0.9661, 0.1256),
  stringsAsFactors = FALSE
)
fs1_value <- mapply(function(p, column) value_of("Weibull", p, column),
                    fs1$parameter, fs1$column)
result <- verdict_table(
  script, paste("fs1 (Weibull):", fs1$parameter, fs1$column),
  fs1$published, fs1_value, printed_tolerance(4)
)
result <- rbind(result, verdict_table(
  script, c("fs1 (Weibull): log-likelihood", "fs1 (Weibull): AIC",
            "fs1 (Weibull): number of patients", "fs1 (Weibull): events"),
  c(-811.9419, 1631.884, 686, 299),
  c(fit_of("Weibull", "loglik"), fit_of("Weibull", "AIC"),
    fit_of("Weibull", "n"), fit_of("Weibull", "events")),
  c(printed_tolerance(4), printed_tolerance(3), 0, 0)
))

# Table 3
table3 <- data.frame(
  label = c("fs1 (Weibull)", "fs2 (generalized gamma)",
            "fs4 (log-logistic)", "sp1 (spline, odds, 1 knot)",
            "sp3 (spline, hazard, 1 knot)"),
  model = c("Weibull", "Generalized Gamma", "Log-logistic",
            "Odds, 1 knot", "Hazards, 1 knot"),
  minus2ll = c(1623.8, 1575.2, 1598.2, 1578.0, 1585.8),
  npars = c(4, 5, 4, 5, 5),
  aic = c(1631.9, 1585.1, 1606.1, 1588.0, 1595.7),
  stringsAsFactors = FALSE
)
loglik <- vapply(table3$model, fit_of, numeric(1), column = "loglik")
# The article prints -2 * round(loglik, 1)
result <- rbind(result, verdict_table(
  script, paste("Table 3,", table3$label, "-2 log-likelihood"),
  table3$minus2ll, -2 * round(loglik, 1), 1e-9
), verdict_table(
  script, paste("Table 3,", table3$label, "parameters"),
  table3$npars, vapply(table3$model, fit_of, numeric(1), column = "npars"),
  0
), verdict_table(
  script, paste("Table 3,", table3$label, "AIC"),
  table3$aic, vapply(table3$model, fit_of, numeric(1), column = "AIC"),
  printed_tolerance(1)
))

# Section 5.1: sp3 group effects and the 5-year survival under sp1
result <- rbind(result, verdict_table(
  script, c("sp3: ARMMedium est", "sp3: ARMMedium se",
            "sp3: ARMPoor est", "sp3: ARMPoor se"),
  c(0.8334517, 0.1712042, 1.6111788, 0.1640933),
  c(value_of("Hazards, 1 knot", "ARMMedium", "est"),
    value_of("Hazards, 1 knot", "ARMMedium", "se"),
    value_of("Hazards, 1 knot", "ARMPoor", "est"),
    value_of("Hazards, 1 knot", "ARMPoor", "se")),
  printed_tolerance(7)
))
pred <- predict_survival(fits, times = 5)
s5 <- pred$est[pred$model == "Odds, 1 knot" & pred$arm == "Good"]
eta <- log((1 - s5) / s5)
# Survival of sp1 implied by the printed exp(-exp(eta)), with exp(eta)
# equal to -log(0.6896969). The rounding of the printed value moves it by
# at most 0.77 times printed_tolerance(7), which is used as the tolerance.
result <- rbind(result, verdict_table(
  script,
  c(paste("sp1: 5-year survival of the Good group (odds scale, from the",
          "printed value)"),
    "sp1: 5-year survival as computed in the article (hazard scale)"),
  c(1 / (1 - log(0.6896969)), 0.6896969), c(s5, exp(-exp(eta))),
  printed_tolerance(7)
))

write_verdicts(result, file.path("inst", "reproduce", "results"))
