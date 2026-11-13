# The package against closed forms (base R only, no flexsurv function).
#
# 1. Exponential model: the maximum likelihood estimates are closed form.
#    For each arm, rate = d / T (d events, T total time at risk), standard
#    error of log(rate) = 1 / sqrt(d) and log-likelihood d log(rate) - d.
#    For the arms fitted together, the effect of the arm on log(rate) is
#    log((d1 / T1) / (d0 / T0)) with standard error sqrt(1 / d0 + 1 / d1).
# 2. With the estimated parameters of each model, the survival, hazard,
#    median, mean and restricted mean of the package against their formulas
#    (exponential, Weibull, lognormal, log-logistic, Gompertz, gamma); the
#    restricted mean by stats::integrate() of the closed-form survival.
# Data: inst/extdata/fit_example.csv (two arms of 300 patients).
#
# Run from the root of the package after devtools::load_all():
#   source(file.path("inst", "validation", "validate_closed_forms.R"))
# Output: inst/validation/results/validate_closed_forms.csv

source(file.path("inst", "validation", "verdict.R"))
script <- "validate_closed_forms"

dat <- processing_dataset(
  file.path("inst", "extdata", "fit_example.csv"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
distributions <- c("exp", "weibull", "lnorm", "llogis", "gompertz", "gamma")
fits <- fit_survival_models(dat, distributions = distributions,
                            spline_scales = NULL)
par <- parameter_table(fits, scale = "natural")
par_t <- parameter_table(fits, scale = "transformed")
gof <- gof_table(fits)
arms <- levels(dat$ARM)

# 1. Exponential model
result <- NULL
d <- tapply(dat$EVENT, dat$ARM, sum)
total <- tapply(dat$SURVTIME, dat$ARM, sum)
rate <- d / total
pick <- function(table, fit_type, arm, model, parameter, column) {
  table[[column]][table$fit_type == fit_type & table$arm == arm &
                    table$model == model & table$parameter == parameter]
}
for (a in arms) {
  result <- rbind(result, verdict_table(
    script,
    paste(a, "exponential:", c("rate", "se of log(rate)",
                               "log-likelihood")),
    c(rate[[a]], 1 / sqrt(d[[a]]), d[[a]] * log(rate[[a]]) - d[[a]]),
    c(pick(par, "independent", a, "Exponential", "rate", "est"),
      pick(par_t, "independent", a, "Exponential", "rate", "se"),
      gof$loglik[gof$fit_type == "independent" & gof$arm == a &
                   gof$model == "Exponential"]),
    c(1e-6 * rate[[a]], 1e-4 / sqrt(d[[a]]), 1e-6)
  ))
}
effect <- paste0("ARM", arms[2L])
result <- rbind(result, verdict_table(
  script,
  paste("Both arms exponential:", c("effect on log(rate)", "its se")),
  c(log(rate[[2L]] / rate[[1L]]), sqrt(1 / d[[1L]] + 1 / d[[2L]])),
  c(pick(par, "dependent", "Both arms", "Exponential", effect, "est"),
    pick(par, "dependent", "Both arms", "Exponential", effect, "se")),
  c(1e-6, 1e-4)
))

# 2. Survival, hazard, median, mean and restricted mean from the formulas
formulas <- list(
  Exponential = function(p) list(
    surv = function(t) exp(-p$rate * t),
    haz = function(t) rep(p$rate, length(t)),
    median = log(2) / p$rate, mean = 1 / p$rate),
  Weibull = function(p) list(
    surv = function(t) exp(-(t / p$scale)^p$shape),
    haz = function(t) p$shape / p$scale * (t / p$scale)^(p$shape - 1),
    median = p$scale * log(2)^(1 / p$shape),
    mean = p$scale * gamma(1 + 1 / p$shape)),
  Lognormal = function(p) list(
    surv = function(t) stats::pnorm((log(t) - p$meanlog) / p$sdlog,
                                    lower.tail = FALSE),
    haz = function(t) {
      z <- (log(t) - p$meanlog) / p$sdlog
      stats::dnorm(z) / (t * p$sdlog * stats::pnorm(z, lower.tail = FALSE))
    },
    median = exp(p$meanlog), mean = exp(p$meanlog + p$sdlog^2 / 2)),
  "Log-logistic" = function(p) list(
    surv = function(t) 1 / (1 + (t / p$scale)^p$shape),
    haz = function(t) {
      (p$shape / p$scale) * (t / p$scale)^(p$shape - 1) /
        (1 + (t / p$scale)^p$shape)
    },
    median = p$scale,
    mean = if (p$shape > 1) {
      p$scale * (pi / p$shape) / sin(pi / p$shape)
    } else {
      NA
    }),
  Gompertz = function(p) list(
    surv = function(t) exp(-p$rate / p$shape * (exp(p$shape * t) - 1)),
    haz = function(t) p$rate * exp(p$shape * t),
    median = log(1 + p$shape * log(2) / p$rate) / p$shape,
    mean = NA),
  Gamma = function(p) list(
    surv = function(t) stats::pgamma(t, p$shape, p$rate, lower.tail = FALSE),
    haz = function(t) {
      stats::dgamma(t, p$shape, p$rate) /
        stats::pgamma(t, p$shape, p$rate, lower.tail = FALSE)
    },
    median = stats::qgamma(0.5, p$shape, p$rate), mean = p$shape / p$rate)
)
times <- c(1, 6, 12, 24, 60, 120)
horizon <- 60
pred_s <- predict_survival(fits, times, type = "survival")
pred_h <- predict_survival(fits, times, type = "hazard")
summary <- survival_summary(fits, time_horizon = horizon)
for (model in names(formulas)) {
  for (a in arms) {
    rows <- par$fit_type == "independent" & par$arm == a & par$model == model
    p <- as.list(stats::setNames(par$est[rows], par$parameter[rows]))
    f <- formulas[[model]](p)
    s_rows <- pred_s$fit_type == "independent" & pred_s$arm == a &
      pred_s$model == model
    h_rows <- pred_h$fit_type == "independent" & pred_h$arm == a &
      pred_h$model == model
    m_row <- summary$fit_type == "independent" & summary$arm == a &
      summary$model == model
    rmst <- stats::integrate(f$surv, 0, horizon, rel.tol = 1e-10)$value
    prefix <- paste0(a, " ", model, ": ")
    result <- rbind(
      result,
      verdict_table(script, paste0(prefix, "S(", times, ")"),
                    f$surv(times), pred_s$est[s_rows][match(times,
                                                       pred_s$time[s_rows])],
                    1e-8),
      verdict_table(script, paste0(prefix, "h(", times, ")"),
                    f$haz(times), pred_h$est[h_rows][match(times,
                                                      pred_h$time[h_rows])],
                    1e-8 * pmax(1, f$haz(times))),
      verdict_table(script, paste0(prefix, c("median", "mean",
                                             "restricted mean to 60")),
                    c(f$median, f$mean, rmst),
                    c(summary$median[m_row], summary$mean[m_row],
                      summary$rmst[m_row]),
                    1e-5 * pmax(1, abs(c(f$median, f$mean, rmst)),
                                na.rm = TRUE))
    )
  }
}

write_verdicts(result, file.path("inst", "validation", "results"))
