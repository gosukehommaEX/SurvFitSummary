# Parametric models of the package against survival::survreg().
#
# survreg() fits the exponential, Weibull, lognormal and log-logistic models
# as accelerated failure time models with its own parameterization (an
# intercept, the covariate effects and log(scale)). The fits of the package
# (flexsurv) are converted to that parameterization and compared: the
# log-likelihood, the parameters and their standard errors on the scale of
# estimation, for each arm alone (fit_type = "independent") and for the
# arms fitted together (fit_type = "dependent").
#   exponential: log(rate) = -intercept; arm effect on log(rate) = -beta
#   Weibull, log-logistic: log(shape) = -log(scale of survreg);
#     log(scale) = intercept; arm effect on log(scale) = beta
#   lognormal: meanlog = intercept; log(sdlog) = log(scale of survreg);
#     arm effect on meanlog = beta
# Data: survival::lung (arms: sex) and flexsurv::bc (arms: prognostic
# group). Tolerance: 1e-4 times max(1, |value|).
#
# Run from the root of the package after devtools::load_all():
#   source(file.path("inst", "validation", "validate_survreg.R"))
# Output: inst/validation/results/validate_survreg.csv

source(file.path("inst", "validation", "verdict.R"))
script <- "validate_survreg"

datasets <- list(
  lung = data.frame(ARM = factor(survival::lung$sex, levels = 1:2,
                                 labels = c("Male", "Female")),
                    TIME = survival::lung$time,
                    EVENT = as.integer(survival::lung$status == 2)),
  bc = data.frame(ARM = flexsurv::bc$group, TIME = flexsurv::bc$recyrs,
                  EVENT = flexsurv::bc$censrec)
)
distributions <- c(exp = "exponential", weibull = "weibull",
                   lnorm = "lognormal", llogis = "loglogistic")
labels <- c(exp = "Exponential", weibull = "Weibull", lnorm = "Lognormal",
            llogis = "Log-logistic")

# survreg() estimates in the parameterization of flexsurv on the scale of
# estimation, named as in parameter_table(scale = "transformed")
from_survreg <- function(fit, distribution) {
  b <- stats::coef(fit)
  v <- stats::vcov(fit)
  effects <- names(b)[-1L]
  log_scale <- if (distribution == "exp") NULL else "Log(scale)"
  # Linear map from (coefficients, log(scale)) to the flexsurv parameters
  names_in <- c(names(b), log_scale)
  est_in <- c(b, if (!is.null(log_scale)) log(fit$scale))
  if (distribution == "exp") {
    out <- c("rate", effects)
    a <- -diag(length(names_in))
  } else if (distribution %in% c("weibull", "llogis")) {
    out <- c("shape", "scale", effects)
    a <- matrix(0, length(out), length(names_in))
    a[1L, length(names_in)] <- -1
    a[2L, 1L] <- 1
    for (j in seq_along(effects)) a[2L + j, 1L + j] <- 1
  } else {
    out <- c("meanlog", "sdlog", effects)
    a <- matrix(0, length(out), length(names_in))
    a[1L, 1L] <- 1
    a[2L, length(names_in)] <- 1
    for (j in seq_along(effects)) a[2L + j, 1L + j] <- 1
  }
  est <- as.vector(a %*% est_in)
  se <- sqrt(diag(a %*% v[names_in, names_in] %*% t(a)))
  data.frame(parameter = out, est = est, se = se, stringsAsFactors = FALSE)
}

result <- NULL
for (data_name in names(datasets)) {
  raw <- datasets[[data_name]]
  dat <- processing_dataset(raw, column_arm = "ARM", column_survtime = "TIME",
                            column_event = "EVENT",
                            control_arm = levels(raw$ARM)[1L])
  fits <- fit_survival_models(dat, distributions = names(distributions),
                              spline_scales = NULL)
  par <- parameter_table(fits, scale = "transformed")
  gof <- gof_table(fits)
  groups <- c(lapply(levels(dat$ARM), function(a) {
    list(fit_type = "independent", arm = a, data = dat[dat$ARM == a, ],
         formula = survival::Surv(SURVTIME, EVENT) ~ 1)
  }), list(list(fit_type = "dependent", arm = "Both arms", data = dat,
                formula = survival::Surv(SURVTIME, EVENT) ~ ARM)))
  for (g in groups) {
    for (d in names(distributions)) {
      ref <- survival::survreg(g$formula, data = g$data,
                               dist = distributions[[d]])
      expected <- from_survreg(ref, d)
      rows <- par$fit_type == g$fit_type & par$arm == g$arm &
        par$model == labels[[d]]
      value <- par[rows, ]
      value <- value[match(expected$parameter, value$parameter), ]
      prefix <- paste0(data_name, ", ", g$arm, ", ", labels[[d]], ": ")
      ll <- gof$loglik[gof$fit_type == g$fit_type & gof$arm == g$arm &
                         gof$model == labels[[d]]]
      result <- rbind(
        result,
        verdict_table(script, paste0(prefix, "log-likelihood"),
                      ref$loglik[2L], ll, 1e-4 * max(1, abs(ll))),
        verdict_table(script, paste0(prefix, expected$parameter, " est"),
                      expected$est, value$est,
                      1e-4 * pmax(1, abs(expected$est))),
        verdict_table(script, paste0(prefix, expected$parameter, " se"),
                      expected$se, value$se,
                      1e-4 * pmax(1, abs(expected$se)))
      )
    }
  }
}

write_verdicts(result, file.path("inst", "validation", "results"))
