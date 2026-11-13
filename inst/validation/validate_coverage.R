# Coverage of the confidence intervals of the package, by simulation (INFO).
#
# Samples of 300 patients from three distributions with time in months,
# each with censoring times uniform on 12 to 48 months: Weibull (shape 1.3,
# scale 12), lognormal (meanlog 2.3, sdlog 0.9) and log-logistic (shape
# 1.8, scale 10). Each sample is fitted with its own distribution and, for
# the Weibull samples, also with the spline model on the hazard scale with
# 1 knot, which contains the Weibull model (gamma0 = -shape * log(scale),
# gamma1 = shape, gamma2 = 0, whatever the knot). For each sample the
# script records whether the 95% confidence intervals contain the true
# values: the parameters on the natural scale (parameter_table(): Wald
# intervals on the transformed scale) and the survival at 6, 12 and 24
# months (predict_survival(conf_int = TRUE): limits from 1000 simulated
# parameter vectors, seed = index of the sample). The coverage of each
# interval over the samples is reported as INFO, with the Monte Carlo
# standard error sqrt(0.95 * 0.05 / 500) = 0.0097 of a coverage of 0.95:
# the intervals are asymptotic, so values near 0.95 are expected but not
# required. A sample whose fit fails or warns is left out and counted.
#
# Run from the root of the package after devtools::load_all():
#   source(file.path("inst", "validation", "validate_coverage.R"))
# It takes several minutes. Output:
# inst/validation/results/validate_coverage.csv

source(file.path("inst", "validation", "verdict.R"))
script <- "validate_coverage"
n_samples <- 500
n_patients <- 300
times <- c(6, 12, 24)
truths <- list(
  list(label = "Weibull", code = "weibull",
       draw = function(n) stats::rweibull(n, shape = 1.3, scale = 12),
       surv = function(t) {
         stats::pweibull(t, shape = 1.3, scale = 12, lower.tail = FALSE)
       },
       pars = c(shape = 1.3, scale = 12),
       spline_pars = c(gamma0 = -1.3 * log(12), gamma1 = 1.3, gamma2 = 0)),
  list(label = "Lognormal", code = "lnorm",
       draw = function(n) stats::rlnorm(n, meanlog = 2.3, sdlog = 0.9),
       surv = function(t) {
         stats::plnorm(t, meanlog = 2.3, sdlog = 0.9, lower.tail = FALSE)
       },
       pars = c(meanlog = 2.3, sdlog = 0.9), spline_pars = NULL),
  list(label = "Log-logistic", code = "llogis",
       draw = function(n) flexsurv::rllogis(n, shape = 1.8, scale = 10),
       surv = function(t) {
         flexsurv::pllogis(t, shape = 1.8, scale = 10, lower.tail = FALSE)
       },
       pars = c(shape = 1.8, scale = 10), spline_pars = NULL)
)

set.seed(20261113)
result <- NULL
for (truth in truths) {
  models <- c(truth$label, if (!is.null(truth$spline_pars)) "Hazards, 1 knot")
  hits <- list()
  failed <- stats::setNames(integer(length(models)), models)
  for (b in seq_len(n_samples)) {
    event_time <- truth$draw(n_patients)
    censor_time <- stats::runif(n_patients, 12, 48)
    sample_data <- processing_dataset(
      data.frame(ARM = "Sample", TIME = pmin(event_time, censor_time),
                 EVENT = as.integer(event_time <= censor_time)),
      column_arm = "ARM", column_survtime = "TIME", column_event = "EVENT",
      control_arm = "Sample"
    )
    fits <- fit_survival_models(
      sample_data, distributions = truth$code,
      spline_scales = if (is.null(truth$spline_pars)) NULL else "hazard",
      spline_knots = 1L, fit_type = "independent"
    )
    par <- parameter_table(fits, scale = "natural")
    pred <- predict_survival(fits, times, conf_int = TRUE, n_sim = 1000,
                             seed = b)
    for (m in models) {
      if (fits$models$status[fits$models$model == m] != "ok") {
        failed[[m]] <- failed[[m]] + 1L
        next
      }
      true_pars <- if (m == truth$label) truth$pars else truth$spline_pars
      p <- par[par$model == m, ]
      p <- p[match(names(true_pars), p$parameter), ]
      s <- pred[pred$model == m, ]
      s <- s[match(times, s$time), ]
      true_surv <- truth$surv(times)
      hits[[m]] <- rbind(hits[[m]], c(
        stats::setNames(p$lcl <= true_pars & true_pars <= p$ucl,
                        names(true_pars)),
        stats::setNames(s$lcl <= true_surv & true_surv <= s$ucl,
                        paste0("S(", times, ")"))
      ))
    }
    if (b %% 100L == 0L) {
      message(truth$label, ": ", b, " of ", n_samples, " samples")
    }
  }
  for (m in models) {
    prefix <- paste0(truth$label, " samples, ", m, " model: ")
    result <- rbind(result, verdict_table(
      script,
      c(paste0(prefix, "coverage of the 95% CI of ", colnames(hits[[m]])),
        paste0(prefix, "samples left out (fit failed or warned)")),
      NA, c(colMeans(hits[[m]]), failed[[m]]), NA
    ))
  }
}
result <- rbind(result, verdict_table(
  script, paste("Monte Carlo standard error of a coverage of 0.95 with",
                n_samples, "samples"),
  NA, sqrt(0.95 * 0.05 / n_samples), NA
))

write_verdicts(result, file.path("inst", "validation", "results"))
