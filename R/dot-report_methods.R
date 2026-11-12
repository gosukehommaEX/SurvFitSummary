#' Methods Section of the Word Report
#'
#' @description
#' Internal helper: the Methods section of the report template, with the
#' methods, settings and software versions of the analysis filled in. The
#' text of the template is used, with these corrections: the extrapolation
#' of restricted cubic splines beyond the last knot is linear in log time
#' on the scale of the model (log cumulative hazard, log cumulative odds or
#' normal quantile), not on the scale of the hazard; AIC and BIC are not
#' residuals.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param s List returned by \code{.report_settings()}.
#'
#' @return A list of elements made by \code{.report_block()}.
#'
#' @keywords internal
#' @noRd
.report_methods <- function(x, s) {
  b <- .report_block
  level <- paste0(format(100 * x$conf_level), "%")
  unit_one <- sub("s$", "", s$unit)
  two_arms <- length(x$arms) >= 2L
  m <- x$models
  blocks <- list(
    b("heading", "Methods", level = 1L),
    b("heading", "Kaplan-Meier analysis", level = 2L),
    b("paragraph", paste0(
      "The Kaplan-Meier (KM) method is used to estimate the KM ",
      "probabilities and median time to event. ", level, " confidence ",
      "intervals (CI) of the survival probabilities are computed with the ",
      "log transformation and the Greenwood variance (the default of the ",
      "survival R package), and the CI of the median is obtained from the ",
      "pointwise CI of the KM curve (Brookmeyer and Crowley, 1982)."
    )),
    b("paragraph", paste0(
      "Summary data from the KM analysis are presented at ",
      format(s$interval), "-", unit_one, " intervals, including the number ",
      "of events and the number of censored observations in each ",
      "interval, the number at risk, the survival probability, its ",
      "standard error (SE) and its ", level, " CI."
    ))
  )

  if (two_arms) {
    hr <- s$hazard_ratios
    specs <- unique(hr[c("model", "covariates", "strata")])
    described <- vapply(seq_len(nrow(specs)), function(i) {
      details <- c(
        if (specs$covariates[i] != "None") {
          paste("covariates:", specs$covariates[i])
        },
        if (specs$strata[i] != "None") paste("strata:", specs$strata[i])
      )
      if (length(details) == 0L) return(specs$model[i])
      paste0(specs$model[i], " (", paste(details, collapse = "; "), ")")
    }, character(1))
    blocks <- c(blocks, list(
      b("heading", "Hazard Ratios", level = 2L),
      b("paragraph", paste0(
        "The hazard ratio (HR) of each treatment arm versus the control ",
        "arm (", x$arms[1L], ") is estimated using a Cox proportional ",
        "hazards model, with its ", level, " Wald CI and the p-value of ",
        "the Wald test. The following models are reported: ",
        paste(described, collapse = "; "), "."
      ))
    ))
  }

  characteristics <- data.frame(
    "Distribution" = c("Exponential", "Weibull/Gompertz", "Gamma",
                       "Log-logistic/lognormal", "Generalised gamma"),
    "Characteristic of hazard function" = c(
      "Constant over time",
      paste("Can decrease or increase monotonically, but cannot change",
            "direction"),
      paste("Generalisation of exponential distribution. Monotonic, tends",
            "to constant hazard asymptotically"),
      "Can be non-monotonic",
      paste("Can simplify to Weibull, exponential, and lognormal",
            "distributions. Can be non-monotonic")
    ),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  blocks <- c(blocks, list(
    b("heading", "Diagnostics", level = 2L),
    b("heading", "Empirical hazard rates over time", level = 3L),
    b("paragraph", paste0(
      "Smoothed empirical hazard rates are estimated without a parametric ",
      "assumption with the bshazard R package (Rebora, Salim and Reilly, ",
      "2014), which smooths the hazard with B-splines whose smoothness is ",
      "estimated from the data through a Poisson generalised linear mixed ",
      "model, with ", level, " pointwise CI."
    )),
    b("paragraph", paste(
      "The hazard plots can provide insights on the appropriate survival",
      "model, depending on their profile over time as per the table below."
    )),
    b("table", paste("Characteristics of the hazard functions of the",
                     "standard parametric models"),
      table = characteristics),
    b("heading", "Distributional assumptions", level = 3L),
    b("paragraph", paste(
      "The appropriateness of standard parametric models can be assessed",
      "with parametric distribution diagnostic plots, in which a",
      "transformation of the KM estimate S(t) is close to a straight line",
      "when the data follow the distribution: log(-log(S(t))) against log",
      "time for the Weibull distribution (with slope 1 for the exponential",
      "distribution), log((1 - S(t)) / S(t)) against log time for the",
      "log-logistic distribution, the standard normal quantile of",
      "1 - S(t) against log time for the lognormal distribution, and the",
      "logarithm of the smoothed empirical hazard against time for the",
      "Gompertz distribution."
    ))
  ))
  if (two_arms) {
    blocks <- c(blocks, list(
      b("heading", "Proportional hazards assessment", level = 3L),
      b("paragraph", paste(
        "The proportional hazards assumption is examined using the plot of",
        "the scaled Schoenfeld residuals, the Grambsch-Therneau test (with",
        "time transformed by the KM estimate), the plot of log(cumulative",
        "hazard) versus log(time), the plot of the smoothed empirical",
        "hazard versus time, and the quantile-quantile (QQ) plot of the",
        "times of survival percentiles."
      ))
    ))
  }

  fitted <- !vapply(x$fits, is.null, logical(1))
  dependent <- any(m$fit_type == "dependent")
  parametric <- unique(m$model[m$family == "parametric" & fitted])
  if (length(parametric) > 0L) {
    blocks <- c(blocks, list(
      b("heading", "Standard parametric models", level = 2L),
      b("paragraph", paste0(
        "Time to event analysis is conducted using the standard parametric ",
        "models, as described in NICE TSD 14: ", .and_list(parametric),
        ". The models are fitted by maximum likelihood with the flexsurv R ",
        "package to each arm separately",
        if (dependent) {
          paste(" and to both arms together, with the treatment arm as a",
                "covariate on the location parameter of the distribution")
        },
        ". The fits are compared with the Akaike information criterion ",
        "(AIC) and the Bayesian information criterion (BIC, with the number ",
        "of patients as the sample size)."
      ))
    ))
  }

  spline_rows <- m$family == "spline" & fitted
  if (any(spline_rows)) {
    n_knots <- sort(unique(x$knots$n_knots))
    words <- c("one", "two", "three")[n_knots]
    default_knots <- nrow(x$knots) > 0L &&
      all(abs(x$knots$quantile - x$knots$knot / (x$knots$n_knots + 1)) <
            1e-8)
    placement <- c("one knot at the median",
                   "two knots at the 33.3% and 66.7% quantiles",
                   "three knots at the 25%, 50% and 75% quantiles")
    knot_text <- if (default_knots) {
      paste0(
        "The internal knots are placed at quantiles of the log event times ",
        "(the log times of the events, without the censored times): ",
        .and_list(placement[n_knots]), "."
      )
    } else {
      "The internal knots are placed at times chosen by the user."
    }
    scales <- unique(m$scale[spline_rows])
    scale_names <- c(hazard = "hazard", odds = "odds", normal = "normal")
    modelled <- c(hazard = "the log cumulative hazard",
                  odds = "the log cumulative odds",
                  normal = paste("the standard normal quantile of the",
                                 "distribution function"))
    blocks <- c(blocks, list(
      b("heading", "Spline-based extrapolation", level = 2L),
      b("paragraph", paste(
        "Spline-based extrapolation consists of fitting piecewise",
        "polynomials (and specifically restricted cubic splines) to",
        "survival data. Introduced by Royston and Parmar (2002),",
        "spline-based extrapolation allows the shape of the hazard function",
        "to vary over time without making strong parametric assumptions",
        "about its form, resulting in a more adaptable model, whilst",
        "maintaining smoothness."
      )),
      b("paragraph", paste(
        "A cubic spline is a piecewise polynomial function of degree three.",
        "It is continuous and has continuous first and second derivatives",
        "at the points where the pieces join, known as knots. The function",
        "is forced to be smooth by imposing constraints such that the",
        "gradient and the rate of change of the gradient of the function",
        "also agree at the knots. The number and placement of these knots",
        "determines the flexibility of the spline, allowing the model to",
        "adapt to different patterns in the data."
      )),
      b("paragraph", paste(
        "Restricted cubic splines enforce that the spline is linear beyond",
        "the boundary knots, avoiding unrealistic fluctuations at the tails",
        "of the distribution. The boundary knots are placed at the smallest",
        "and largest log event times. Beyond the last event time, the",
        "modelled function is therefore linear in log time, so that the",
        "projection beyond the data follows a Weibull (hazard scale),",
        "log-logistic (odds scale) or lognormal (normal scale) form."
      )),
      b("paragraph", paste(
        "The number of knots and their placement are critical for fitting",
        "cubic splines. Royston and Parmar recommended using 3 to 5 knots,",
        "with their positions typically set at percentiles of the observed",
        "survival times (e.g., at the 25th, 50th, and 75th percentiles).",
        "The choice of knot locations balances the trade-off between",
        "flexibility and overfitting."
      )),
      b("paragraph", paste0(
        "For this analysis, models with ", .and_list(words), " internal ",
        if (length(n_knots) == 1L && n_knots == 1L) "knot are" else
          "knots are",
        " fitted. ", knot_text
      )),
      b("paragraph", paste0(
        "On the ", .and_list(scale_names[scales]),
        if (length(scales) == 1L) " scale, " else " scales, ",
        .and_list(modelled[scales]),
        if (length(scales) == 1L) {
          " is modelled as a restricted cubic spline function of log time"
        } else {
          paste(", respectively, are modelled as restricted cubic spline",
                "functions of log time")
        },
        " (Royston and Parmar, 2002)."
      )),
      b("paragraph", paste0(
        "After specifying each model, the parameters (the spline ",
        "coefficients",
        if (dependent) {
          paste(" and, for the models of both arms, the effect of the",
                "treatment arm on the first spline coefficient, gamma0")
        },
        ") are estimated by maximum likelihood with the flexsurv R ",
        "package. The fits are assessed with the AIC and the BIC, by ",
        "comparing the predicted survival curves with the Kaplan-Meier ",
        "curves, and by comparing the predicted hazards with the smoothed ",
        "empirical hazards."
      ))
    ))
  }

  versions <- vapply(c("SurvFitSummary", "flexsurv", "survival", "bshazard",
                       "ggplot2", "officer"),
                     function(p) as.character(utils::packageVersion(p)),
                     character(1))
  blocks <- c(blocks, list(
    b("heading", "Software", level = 2L),
    b("paragraph", paste0(
      "All analyses were performed using R version ",
      paste(R.version$major, R.version$minor, sep = "."),
      " (https://www.r-project.org/). The following R packages were used: ",
      .and_list(paste(names(versions), unname(versions))), "."
    ))
  ))
  blocks
}
