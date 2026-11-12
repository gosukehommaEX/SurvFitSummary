#' Methods Section of the Word Report
#'
#' @description
#' Internal helper: the Methods section of the report, with the methods,
#' settings and software versions of the analysis filled in: Kaplan-Meier
#' analysis, hazard ratios, diagnostics (with a table of the shapes of the
#' hazard functions of the standard parametric models), standard parametric
#' models, spline-based models and software.
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
      "Survival probabilities and the median time to event are estimated ",
      "with the Kaplan-Meier (KM) method. ", level, " confidence intervals ",
      "(CI) of the survival probabilities use the log transformation and ",
      "the Greenwood variance (the default of the survival R package); the ",
      "CI of the median is obtained from the pointwise CI of the KM curve ",
      "(Brookmeyer and Crowley, 1982)."
    )),
    b("paragraph", paste0(
      "The KM estimates are tabulated every ", format(s$interval), " ",
      s$unit, ": the number at risk, the numbers of events and of censored ",
      "observations since the previous time, and the survival probability ",
      "with its standard error (SE) and ", level, " CI."
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
        "Hazard ratios (HR) of each treatment arm versus the control arm (",
        x$arms[1L], ") come from Cox proportional hazards models, with ",
        level, " Wald CI and Wald test p-values. Models: ",
        paste(described, collapse = "; "), "."
      ))
    ))
  }

  shapes <- data.frame(
    "Distribution" = c("Exponential", "Weibull", "Gompertz", "Gamma",
                       "Log-logistic", "Lognormal", "Generalized gamma"),
    "Shape of the hazard function" = c(
      "Constant",
      "Constant, or increasing or decreasing over the whole time",
      "Constant, or increasing or decreasing exponentially",
      paste("Constant, or increasing or decreasing towards a constant",
            "value"),
      "Decreasing, or first increasing and then decreasing",
      "First increasing and then decreasing",
      paste("Increasing, decreasing, first increasing and then decreasing,",
            "or the reverse; includes the Weibull, exponential, lognormal",
            "and gamma distributions")
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
      "The shape of the empirical hazard over time helps to choose among",
      "the standard parametric models, whose hazard functions can take the",
      "shapes listed below."
    )),
    b("table", paste("Shapes of the hazard functions of the standard",
                     "parametric models"), table = shapes),
    b("heading", "Distributional assumptions", level = 3L),
    b("paragraph", paste(
      "Diagnostic plots show whether a distribution is compatible with the",
      "data: a transformation of the KM estimate S(t) is close to a",
      "straight line when the data follow the distribution. The plots are",
      "log(-log(S(t))) against log time for the Weibull distribution (with",
      "slope 1 for the exponential distribution), log((1 - S(t)) / S(t))",
      "against log time for the log-logistic distribution, the standard",
      "normal quantile of 1 - S(t) against log time for the lognormal",
      "distribution, and the logarithm of the smoothed empirical hazard",
      "against time for the Gompertz distribution."
    ))
  ))
  if (two_arms) {
    blocks <- c(blocks, list(
      b("heading", "Proportional hazards assessment", level = 3L),
      b("paragraph", paste(
        "Proportional hazards are assessed with the scaled Schoenfeld",
        "residuals of the treatment effect and the Grambsch-Therneau test",
        "(with time transformed by the KM estimate), the plot of",
        "log(cumulative hazard) against log(time), the smoothed empirical",
        "hazards of the arms over time, and the quantile-quantile (QQ) plot",
        "of the survival times."
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
        "The parametric models fitted are: ", .and_list(parametric),
        " (NICE DSU Technical Support Document 14; Latimer, 2011). They are ",
        "fitted by maximum likelihood with the flexsurv R package to each ",
        "arm separately",
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
        "Spline-based models (Royston and Parmar, 2002) describe a",
        "transformation of the survival function as a smooth function of",
        "log time, so that the shape of the hazard is estimated from the",
        "data instead of being fixed by a parametric family."
      )),
      b("paragraph", paste0(
        "On the ", .and_list(scale_names[scales]),
        if (length(scales) == 1L) " scale, " else " scales, ",
        .and_list(modelled[scales]),
        if (length(scales) == 1L) {
          " is modelled as a restricted cubic spline function of log time."
        } else {
          paste(", respectively, are modelled as restricted cubic spline",
                "functions of log time.")
        }
      )),
      b("paragraph", paste(
        "A restricted cubic spline consists of cubic polynomials between",
        "the knots, joined so that the function and its first and second",
        "derivatives are continuous at the knots, and is linear beyond the",
        "boundary knots, which are placed at the smallest and largest log",
        "event times. More knots give a more flexible function, at the risk",
        "of overfitting. Because the function is linear in log time after",
        "the last event time, the extrapolation follows a Weibull (hazard",
        "scale), log-logistic (odds scale) or lognormal (normal scale) form."
      )),
      b("paragraph", paste0(
        "Models with ", .and_list(words), " internal ",
        if (length(n_knots) == 1L && n_knots == 1L) "knot are" else
          "knots are",
        " fitted. ", knot_text
      )),
      b("paragraph", paste0(
        "The spline coefficients",
        if (dependent) {
          paste(" (and, for the models of both arms, the effect of the",
                "treatment arm on the first coefficient, gamma0)")
        },
        " are estimated by maximum likelihood with the flexsurv R package. ",
        "The models are compared with the AIC and the BIC, and their ",
        "survival and hazard curves with the KM curves and the smoothed ",
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
