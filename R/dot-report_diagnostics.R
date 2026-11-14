#' Diagnostics in the Word Report
#'
#' @description
#' Internal helper: the Diagnostics part of the Results section of the
#' report, with the smoothed empirical hazards (all arms and each arm, with
#' and/or without confidence intervals), the distribution diagnostic plots,
#' and, with two or more arms, the scaled Schoenfeld residuals, the
#' Grambsch-Therneau test and the QQ plot, each with a text on how to read
#' it. A straight line through the origin in the QQ plot indicates a
#' constant ratio of the survival times (accelerated failure time), which
#' is not the same as proportional hazards.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param s List returned by \code{.report_settings()}.
#' @param conf_int Logical vector: versions of the figures without
#'   (\code{FALSE}) and/or with (\code{TRUE}) confidence intervals.
#'
#' @return A list of elements made by \code{.report_block()}.
#'
#' @keywords internal
#' @noRd
.report_diagnostics <- function(x, s, conf_int) {
  b <- .report_block
  level <- paste0(format(100 * x$conf_level), "%")
  ci_text <- function(ci) {
    if (ci) paste0(", with ", level, " confidence intervals") else ""
  }
  two_arms <- length(x$arms) >= 2L
  hazard_args <- function(arms, ci) {
    list(fun = plot_empirical_hazard,
         args = list(data = x$data, arms = arms, conf_int = ci,
                     conf_level = x$conf_level, time_unit = x$time_unit))
  }
  blocks <- list(
    b("heading", "Diagnostics", level = 2L),
    b("heading", "Empirical hazard rates over time", level = 3L)
  )
  # As the Kaplan-Meier curves: a set of figures without confidence
  # intervals (all arms, then each arm), then a set with them
  for (ci in conf_int) {
    if (two_arms) {
      blocks <- c(blocks, list(b(
        "figure", paste0("Smoothed empirical hazards of all arms",
                         ci_text(ci)),
        plot = hazard_args(NULL, ci), height = 5.5
      )))
    }
    for (a in x$arms) {
      blocks <- c(blocks, list(b(
        "figure", paste0("Smoothed empirical hazard of ", a, ci_text(ci)),
        plot = hazard_args(a, ci), height = 5.5
      )))
    }
  }
  blocks <- c(blocks, list(
    b("heading", "Distributional assumptions", level = 3L),
    b("figure", "Diagnostic plots of the standard parametric distributions",
      plot = list(fun = plot_distribution_diagnostics,
                  args = list(data = x$data, time_unit = x$time_unit,
                              legend_position = "right")),
      height = 8),
    b("paragraph", paste(
      "In each plot, a transformation of the Kaplan-Meier estimate is",
      "plotted against log time (Weibull, log-logistic and lognormal",
      "plots), or the logarithm of the smoothed empirical hazard against",
      "time (Gompertz plot), with a dashed least-squares line for each arm.",
      "A curve close to a straight line supports the distribution of the",
      "plot; in the Weibull plot, a straight line with slope 1 supports the",
      "exponential distribution. Systematic curvature indicates that the",
      "distribution may not describe the data well. Parallel lines of the",
      "arms indicate proportional hazards in the Weibull plot and",
      "proportional odds in the log-logistic plot."
    ))
  ))
  if (two_arms) {
    ph <- test_proportional_hazards(x$data)
    p_text <- .format_p(ph$p_value)
    p_text <- ifelse(substr(p_text, 1L, 1L) == "<",
                     paste("p <", substring(p_text, 2L)),
                     paste("p =", p_text))
    blocks <- c(blocks, list(
      b("heading", "Proportional hazards assessment", level = 3L),
      b("figure", "Scaled Schoenfeld residuals of the treatment effect",
        plot = list(fun = plot_schoenfeld,
                    args = list(data = x$data, conf_level = x$conf_level,
                                time_unit = x$time_unit)),
        height = 6),
      b("paragraph", paste0(
        "The plot shows the scaled Schoenfeld residuals of the treatment ",
        "effect at the event times (points), the log hazard ratio (HR) of ",
        "the Cox model (pink line), a natural spline regression of the ",
        "residuals on time (black line) and its ", level, " confidence ",
        "interval (dotted lines). A regression line close to horizontal ",
        "indicates a treatment effect that does not change over time, ",
        "consistent with proportional hazards."
      )),
      b("paragraph", paste0(
        "The Grambsch-Therneau test of the proportional hazards assumption ",
        "gives ", .and_list(paste0(p_text, " (", ph$comparison, ")")),
        "; a small p-value indicates that the hazard ratio changes over ",
        "time."
      )),
      b("paragraph", paste(
        "In the plot of log(cumulative hazard) versus log(time) (the",
        "Weibull plot of the distribution diagnostics), parallel curves of",
        "the arms indicate proportional hazards. In the plot of the",
        "smoothed empirical hazards of all arms, a constant ratio of the",
        "hazards over time indicates proportional hazards."
      )),
      b("figure", "Quantile-quantile plot of the survival times",
        plot = list(fun = plot_qq_survival,
                    args = list(data = x$data, time_unit = x$time_unit)),
        # The square plot is drawn 7 inches wide and shown at the scale of
        # the other figures (6.2 / 9), so that its text has the same size
        width = 7, height = 7.5, display_width = 7 * 6.2 / 9),
      b("paragraph", paste(
        "In the QQ plot, the quantiles of the survival times of each",
        "treatment arm are plotted against those of the control arm, with",
        "a least-squares line. When the points follow a straight line",
        "through the origin, the survival times of the arms differ by a",
        "constant factor (accelerated failure time), estimated by the slope",
        "of the line. Proportional hazards and accelerated failure time hold",
        "together only for the Weibull distribution, so the QQ plot",
        "complements, but does not test, the proportional hazards",
        "assumption."
      ))
    ))
  }
  blocks
}
