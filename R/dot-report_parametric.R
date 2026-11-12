#' Standard Parametric Models in the Word Report
#'
#' @description
#' Internal helper: the section of the report on the standard
#' parametric models. For each arm fitted separately and for both arms
#' fitted together: the smoothed empirical hazard and the Kaplan-Meier
#' curve overlaid with all the models, the table of AIC and BIC in
#' ascending order of AIC, and for each model its hazard and survival
#' plots, with and/or without confidence intervals (both arms in one plot
#' for the models of both arms).
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param s List returned by \code{.report_settings()}.
#' @param conf_int Logical vector: versions of the figures without
#'   (\code{FALSE}) and/or with (\code{TRUE}) confidence intervals.
#' @param n_sim,seed Simulation of the confidence intervals of the models.
#'
#' @return A list of elements made by \code{.report_block()}; an empty list
#'   when no parametric model was fitted.
#'
#' @keywords internal
#' @noRd
.report_parametric <- function(x, s, conf_int, n_sim, seed) {
  groups <- .report_groups(x, "parametric")
  if (length(groups) == 0L) return(list())
  b <- .report_block
  level <- paste0(format(100 * x$conf_level), "%")
  ci_text <- function(ci) {
    if (ci) paste0(", with ", level, " confidence intervals") else ""
  }
  fit_args <- function(fun, g, models, arms, ci) {
    list(fun = fun,
         args = list(x = x, models = models, arms = arms,
                     fit_type = g$fit_type, conf_int = ci,
                     time_horizon = s$time_horizon, n_sim = n_sim,
                     seed = seed))
  }
  gof <- gof_table(x)
  blocks <- list(
    b("heading", "Standard parametric models", level = 2L),
    b("paragraph", paste0(
      "For each arm, and for both arms fitted together, the smoothed ",
      "empirical hazards and the Kaplan-Meier curves are overlaid with the ",
      "fitted models, extrapolated to ", format(s$time_horizon), " ",
      s$unit, ", and the goodness of fit statistics are listed in ",
      "ascending order of AIC."
    ))
  )
  for (g in groups) {
    by_both <- g$fit_type == "dependent"
    fitted_to <- if (by_both) " fitted to both arms" else ""
    blocks <- c(blocks, list(b("heading", g$arm, level = 3L)))
    for (a in g$arms) {
      blocks <- c(blocks, list(
        b("figure", paste0(
          "Smoothed empirical hazard of ", a, " and the hazards of the ",
          "standard parametric models", fitted_to
        ), plot = fit_args(plot_hazard_fits, g, g$models, a, FALSE),
        height = 5.5),
        b("figure", paste0(
          "Kaplan-Meier curve of ", a, " and the survival curves of the ",
          "standard parametric models", fitted_to
        ), plot = fit_args(plot_survival_fits, g, g$models, a, FALSE),
        height = 6.5)
      ))
    }
    rows <- gof[gof$fit_type == g$fit_type & gof$arm == g$arm &
                  gof$model %in% g$models, , drop = FALSE]
    gof_text <- data.frame(
      "Distribution" = rows$model,
      "AIC" = .format_number(rows$AIC, 1L),
      "BIC" = .format_number(rows$BIC, 1L),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    blocks <- c(blocks, list(
      b("table", paste0(
        "Goodness of fit statistics of the standard parametric models (",
        g$arm, "), in ascending order of AIC"
      ), table = gof_text),
      b("paragraph", paste(
        "Abbreviations: AIC, Akaike information criterion; BIC, Bayesian",
        "information criterion."
      ))
    ))
    arms_text <- .and_list(g$arms)
    for (d in g$models) {
      blocks <- c(blocks, list(
        b("heading", paste("Distribution:", d), level = 3L),
        b("heading", "Smoothed empirical hazard plots", level = 4L)
      ))
      for (ci in conf_int) {
        blocks <- c(blocks, list(b("figure", paste0(
          "Smoothed empirical hazard of ", arms_text, " and the hazard of ",
          "the ", d, " model", fitted_to, ci_text(ci)
        ), plot = fit_args(plot_hazard_fits, g, d, g$arms, ci),
        height = 5.5)))
      }
      blocks <- c(blocks, list(b("heading", "Kaplan-Meier curves",
                                 level = 4L)))
      for (ci in conf_int) {
        blocks <- c(blocks, list(b("figure", paste0(
          "Kaplan-Meier curve of ", arms_text, " and the survival curve of ",
          "the ", d, " model", fitted_to, ci_text(ci)
        ), plot = fit_args(plot_survival_fits, g, d, g$arms, ci),
        height = 6.5)))
      }
    }
  }
  blocks
}
