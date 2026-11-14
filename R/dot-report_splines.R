#' Spline Models in the Word Report
#'
#' @description
#' Internal helper: the section of the report on spline-based
#' extrapolation. For each arm fitted separately and for both arms fitted
#' together: the placement of the knots, the AIC, BIC, median, mean and
#' predicted survival at the landmark times of each model (one table per
#' arm for the models of both arms, so that the table fits the page), and
#' the smoothed empirical hazard and the Kaplan-Meier curve overlaid with
#' the models. As in the example figures of the template of the report,
#' the plots also show the standard parametric model with the lowest AIC of
#' the same fit (dashed), the plots of the hazard show the confidence
#' interval of the smoothed hazard when the report has figures with
#' confidence intervals, and the plots of the Kaplan-Meier curve mark the
#' range of the predicted survival at 5 and 10 years (within the time
#' horizon).
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param s List returned by \code{.report_settings()}.
#' @param conf_int Logical vector: versions of the figures without
#'   (\code{FALSE}) and/or with (\code{TRUE}) confidence intervals of the
#'   report.
#' @param n_sim,seed Simulation of the confidence intervals of the models
#'   (passed to the plotting functions).
#'
#' @return A list of elements made by \code{.report_block()}; an empty list
#'   when no spline model was fitted.
#'
#' @keywords internal
#' @noRd
.report_splines <- function(x, s, conf_int, n_sim, seed) {
  groups <- .report_groups(x, "spline")
  if (length(groups) == 0L) return(list())
  b <- .report_block
  level <- paste0(format(100 * x$conf_level), "%")
  gof <- gof_table(x)
  parametric <- unique(x$models$model[x$models$family == "parametric"])
  # Standard parametric model with the lowest AIC of a fit, shown with the
  # spline models for comparison; NULL without parametric models
  best_of <- function(g) {
    rows <- gof[gof$fit_type == g$fit_type & gof$arm == g$arm &
                  gof$model %in% parametric, , drop = FALSE]
    if (nrow(rows) == 0L) NULL else rows$model[which.min(rows$AIC)]
  }
  bests <- lapply(groups, best_of)
  # The confidence interval of the smoothed hazard when the report has
  # figures with confidence intervals (with one model, the plot would also
  # show the interval of the model)
  band <- TRUE %in% conf_int
  bands <- vapply(seq_along(groups), function(i) {
    band && length(c(groups[[i]]$models, bests[[i]])) > 1L
  }, logical(1))
  band_text <- function(band) {
    if (band) paste0(" (with its ", level, " confidence interval)") else ""
  }
  # Times of the ranges of the predicted survival in the plots of the
  # Kaplan-Meier curves: 5 and 10 years, within the time horizon
  years <- .landmark_times(s$unit)
  years <- years[years$label %in% c("5 years", "10 years") &
                   years$time <= s$time_horizon, , drop = FALSE]
  fit_args <- function(fun, g, arm, best, ci, range_times = NULL) {
    # One arm with its models: the legend on the right, as agreed for the
    # report; the standard parametric model dashed
    args <- list(x = x, models = c(g$models, best), arms = arm,
                 fit_type = g$fit_type, conf_int = ci,
                 time_horizon = s$time_horizon, n_sim = n_sim, seed = seed,
                 legend_position = "right")
    if (!is.null(best)) args$linetypes <- stats::setNames("dashed", best)
    if (length(range_times) > 0L) args$range_times <- range_times
    list(fun = fun, args = args)
  }
  knots <- spline_knot_table(x)
  ss <- survival_summary(x)
  pred <- predict_survival(x, s$landmarks$time)
  intro <- paste0(
    "For each arm, and for both arms fitted together, the placement of the ",
    "knots, the fit statistics and the predicted survival of the spline ",
    "models are tabulated, and the smoothed empirical hazards",
    if (any(bands)) paste0(" (with their ", level, " confidence intervals)"),
    " and the Kaplan-Meier curves are overlaid with the spline models, ",
    "extrapolated to ", format(s$time_horizon), " ", s$unit, ".",
    if (!all(vapply(bests, is.null, logical(1)))) {
      paste(" For comparison, the plots also show the standard parametric",
            "model with the lowest AIC (dashed line).")
    },
    if (nrow(years) > 0L) {
      paste0(" In the plots of the Kaplan-Meier curves, a vertical line at ",
             .and_list(years$label), " shows the range of the survival ",
             "predicted by the models, labelled in percent.")
    }
  )
  blocks <- list(b("heading", "Spline-based extrapolation", level = 2L),
                 b("paragraph", intro))
  for (i in seq_along(groups)) {
    g <- groups[[i]]
    best <- bests[[i]]
    by_both <- g$fit_type == "dependent"
    blocks <- c(blocks, list(b("heading", g$arm, level = 3L)))
    k <- knots[knots$fit_type == g$fit_type & knots$arm == g$arm, ,
               drop = FALSE]
    knot_table <- data.frame(
      n_knots = as.character(k$n_knots),
      knot = as.character(k$knot),
      quantile = .format_number(k$quantile, 3L),
      log_time = .format_number(k$log_time, 3L),
      time = .format_number(k$time, 2L),
      n_risk = as.character(k$n_risk),
      surv_est = .format_number(k$surv_est, 3L),
      surv_se = .format_number(k$surv_se, 3L),
      stringsAsFactors = FALSE
    )
    fit_tables <- list()
    for (a in g$arms) {
      rows <- list()
      for (m in g$models) {
        fit_row <- gof$fit_type == g$fit_type & gof$arm == g$arm &
          gof$model == m
        sum_row <- ss$fit_type == g$fit_type & ss$arm == a & ss$model == m
        pred_rows <- pred$fit_type == g$fit_type & pred$arm == a &
          pred$model == m
        survival <- pred$est[pred_rows][match(s$landmarks$time,
                                              pred$time[pred_rows])]
        rows[[length(rows) + 1L]] <- c(
          m,
          .format_number(gof$AIC[fit_row], 1L),
          .format_number(gof$BIC[fit_row], 1L),
          .format_number(ss$median[sum_row], 1L),
          .format_number(ss$mean[sum_row], 1L),
          .format_number(survival, 3L)
        )
      }
      fit_table <- as.data.frame(do.call(rbind, rows),
                                 stringsAsFactors = FALSE)
      names(fit_table) <- c("Model", "AIC", "BIC",
                            paste0("Median (", s$unit, ")"),
                            paste0("Mean (", s$unit, ")"),
                            paste("Survival at", s$landmarks$label))
      title <- if (by_both) paste0(g$arm, ", ", a) else g$arm
      fit_tables <- c(fit_tables, list(
        b("table", paste0("Tests of fit and predicted survival ",
                          "probability (", title, ")"), table = fit_table)
      ))
    }
    blocks <- c(blocks, list(
      b("table", paste0("Placement of knots (", g$arm, ")"),
        table = knot_table),
      b("paragraph", paste0(
        "Abbreviations: n_knots, number of internal knots; knot, index of ",
        "the knot; quantile, quantile of the log event times; log_time, ",
        "log time of the knot; time, time of the knot (", s$unit, "); ",
        "n_risk, number at risk; surv_est, estimated Kaplan-Meier survival ",
        "probability; surv_se, standard error of the estimated survival ",
        "probability."
      ))
    ))
    blocks <- c(blocks, fit_tables, list(
      b("paragraph", paste(
        "Abbreviations: AIC, Akaike information criterion; BIC, Bayesian",
        "information criterion; NE, not estimable."
      ))
    ))
    fitted_to <- if (by_both) " fitted to both arms" else ""
    with_best <- if (is.null(best)) {
      ""
    } else {
      paste0(" and of the ", best, " model (dashed)")
    }
    for (a in g$arms) {
      blocks <- c(blocks, list(
        b("figure", paste0(
          "Smoothed empirical hazard of ", a, band_text(bands[[i]]),
          " and the hazards of the spline models", fitted_to, with_best
        ), plot = fit_args(plot_hazard_fits, g, a, best, bands[[i]]),
        height = 5.5),
        b("figure", paste0(
          "Kaplan-Meier curve of ", a, " and the survival curves of the ",
          "spline models", fitted_to, with_best,
          if (nrow(years) > 0L) {
            paste0(", with the range of the predicted survival at ",
                   .and_list(years$label))
          }
        ), plot = fit_args(plot_survival_fits, g, a, best, FALSE,
                           years$time),
        height = 6.5)
      ))
    }
  }
  blocks
}
