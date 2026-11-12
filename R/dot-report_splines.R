#' Spline Models in the Word Report
#'
#' @description
#' Internal helper: the section of the report template on spline-based
#' extrapolation. For each arm fitted separately and for both arms fitted
#' together: the placement of the knots, the AIC, BIC, median, mean and
#' predicted survival at the landmark times of each model (one table per
#' arm for the models of both arms, so that the table fits the page), and
#' the smoothed empirical hazard and the Kaplan-Meier curve overlaid with
#' the models.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param s List returned by \code{.report_settings()}.
#' @param n_sim,seed Simulation of the confidence intervals of the models
#'   (passed to the plotting functions).
#'
#' @return A list of elements made by \code{.report_block()}; an empty list
#'   when no spline model was fitted.
#'
#' @keywords internal
#' @noRd
.report_splines <- function(x, s, n_sim, seed) {
  groups <- .report_groups(x, "spline")
  if (length(groups) == 0L) return(list())
  b <- .report_block
  fit_args <- function(fun, g, arm) {
    list(fun = fun,
         args = list(x = x, models = g$models, arms = arm,
                     fit_type = g$fit_type, conf_int = FALSE,
                     time_horizon = s$time_horizon, n_sim = n_sim,
                     seed = seed))
  }
  gof <- gof_table(x)
  knots <- spline_knot_table(x)
  ss <- survival_summary(x)
  pred <- predict_survival(x, s$landmarks$time)
  blocks <- list(b("heading", "Spline-based extrapolation", level = 2L))
  for (g in groups) {
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
    for (a in g$arms) {
      blocks <- c(blocks, list(
        b("figure", paste0(
          "Smoothed empirical hazard of ", a, " and the hazards of the ",
          "spline models", fitted_to
        ), plot = fit_args(plot_hazard_fits, g, a), height = 5.5),
        b("figure", paste0(
          "Kaplan-Meier curve of ", a, " and the survival curves of the ",
          "spline models", fitted_to
        ), plot = fit_args(plot_survival_fits, g, a), height = 6.5)
      ))
    }
  }
  blocks
}
