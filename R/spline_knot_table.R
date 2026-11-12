#' Placement of the Knots of the Spline Models
#'
#' @description
#' Reports the internal knots of the Royston-Parmar spline models of each
#' fit (each arm of the independent fits and both arms together for the
#' dependent fit) with the Kaplan-Meier estimate at each knot, as in the
#' table of knot placement of the survival analysis report template.
#'
#' @details
#' The knots are the same for the hazard, odds and normal scales. The number
#' at risk and the Kaplan-Meier estimate are those of the data fitted (both
#' arms pooled for the dependent fit) at the knot time; a knot at an
#' observed time counts the patients with that time as at risk and their
#' events in the estimate. The standard error is that of Greenwood.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#'
#' @return A data frame with the columns \code{fit_type}, \code{arm},
#'   \code{n_knots}, \code{knot} (index of the knot), \code{quantile}
#'   (probability of the quantile of the log event times; for knots given
#'   by the user, the proportion of log event times at or below the knot),
#'   \code{log_time}, \code{time}, \code{n_risk}, \code{surv_est} and
#'   \code{surv_se}. It has no rows when no spline model was fitted.
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' fits <- fit_survival_models(dat, distributions = NULL,
#'                             spline_scales = "hazard")
#' spline_knot_table(fits)
#'
#' @export
spline_knot_table <- function(x) {
  .check_survfitsummary(x)
  knots <- x$knots
  if (nrow(knots) == 0L) {
    return(data.frame(fit_type = character(0), arm = character(0),
                      n_knots = integer(0), knot = integer(0),
                      quantile = numeric(0), log_time = numeric(0),
                      time = numeric(0), n_risk = integer(0),
                      surv_est = numeric(0), surv_se = numeric(0),
                      stringsAsFactors = FALSE))
  }
  n_risk <- integer(nrow(knots))
  surv_est <- numeric(nrow(knots))
  surv_se <- numeric(nrow(knots))
  group <- paste(knots$fit_type, knots$arm, sep = "\t")
  for (g in unique(group)) {
    rows <- which(group == g)
    d <- if (knots$fit_type[rows[1L]] == "dependent") {
      x$data
    } else {
      x$data[x$data$ARM == knots$arm[rows[1L]], , drop = FALSE]
    }
    km <- .km_at(d$SURVTIME, d$EVENT, knots$log_time[rows],
                 log_scale = TRUE)
    n_risk[rows] <- km$n_risk
    surv_est[rows] <- km$surv
    surv_se[rows] <- km$se
  }
  out <- data.frame(
    fit_type = knots$fit_type, arm = knots$arm, n_knots = knots$n_knots,
    knot = knots$knot, quantile = knots$quantile,
    log_time = knots$log_time, time = exp(knots$log_time),
    n_risk = n_risk, surv_est = surv_est, surv_se = surv_se,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}
