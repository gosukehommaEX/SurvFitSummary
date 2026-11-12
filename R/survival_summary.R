#' Median, Mean and Restricted Mean Survival of the Fitted Models
#'
#' @description
#' Reports, for each fitted model and each arm it predicts, the median
#' survival time, the mean survival time and, when \code{time_horizon} is
#' given, the restricted mean survival time up to the time horizon, as
#' computed by \code{summary.flexsurvreg()}.
#'
#' @details
#' The mean is the integral of the fitted survival function from 0 to
#' infinity; it is infinite for heavy-tailed fits (for example a
#' log-logistic shape of 1 or less) and for a Gompertz model with a negative
#' shape, which leaves a proportion of patients who never have the event.
#' Values that are not finite or cannot be computed are reported as
#' \code{NA}. The confidence limits, when requested, are simulated as in
#' \code{\link{predict_survival}}.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param time_horizon \code{NULL} (the default) or a single positive time
#'   horizon for the restricted mean, in the time unit of the data.
#' @param conf_int Logical; \code{TRUE} to add simulated confidence limits.
#'   The default is \code{FALSE}.
#' @param n_sim Number of simulated parameter vectors for the confidence
#'   limits. The default is 1000.
#' @param seed Seed of the simulation, a single whole number. The default
#'   is 1.
#'
#' @return A data frame with the columns \code{fit_type}, \code{arm} (the
#'   arm predicted), \code{model}, \code{median}, \code{median_lcl},
#'   \code{median_ucl}, \code{mean}, \code{mean_lcl}, \code{mean_ucl},
#'   \code{rmst}, \code{rmst_lcl} and \code{rmst_ucl} (\code{NA} without
#'   \code{time_horizon}; the limits are \code{NA} without
#'   \code{conf_int}).
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' fits <- fit_survival_models(dat, distributions = c("exp", "weibull"),
#'                             spline_scales = NULL)
#' survival_summary(fits, time_horizon = 60)
#'
#' @export
survival_summary <- function(x,
                             time_horizon = NULL,
                             conf_int = FALSE,
                             n_sim = 1000,
                             seed = 1) {
  .check_survfitsummary(x)
  if (!is.null(time_horizon) &&
      (!is.numeric(time_horizon) || length(time_horizon) != 1L ||
       !is.finite(time_horizon) || time_horizon <= 0)) {
    stop("'time_horizon' must be NULL or a single positive number.",
         call. = FALSE)
  }
  if (!is.logical(conf_int) || length(conf_int) != 1L || is.na(conf_int)) {
    stop("'conf_int' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(n_sim) || length(n_sim) != 1L || !is.finite(n_sim) ||
      n_sim < 2 || n_sim != round(n_sim)) {
    stop("'n_sim' must be a single whole number of at least 2.",
         call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != round(seed)) {
    stop("'seed' must be a single whole number.", call. = FALSE)
  }
  # Estimate and limits of one quantity; NA when not finite or failed
  quantity <- function(fit, newdata, type, t = NULL) {
    out <- tryCatch(
      .summary_one(fit, newdata, type = type, t = t, conf_int = conf_int,
                   n_sim = n_sim, seed = seed, conf_level = x$conf_level),
      error = function(e) NULL,
      warning = function(w) NULL
    )
    if (is.null(out)) return(rep(NA_real_, 3L))
    value <- c(out$est[1L], out$lcl[1L], out$ucl[1L])
    value[!is.finite(value)] <- NA_real_
    value
  }
  rows <- list()
  for (i in which(!vapply(x$fits, is.null, logical(1)))) {
    for (pa in .model_arms(x, i)) {
      fit <- x$fits[[i]]
      med <- quantity(fit, pa$newdata, "median")
      avg <- quantity(fit, pa$newdata, "mean")
      res <- if (is.null(time_horizon)) {
        rep(NA_real_, 3L)
      } else {
        quantity(fit, pa$newdata, "rmst", t = time_horizon)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        fit_type = x$models$fit_type[i], arm = pa$arm,
        model = x$models$model[i],
        median = med[1L], median_lcl = med[2L], median_ucl = med[3L],
        mean = avg[1L], mean_lcl = avg[2L], mean_ucl = avg[3L],
        rmst = res[1L], rmst_lcl = res[2L], rmst_ucl = res[3L],
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) {
    empty <- numeric(0)
    return(data.frame(fit_type = character(0), arm = character(0),
                      model = character(0), median = empty,
                      median_lcl = empty, median_ucl = empty, mean = empty,
                      mean_lcl = empty, mean_ucl = empty, rmst = empty,
                      rmst_lcl = empty, rmst_ucl = empty,
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
