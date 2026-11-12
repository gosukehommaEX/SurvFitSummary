#' Kaplan-Meier Summary Table at Regular Times
#'
#' @description
#' Reports the Kaplan-Meier estimate of each treatment arm at the times
#' \code{0, interval, 2 * interval, ...} up to the largest observed time of
#' the arm, with the number at risk, the numbers of events and censored
#' observations since the previous time, the standard error (Greenwood) and
#' the pointwise confidence interval of \code{survival::survfit()}.
#'
#' @param data Data frame returned by \code{\link{processing_dataset}}.
#' @param interval Positive number, the spacing of the reporting times in
#'   the time unit of the data.
#' @param conf_level Confidence level, a single number in (0, 1). The
#'   default is 0.95.
#' @param conf_type Transformation for the pointwise confidence intervals:
#'   \code{"log"} (the default of \code{survival::survfit()}),
#'   \code{"log-log"} or \code{"plain"}.
#'
#' @return A data frame with the columns \code{ARM}, \code{time},
#'   \code{n_risk} (number with an observed time at or after \code{time}),
#'   \code{n_event} and \code{n_censor} (numbers of events and censored
#'   observations after the previous reporting time and up to \code{time}),
#'   \code{surv} (Kaplan-Meier estimate), \code{se} (standard error of
#'   \code{surv}), \code{lcl} and \code{ucl} (confidence limits; both equal
#'   to 1 while \code{surv} is 1).
#'
#' @examples
#' raw <- generate_dummy_survival_data(
#'   arm = c("Control", "Treatment"), n = c(100, 100),
#'   hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 123
#' )
#' dat <- processing_dataset(raw, column_arm = "ARM",
#'                           column_survtime = "SURVTIME",
#'                           column_cnsr = "CNSR")
#' km_table(dat, interval = 6)
#'
#' @export
km_table <- function(data,
                     interval,
                     conf_level = 0.95,
                     conf_type = c("log", "log-log", "plain")) {
  .check_survival_data(data)
  conf_type <- match.arg(conf_type)
  if (!is.numeric(interval) || length(interval) != 1L || interval <= 0) {
    stop("'interval' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  data$ARM <- droplevels(data$ARM)
  rows <- lapply(split(data, data$ARM), function(d) {
    fit <- survival::survfit(survival::Surv(SURVTIME, EVENT) ~ 1, data = d,
                             conf.int = conf_level, conf.type = conf_type)
    # Reporting times; the small allowance avoids losing the last multiple
    # of 'interval' to rounding
    grid <- interval * seq(0, floor(max(d$SURVTIME) / interval + 1e-10))
    # Right-continuous step function of the Kaplan-Meier estimate
    idx <- findInterval(grid, fit$time) + 1L
    surv <- c(1, fit$surv)[idx]
    se_log <- c(0, fit$std.err)[idx]
    lcl <- c(1, fit$lower)[idx]
    ucl <- c(1, fit$upper)[idx]
    lcl[surv == 1] <- 1
    ucl[surv == 1] <- 1
    se <- surv * se_log
    se[!is.finite(se)] <- NA_real_
    previous <- c(-Inf, grid[-length(grid)])
    in_window <- function(k, status) {
      sum(d$SURVTIME > previous[k] & d$SURVTIME <= grid[k] &
            d$EVENT == status)
    }
    k <- seq_along(grid)
    data.frame(
      ARM = as.character(d$ARM[1L]),
      time = grid,
      n_risk = vapply(grid, function(g) sum(d$SURVTIME >= g), integer(1)),
      n_event = vapply(k, in_window, integer(1), status = 1L),
      n_censor = vapply(k, in_window, integer(1), status = 0L),
      surv = surv,
      se = se,
      lcl = lcl,
      ucl = ucl
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
