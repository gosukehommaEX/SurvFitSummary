#' Summarise Survival Data by Treatment Arm
#'
#' @description
#' Counts the patients, events and censored observations in each treatment
#' arm and estimates the median survival time with its confidence interval
#' from the Kaplan-Meier estimate (\code{survival::survfit()}). The median is
#' the first time at which the Kaplan-Meier curve is at or below 0.5, and its
#' confidence limits are the first times at which the pointwise lower and
#' upper confidence limits of the curve are at or below 0.5.
#'
#' @param data Data frame returned by \code{\link{processing_dataset}}.
#' @param conf_level Confidence level, a single number in (0, 1). The
#'   default is 0.95.
#' @param conf_type Transformation for the pointwise confidence intervals of
#'   the survival probability: \code{"log"} (the default of
#'   \code{survival::survfit()}), \code{"log-log"} or \code{"plain"}.
#'
#' @return A data frame with one row per treatment arm and, when there are
#'   two or more arms, a last row \code{"Overall"} for the pooled data. The
#'   columns are \code{ARM}, \code{n}, \code{events}, \code{censored},
#'   \code{median}, \code{median_lcl}, \code{median_ucl} (\code{NA} when the
#'   curve or the confidence limit does not reach 0.5) and
#'   \code{max_followup} (largest observed time).
#'
#' @examples
#' raw <- generate_dummy_survival_data(
#'   arm = c("Control", "Treatment"), n = c(100, 100),
#'   hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 123
#' )
#' dat <- processing_dataset(raw, column_arm = "ARM",
#'                           column_survtime = "SURVTIME",
#'                           column_cnsr = "CNSR")
#' summarise_data(dat)
#'
#' @export
summarise_data <- function(data,
                           conf_level = 0.95,
                           conf_type = c("log", "log-log", "plain")) {
  .check_survival_data(data)
  conf_type <- match.arg(conf_type)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  data$ARM <- droplevels(data$ARM)
  groups <- split(data, data$ARM)
  if (length(groups) > 1L) {
    groups[["Overall"]] <- data
  }
  rows <- lapply(names(groups), function(name) {
    d <- groups[[name]]
    fit <- survival::survfit(survival::Surv(SURVTIME, EVENT) ~ 1, data = d,
                             conf.int = conf_level, conf.type = conf_type)
    q <- stats::quantile(fit, probs = 0.5, conf.int = TRUE)
    data.frame(
      ARM = name,
      n = nrow(d),
      events = sum(d$EVENT),
      censored = sum(1L - d$EVENT),
      median = unname(q$quantile),
      median_lcl = unname(q$lower),
      median_ucl = unname(q$upper),
      max_followup = max(d$SURVTIME)
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
