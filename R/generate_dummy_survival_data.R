#' Generate Dummy Survival Data
#'
#' @description
#' Simulates exponential event times for each treatment arm and exponential
#' dropout times, and returns the observed time and the censoring indicator.
#' The random number generator state of the session is restored afterwards.
#'
#' @param arm Character vector with the names of the treatment arms.
#' @param n Integer vector with the number of subjects in each arm.
#' @param hazards Numeric vector with the hazard rate of the event in each
#'   arm, per month.
#' @param dropout_per_year Proportion of subjects who drop out within one
#'   year, a single number in [0, 1). The dropout hazard per month is
#'   \code{-log(1 - dropout_per_year) / 12}.
#' @param seed Single integer, the random seed.
#'
#' @return A data frame with the columns \code{SUBJID}, \code{ARM},
#'   \code{SURVTIME} (minimum of the event and dropout times, in months) and
#'   \code{CNSR} (1 = censored by dropout, 0 = event).
#'
#' @examples
#' dat <- generate_dummy_survival_data(
#'   arm = c("Control", "Treatment"), n = c(100, 100),
#'   hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 123
#' )
#' head(dat)
#'
#' @export
generate_dummy_survival_data <- function(arm,
                                         n,
                                         hazards,
                                         dropout_per_year,
                                         seed) {
  if (length(arm) != length(n)) {
    stop("Length of 'arm' must match length of 'n'.", call. = FALSE)
  }
  if (length(arm) != length(hazards)) {
    stop("Length of 'arm' must match length of 'hazards'.", call. = FALSE)
  }
  if (!is.numeric(dropout_per_year) || length(dropout_per_year) != 1L ||
      dropout_per_year < 0 || dropout_per_year >= 1) {
    stop("'dropout_per_year' must be a single number in [0, 1).",
         call. = FALSE)
  }
  total_n <- sum(n)
  withr::with_seed(seed, {
    event_time <- unlist(lapply(seq_along(arm), function(i) {
      stats::rexp(n[i], hazards[i])
    }))
    drop_time <- stats::rexp(total_n, rate = -log(1 - dropout_per_year) / 12)
  })
  data.frame(
    SUBJID = paste0("ID", seq_len(total_n)),
    ARM = rep(arm, n),
    SURVTIME = pmin(event_time, drop_time),
    CNSR = as.numeric(drop_time < event_time),
    stringsAsFactors = FALSE
  )
}
