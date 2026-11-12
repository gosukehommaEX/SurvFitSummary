#' Summary of One Fitted Model for One Arm
#'
#' @description
#' Internal helper that calls \code{summary.flexsurvreg()} for one model and
#' one arm and returns the estimate and confidence limits in the order of
#' \code{t}. The confidence limits are simulated by flexsurv from the
#' asymptotic normal distribution of the estimates, with the random number
#' generator seeded by \code{seed} and restored afterwards.
#'
#' @param fit A flexsurvreg object.
#' @param newdata \code{NULL} or a data frame with the factor \code{ARM}.
#' @param type Quantity passed to \code{summary.flexsurvreg()}, for example
#'   \code{"survival"}, \code{"hazard"}, \code{"median"}.
#' @param t \code{NULL} or the times.
#' @param conf_int Logical; \code{TRUE} to simulate confidence limits.
#' @param n_sim Number of simulated parameter vectors.
#' @param seed Seed of the simulation.
#' @param conf_level Confidence level.
#'
#' @return A list with \code{est}, \code{lcl} and \code{ucl} (\code{NA}
#'   when \code{conf_int} is \code{FALSE}).
#'
#' @keywords internal
#' @noRd
.summary_one <- function(fit, newdata, type, t = NULL, conf_int = FALSE,
                         n_sim = 1000, seed = 1, conf_level = 0.95) {
  args <- list(object = fit, type = type, ci = conf_int, B = n_sim,
               cl = conf_level)
  if (!is.null(t)) args$t <- t
  if (!is.null(newdata)) args$newdata <- newdata
  out <- withr::with_seed(seed, do.call(summary, args))[[1L]]
  idx <- if (!is.null(t) && "time" %in% names(out)) {
    match(t, out[["time"]])
  } else {
    seq_len(nrow(out))
  }
  column <- function(name) {
    if (name %in% names(out)) {
      unname(out[[name]][idx])
    } else {
      rep(NA_real_, length(idx))
    }
  }
  list(est = column("est"), lcl = column("lcl"), ucl = column("ucl"))
}
