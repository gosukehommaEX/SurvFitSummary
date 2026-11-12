#' Estimated Parameters and Covariance Matrix of a Fitted Model
#'
#' @description
#' Internal helper that returns the estimated parameters of a flexsurvreg
#' object on the transformed (real-line) scale of flexsurv, their covariance
#' matrix and which of them are on the log scale. The names are taken from
#' the row names of \code{res.t}, because \code{coefficients} loses its name
#' when the model has a single parameter. The shape of the Rayleigh model,
#' fixed at 2, is not a parameter and is left out; the covariance matrix is
#' taken for the estimated parameters whether or not flexsurv includes a row
#' for the fixed one.
#'
#' @param fit A flexsurvreg object.
#' @param distribution Distribution code of \code{.distribution_info()}, or
#'   \code{NA} for a spline model (whose parameters are not transformed).
#'
#' @return A list with \code{est} (named vector on the transformed scale),
#'   \code{cov} (covariance matrix on the transformed scale, \code{NA} when
#'   not available) and \code{log_scale} (logical vector, \code{TRUE} for
#'   the parameters estimated on the log scale).
#'
#' @keywords internal
#' @noRd
.model_parameters <- function(fit, distribution) {
  est_all <- fit$res.t[, 1L]
  names(est_all) <- rownames(fit$res.t)
  fixed <- if (!is.na(distribution) && distribution == "rayleigh") {
    1L
  } else {
    integer(0)
  }
  free <- setdiff(seq_along(est_all), fixed)
  est <- est_all[free]
  cov <- fit$cov
  if (is.matrix(cov) && length(fixed) > 0L &&
      nrow(cov) == length(est_all)) {
    cov <- cov[free, free, drop = FALSE]
  }
  if (!is.matrix(cov) || nrow(cov) != length(est) ||
      ncol(cov) != length(est)) {
    cov <- matrix(NA_real_, length(est), length(est))
  }
  dimnames(cov) <- list(names(est), names(est))
  log_names <- if (is.na(distribution)) {
    character(0)
  } else {
    info <- .distribution_info()
    strsplit(info$log_parameters[info$code == distribution], ",")[[1L]]
  }
  list(est = est, cov = cov, log_scale = names(est) %in% log_names)
}
