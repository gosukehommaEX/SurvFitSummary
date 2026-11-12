#' Smoothed Empirical Hazard by Treatment Arm
#'
#' @description
#' Estimates the hazard function of each treatment arm without a parametric
#' assumption with \code{bshazard::bshazard()}, which smooths the hazard
#' with B-splines whose smoothness is estimated from the data through a
#' Poisson generalised linear mixed model (Rebora et al., 2014).
#'
#' @param data Data frame returned by \code{\link{processing_dataset}}.
#' @param conf_level Confidence level of the pointwise intervals, a single
#'   number in (0, 1). The default is 0.95.
#' @param nbin Number of time intervals in which the data are split, passed
#'   to \code{bshazard::bshazard()}. The default \code{NULL} uses the
#'   distinct observed times.
#' @param nk Number of knots of the B-splines. The default is 31, as in
#'   \code{bshazard::bshazard()}.
#' @param degree Degree of the B-splines. The default is 1, as in
#'   \code{bshazard::bshazard()}.
#'
#' @return A data frame with the columns \code{ARM}, \code{time},
#'   \code{hazard}, \code{lcl} and \code{ucl}. Rows with a non-finite
#'   hazard are removed.
#'
#' @references
#' Rebora P, Salim A, Reilly M (2014). bshazard: A flexible tool for
#' nonparametric smoothing of the hazard function. The R Journal
#' 6(2):114-122. \doi{10.32614/RJ-2014-028}
#'
#' @examples
#' raw <- generate_dummy_survival_data(
#'   arm = c("Control", "Treatment"), n = c(100, 100),
#'   hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 123
#' )
#' dat <- processing_dataset(raw, column_arm = "ARM",
#'                           column_survtime = "SURVTIME",
#'                           column_cnsr = "CNSR")
#' head(estimate_empirical_hazard(dat))
#'
#' @export
estimate_empirical_hazard <- function(data,
                                      conf_level = 0.95,
                                      nbin = NULL,
                                      nk = 31,
                                      degree = 1) {
  .check_survival_data(data)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  data$ARM <- droplevels(data$ARM)
  rows <- lapply(split(data, data$ARM), function(d) {
    fit <- bshazard::bshazard(survival::Surv(SURVTIME, EVENT) ~ 1, data = d,
                              nbin = nbin, nk = nk, degree = degree,
                              alpha = 1 - conf_level, verbose = FALSE)
    out <- data.frame(
      ARM = as.character(d$ARM[1L]),
      time = fit$time,
      hazard = fit$hazard,
      lcl = fit$lower.ci,
      ucl = fit$upper.ci
    )
    out[is.finite(out$hazard), , drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
