#' Parametric Distributions Available in the Package
#'
#' @description
#' Internal table of the parametric distributions: the code used in
#' \code{fit_survival_models(distributions = )}, the print-ready label used
#' in tables and reports, the distribution passed to
#' \code{flexsurv::flexsurvreg()}, the parameters (in the order of flexsurv),
#' the parameters estimated on the log scale, the location parameter (which
#' carries the arm effect in dependent fits) and whether the distribution is
#' one of the seven standard parametric models of NICE DSU Technical Support
#' Document 14.
#'
#' @return A data frame with the columns \code{code}, \code{label},
#'   \code{flexsurv_dist}, \code{parameters}, \code{log_parameters} and
#'   \code{location} (comma-separated parameter names) and \code{standard}
#'   (logical).
#'
#' @keywords internal
#' @noRd
.distribution_info <- function() {
  data.frame(
    code = c("exp", "weibull", "weibullPH", "rayleigh", "lnorm", "llogis",
             "gompertz", "gamma", "gengamma", "gengamma.orig", "genf",
             "genf.orig"),
    label = c("Exponential", "Weibull", "Weibull (PH)", "Rayleigh",
              "Lognormal", "Log-logistic", "Gompertz", "Gamma",
              "Generalized Gamma", "Generalized Gamma (original)",
              "Generalized F", "Generalized F (original)"),
    flexsurv_dist = c("exp", "weibull", "weibullPH", "weibull", "lnorm",
                      "llogis", "gompertz", "gamma", "gengamma",
                      "gengamma.orig", "genf", "genf.orig"),
    parameters = c("rate", "shape,scale", "shape,scale", "scale",
                   "meanlog,sdlog", "shape,scale", "shape,rate",
                   "shape,rate", "mu,sigma,Q", "shape,scale,k",
                   "mu,sigma,Q,P", "mu,sigma,s1,s2"),
    log_parameters = c("rate", "shape,scale", "shape,scale", "scale",
                       "sdlog", "shape,scale", "rate", "shape,rate",
                       "sigma", "shape,scale,k", "sigma,P",
                       "sigma,s1,s2"),
    location = c("rate", "scale", "scale", "scale", "meanlog", "scale",
                 "rate", "rate", "mu", "scale", "mu", "mu"),
    standard = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
                 FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}
