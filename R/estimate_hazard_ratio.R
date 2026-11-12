#' Hazard Ratios from Cox Proportional Hazards Models
#'
#' @description
#' Fits Cox proportional hazards models (\code{survival::coxph()}) with the
#' treatment arm as the explanatory variable, and reports the hazard ratio
#' of each arm against the control arm (the first level of \code{ARM}).
#' The unadjusted model is always fitted; models adjusted for
#' \code{covariates}, stratified by \code{strata}, or both are added when
#' these are given.
#'
#' @param data Data frame returned by \code{\link{processing_dataset}}, with
#'   two or more treatment arms.
#' @param covariates Character vector with the names of covariate columns of
#'   \code{data} for the adjusted models. The default is \code{NULL}.
#' @param strata Character vector with the names of stratification columns
#'   of \code{data} for the stratified models. The default is \code{NULL}.
#' @param ties Method for tied event times: \code{"efron"} (the default of
#'   \code{survival::coxph()}), \code{"breslow"} or \code{"exact"}.
#' @param conf_level Confidence level, a single number in (0, 1). The
#'   default is 0.95.
#'
#' @return A data frame with one row per model and comparison, and the
#'   columns \code{model} (\code{"Unadjusted, unstratified"},
#'   \code{"Unadjusted, stratified"}, \code{"Adjusted, unstratified"} or
#'   \code{"Adjusted, stratified"}), \code{covariates} and \code{strata} (the
#'   columns used, \code{"None"} if none), \code{comparison} (for example
#'   \code{"Treatment vs Control"}), \code{hr}, \code{lcl}, \code{ucl}
#'   (Wald confidence interval), \code{p_value} (Wald test), \code{n} and
#'   \code{events} (observations and events used by the model).
#'
#' @examples
#' raw <- generate_dummy_survival_data(
#'   arm = c("Control", "Treatment"), n = c(100, 100),
#'   hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 123
#' )
#' dat <- processing_dataset(raw, column_arm = "ARM",
#'                           column_survtime = "SURVTIME",
#'                           column_cnsr = "CNSR", control_arm = "Control")
#' estimate_hazard_ratio(dat)
#'
#' @importFrom survival Surv strata
#' @export
estimate_hazard_ratio <- function(data,
                                  covariates = NULL,
                                  strata = NULL,
                                  ties = c("efron", "breslow", "exact"),
                                  conf_level = 0.95) {
  .check_survival_data(data, min_arms = 2L)
  ties <- match.arg(ties)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  absent <- setdiff(c(covariates, strata), names(data))
  if (length(absent) > 0L) {
    stop("Column(s) not found in 'data': ", paste(absent, collapse = ", "),
         call. = FALSE)
  }
  data$ARM <- droplevels(data$ARM)
  arm_levels <- levels(data$ARM)
  quote_names <- function(x) paste0("`", x, "`")

  # Models to fit: (adjusted, stratified)
  specs <- list(c(FALSE, FALSE))
  if (!is.null(strata)) specs <- c(specs, list(c(FALSE, TRUE)))
  if (!is.null(covariates)) specs <- c(specs, list(c(TRUE, FALSE)))
  if (!is.null(covariates) && !is.null(strata)) {
    specs <- c(specs, list(c(TRUE, TRUE)))
  }

  rows <- lapply(specs, function(spec) {
    rhs <- "ARM"
    if (spec[1]) rhs <- c(rhs, quote_names(covariates))
    if (spec[2]) {
      rhs <- c(rhs, paste0("strata(",
                           paste(quote_names(strata), collapse = ", "), ")"))
    }
    model_formula <- stats::as.formula(
      paste("Surv(SURVTIME, EVENT) ~", paste(rhs, collapse = " + "))
    )
    fit <- survival::coxph(model_formula, data = data, ties = ties)
    fit_summary <- summary(fit, conf.int = conf_level)
    terms <- paste0("ARM", arm_levels[-1L])
    data.frame(
      model = paste0(if (spec[1]) "Adjusted" else "Unadjusted", ", ",
                     if (spec[2]) "stratified" else "unstratified"),
      covariates = if (spec[1]) paste(covariates, collapse = ", ") else "None",
      strata = if (spec[2]) paste(strata, collapse = ", ") else "None",
      comparison = paste(arm_levels[-1L], "vs", arm_levels[1L]),
      hr = unname(fit_summary$conf.int[terms, 1L]),
      lcl = unname(fit_summary$conf.int[terms, 3L]),
      ucl = unname(fit_summary$conf.int[terms, 4L]),
      p_value = unname(fit_summary$coefficients[terms, "Pr(>|z|)"]),
      n = as.integer(fit$n),
      events = as.integer(fit$nevent)
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
