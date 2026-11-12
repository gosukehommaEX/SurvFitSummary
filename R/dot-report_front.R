#' Title and Overview of the Word Report
#'
#' @description
#' Internal helper: the title, the date and the Overview section of the
#' report template. Items not given by the user are left as placeholders in
#' square brackets, as in the template.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param s List returned by \code{.report_settings()}.
#' @param overview List with \code{data_label}, \code{outcome},
#'   \code{author}, \code{version} (each \code{NULL} or a string) and
#'   \code{report_date}.
#'
#' @return A list of elements made by \code{.report_block()}.
#'
#' @keywords internal
#' @noRd
.report_front <- function(x, s, overview) {
  b <- .report_block
  filled <- function(value, placeholder) {
    if (is.null(value)) placeholder else value
  }
  blocks <- list(
    b("title", "Survival Analysis Report"),
    b("subtitle", paste("Generated on:", format(overview$report_date)))
  )
  if (!is.null(overview$author)) {
    blocks <- c(blocks, list(b("paragraph",
                               paste("Prepared by:", overview$author))))
  }
  if (!is.null(overview$version)) {
    blocks <- c(blocks, list(b("paragraph",
                               paste("Version:", overview$version))))
  }
  arms_text <- paste(x$arms, collapse = ", ")
  if (length(x$arms) >= 2L) {
    arms_text <- paste0(arms_text, " (control arm: ", x$arms[1L], ")")
  }
  unit_title <- paste0(toupper(substr(s$unit, 1L, 1L)),
                       substring(s$unit, 2L))
  blocks <- c(blocks, list(
    b("heading", "Overview", level = 1L),
    b("paragraph", paste(
      "Data:", filled(overview$data_label, paste(
        "[user to complete: name of trial or of dataset, including data",
        "cut-off date and if ITT or a subgroup]"
      ))
    )),
    b("paragraph", paste(
      "Outcome:", filled(overview$outcome, paste(
        "[user to complete: PFS, OS, time to treatment discontinuation or",
        "death, or other]"
      ))
    )),
    b("paragraph", paste("Treatment arms:", arms_text)),
    b("paragraph", paste("Time unit:", unit_title)),
    b("paragraph", paste("Time intervals for reporting:",
                         format(s$interval), s$unit))
  ))
  hr <- s$hazard_ratios
  if (!is.null(hr)) {
    covariates <- setdiff(unique(hr$covariates), "None")
    strata <- setdiff(unique(hr$strata), "None")
    listed <- function(v) {
      if (length(v) == 0L) "none" else paste(v, collapse = "; ")
    }
    blocks <- c(blocks, list(b("paragraph", paste0(
      "Explanatory variables in the hazard ratio estimation: ",
      listed(covariates), "; stratification variables: ", listed(strata)
    ))))
  }
  blocks
}
