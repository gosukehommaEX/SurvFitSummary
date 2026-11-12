#' Write the Survival Analysis Report to a Word Document
#'
#' @description
#' Writes a Word report with the structure of the survival analysis report
#' template of the health economics team (NICE TSD 14 models and
#' spline-based extrapolation), with print-ready names of the models:
#' \enumerate{
#'   \item Overview: data, outcome, arms, time unit, reporting interval and
#'     the variables of the hazard ratios. Items not given are left as
#'     placeholders in square brackets.
#'   \item Methods: Kaplan-Meier analysis, hazard ratios, diagnostics
#'     (with the table of the characteristics of the hazard functions),
#'     standard parametric models, spline-based extrapolation and software
#'     versions.
#'   \item Results: data summary, Kaplan-Meier curves (all arms and each
#'     arm) and tables, hazard ratios, smoothed empirical hazards,
#'     distribution diagnostics, proportional hazards assessment
#'     (Schoenfeld residuals, Grambsch-Therneau test, QQ plot), and for
#'     each arm and for both arms fitted together: the standard parametric
#'     models (overlaid plots, AIC and BIC, and hazard and survival plots
#'     of each model) and the spline models (placement of the knots, fit
#'     statistics, median, mean and predicted survival at the landmark
#'     times, and overlaid plots).
#' }
#' The figures are drawn as PNG files and embedded in the document; the
#' versions without and with confidence intervals are both included by
#' default, so that they do not have to be recomputed. The tables of fit
#' statistics and predicted survival of the spline models of both arms are
#' given for each arm, so that they fit the page.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param file Path of the \code{.docx} file to write (replaced if it
#'   exists).
#' @param data_label,outcome,author,version \code{NULL} or a string: the
#'   data (trial, data cut-off, population), the outcome (for example
#'   \code{"Overall survival"}), the author and the version of the report.
#' @param interval Spacing of the reporting times of the Kaplan-Meier
#'   tables. The default \code{NULL} uses the spacing of the time axis of
#'   the plots of the data (for example 12 months).
#' @param landmarks Times of the predicted survival of the spline models.
#'   The default \code{NULL} uses those of the report template: 6 months and
#'   1, 2, 5 and 10 years, in the time unit of the data.
#' @param time_horizon Time to which the models are extrapolated in the
#'   plots. The default \code{NULL} uses the last landmark time.
#' @param hazard_ratios Data frame returned by
#'   \code{\link{estimate_hazard_ratio}} (for example with the covariates
#'   and strata chosen for the analysis). The default \code{NULL} uses the
#'   Cox model of the arm without adjustment.
#' @param conf_int Logical vector: the versions of the figures without
#'   (\code{FALSE}) and with (\code{TRUE}) confidence intervals; both by
#'   default.
#' @param n_sim,seed Number of simulations and seed of the confidence
#'   intervals of the models (see \code{\link{predict_survival}}).
#' @param template \code{NULL} or the path of a \code{.docx} file whose
#'   page set-up and styles are used (for example an empty document with
#'   the styles of the organisation). The default \code{NULL} uses the
#'   template of the package, which follows the styles of the report
#'   template: A4 page, Arial 10 pt, numbered headings in dark red, bold
#'   captions, and tables with a dark red header row (made by
#'   \code{inst/templates/make_report_template.py}). The styles
#'   \code{"Title"}, \code{"Subtitle"}, \code{"heading 1"} to
#'   \code{"heading 4"}, \code{"Caption"}, \code{"Table Caption"},
#'   \code{"Figure"}, \code{"Table Text"} and the table style
#'   \code{"Report Table"} are used when the document defines them.
#' @param dpi Resolution of the figures. The default is 300.
#'
#' @return The path \code{file}, invisibly.
#'
#' @examples
#' \donttest{
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' fits <- fit_survival_models(dat, distributions = c("exp", "weibull"),
#'                             spline_scales = NULL,
#'                             fit_type = "independent")
#' path <- tempfile(fileext = ".docx")
#' write_survival_report(fits, path, outcome = "Overall survival",
#'                       conf_int = FALSE, dpi = 100)
#' unlink(path)
#' }
#'
#' @export
write_survival_report <- function(x,
                                  file,
                                  data_label = NULL,
                                  outcome = NULL,
                                  author = NULL,
                                  version = NULL,
                                  interval = NULL,
                                  landmarks = NULL,
                                  time_horizon = NULL,
                                  hazard_ratios = NULL,
                                  conf_int = c(FALSE, TRUE),
                                  n_sim = 1000,
                                  seed = 1,
                                  template = NULL,
                                  dpi = 300) {
  .check_survfitsummary(x)
  .check_output_file(file, "docx")
  overview <- list(data_label = data_label, outcome = outcome,
                   author = author, version = version,
                   report_date = Sys.Date())
  for (item in c("data_label", "outcome", "author", "version")) {
    value <- overview[[item]]
    if (!is.null(value) &&
        (!is.character(value) || length(value) != 1L || is.na(value))) {
      stop("'", item, "' must be NULL or a single string.", call. = FALSE)
    }
  }
  if (!is.logical(conf_int) || length(conf_int) == 0L ||
      anyNA(conf_int) || anyDuplicated(conf_int) > 0L) {
    stop("'conf_int' must be FALSE, TRUE or c(FALSE, TRUE).",
         call. = FALSE)
  }
  if (!is.null(template) &&
      (!is.character(template) || length(template) != 1L ||
       !grepl("\\.docx$", template, ignore.case = TRUE) ||
       !file.exists(template))) {
    stop("'template' must be NULL or the path of an existing .docx file.",
         call. = FALSE)
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
  if (!is.numeric(dpi) || length(dpi) != 1L || !is.finite(dpi) ||
      dpi <= 0) {
    stop("'dpi' must be a single positive number.", call. = FALSE)
  }
  s <- .report_settings(x, interval, landmarks, time_horizon,
                        hazard_ratios)
  blocks <- .report_blocks(x, s, overview, conf_int, n_sim, seed)
  .render_report(blocks, file, template = template, dpi = dpi)
}
